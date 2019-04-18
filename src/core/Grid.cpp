// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Grid.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/ErrorHandler.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Partition.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/TextProcessing.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>

namespace ovk {

namespace grid_internal {

grid_base::grid_base(core::logger &Logger, core::error_handler &ErrorHandler, core::profiler
  &Profiler, int ID, const std::string &Name, int NumDims, MPI_Comm Comm):
  Exists_(true),
  Logger_(&Logger),
  ErrorHandler_(&ErrorHandler),
  Profiler_(&Profiler),
  ID_(ID),
  Name_(Name.length() > 0 ? Name : "Grid" + std::to_string(ID)),
  NumDims_(NumDims),
  Comm_(Comm)
{}

grid_base::~grid_base() {

  if (Exists_) {
    MPI_Barrier(Comm_);
    Logger_->LogStatus(Comm_.Rank() == 0, 0, "Destroyed grid %s.", Name_);
  }

}

}

grid::grid(int ID, const grid_params &Params, core::logger &Logger, core::error_handler
  &ErrorHandler, core::profiler &Profiler):
  grid_base(Logger, ErrorHandler, Profiler, ID, Params.Name_, Params.NumDims_, Params.Comm_),
  Cart_(NumDims_, {MakeUniformTuple<int>(0), Params.Size_}, Params.Periodic_,
    Params.PeriodicStorage_),
  PeriodicLength_(Params.PeriodicLength_),
  GeometryType_(Params.GeometryType_),
  PartitionHash_(NumDims_, Comm_, Cart_.Range(), Params.LocalRange_),
  Partition_(std::make_shared<core::partition>(Cart_, Comm_, Params.LocalRange_, 1, 1,
    core::DetectNeighbors(Cart_, Comm_, Params.LocalRange_, PartitionHash_), *Profiler_))
{

  if (Comm_.Rank() == 0) {
    PrintSummary_();
  }

  if (OVK_DEBUG) {
    PrintDecomposition_();
  }

  MPI_Barrier(Comm_);

}

grid::~grid() {

  if (Exists_) {
    MPI_Barrier(Comm_);
  }

}

void grid::PrintSummary_() const {

  std::string GlobalSizeString;
  for (int iDim = 0; iDim < NumDims_; ++iDim) {
    GlobalSizeString += core::FormatNumber(Cart_.Range().Size(iDim));
    if (iDim != NumDims_-1) GlobalSizeString += " x ";
  }

  std::string TotalPointsString = core::FormatNumber(Cart_.Range().Count(), "points", "point");
  std::string ProcessesString = core::FormatNumber(Comm_.Size(), "processes", "process");

  Logger_->LogStatus(true, 0, "Created grid %s (ID=%i): %s (%s) on %s.", Name_, ID_,
    GlobalSizeString, TotalPointsString, ProcessesString);

}

void grid::PrintDecomposition_() const {

  const range &LocalRange = Partition_->LocalRange();
  const array<core::partition_info> &Neighbors = Partition_->Neighbors();

  const char *DimNames[3] = {"i", "j", "k"};
  std::string LocalRangeString;
  for (int iDim = 0; iDim < NumDims_; ++iDim) {
    std::string LocalBeginString = core::FormatNumber(LocalRange.Begin(iDim));
    std::string LocalEndString = core::FormatNumber(LocalRange.End(iDim));
    LocalRangeString += std::string(DimNames[iDim]) + '=' + LocalBeginString + ':' + LocalEndString;
    if (iDim != NumDims_-1) LocalRangeString += ", ";
  }

  std::string TotalLocalPointsString = core::FormatNumber(LocalRange.Count(), "points", "point");

  std::string NeighborRanksString;
  for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
    // List separated by commas, so don't add thousands separators
    NeighborRanksString += std::to_string(Neighbors(iNeighbor).Rank);
    if (iNeighbor != Neighbors.Count()-1) NeighborRanksString += ", ";
  }

  Logger_->LogStatus(Comm_.Rank() == 0, 0, "Grid %s decomposition info:", Name_);

  MPI_Barrier(Comm_);

  for (int OtherRank = 0; OtherRank < Comm_.Size(); ++OtherRank) {
    if (OtherRank == Comm_.Rank()) {
      std::string RankString = core::FormatNumber(Comm_.Rank());
      Logger_->LogStatus(true, 1, "Rank %s (global rank @rank@) contains %s (%s).", RankString,
        LocalRangeString, TotalLocalPointsString);
      if (Neighbors.Count() > 0) {
        Logger_->LogStatus(true, 1, "Rank %s has neighbors: %s", RankString, NeighborRanksString);
      }
    }
    MPI_Barrier(Comm_);
  }

}

void CreateGridParams(grid_params &Params, int NumDims) {

  OVK_DEBUG_ASSERT(NumDims == 2 || NumDims == 3, "Invalid dimension.");

  Params.NumDims_ = NumDims;
  Params.Comm_ = MPI_COMM_NULL;

  Params.Size_ = MakeEmptyRange(NumDims).Size();

  Params.Periodic_ = MakeUniformTuple<bool>(false);
  Params.PeriodicStorage_ = periodic_storage::UNIQUE;
  Params.PeriodicLength_ = MakeUniformTuple<double>(0.);
  Params.GeometryType_ = geometry_type::CURVILINEAR;

  Params.LocalRange_ = MakeEmptyRange(NumDims);

}

void DestroyGridParams(grid_params &) {}

void GetGridParamName(const grid_params &Params, std::string &Name) {

  Name = Params.Name_;

}

void SetGridParamName(grid_params &Params, std::string Name) {

  Params.Name_ = std::move(Name);

}

void GetGridParamDimension(const grid_params &Params, int &NumDims) {

  NumDims = Params.NumDims_;

}

void GetGridParamComm(const grid_params &Params, MPI_Comm &Comm) {

  Comm = Params.Comm_;

}

void SetGridParamComm(grid_params &Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Params.Comm_ = Comm;

}

void GetGridParamSize(const grid_params &Params, int *Size) {

  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Size[iDim] = Params.Size_[iDim];
  }

}

void SetGridParamSize(grid_params &Params, const int *Size) {

  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  for (int iDim = 0; iDim < Params.NumDims_; ++iDim) {
    OVK_DEBUG_ASSERT(Size[iDim] > 0, "Size must be greater than 0 in each dimension.");
    Params.Size_[iDim] = Size[iDim];
  }

}

void GetGridParamPeriodic(const grid_params &Params, bool *Periodic) {

  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Periodic[iDim] = Params.Periodic_[iDim];
  }

}

void SetGridParamPeriodic(grid_params &Params, const bool *Periodic) {

  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  for (int iDim = 0; iDim < Params.NumDims_; ++iDim) {
    Params.Periodic_[iDim] = Periodic[iDim];
  }

}

void GetGridParamPeriodicStorage(const grid_params &Params, periodic_storage &PeriodicStorage) {

  PeriodicStorage = Params.PeriodicStorage_;

}

void SetGridParamPeriodicStorage(grid_params &Params, periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(ValidPeriodicStorage(PeriodicStorage), "Invalid periodic storage.");

  OVK_DEBUG_ASSERT(PeriodicStorage == periodic_storage::UNIQUE, "Duplicated periodic storage "
    "is not currently supported.");

  Params.PeriodicStorage_ = PeriodicStorage;

}

void GetGridParamPeriodicLength(const grid_params &Params, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = Params.PeriodicLength_[iDim];
  }

}

void SetGridParamPeriodicLength(grid_params &Params, const double *PeriodicLength) {

  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  for (int iDim = 0; iDim < Params.NumDims_; ++iDim) {
    OVK_DEBUG_ASSERT(PeriodicLength[iDim] >= 0., "Periodic length must be nonnegative.");
    Params.PeriodicLength_[iDim] = PeriodicLength[iDim];
  }

}

void GetGridParamGeometryType(const grid_params &Params, geometry_type &GeometryType) {

  GeometryType = Params.GeometryType_;

}

void SetGridParamGeometryType(grid_params &Params, geometry_type GeometryType) {

  OVK_DEBUG_ASSERT(ValidGeometryType(GeometryType), "Invalid geometry type.");

  Params.GeometryType_ = GeometryType;

}

void GetGridParamLocalRange(const grid_params &Params, range &LocalRange) {

  LocalRange = Params.LocalRange_;

}

void SetGridParamLocalRange(grid_params &Params, const range &LocalRange) {

  Params.LocalRange_ = LocalRange;

}

namespace core {

void CreateGridInfo(grid_info &Info, const grid *Grid, comm_view Comm) {

  bool IsLocal = Grid != nullptr;
  bool IsRoot = false;
  if (IsLocal) {
    IsRoot = Grid->CommRank() == 0;
  }

  int RootRank;
  if (IsRoot) RootRank = Comm.Rank();
  core::BroadcastAnySource(&RootRank, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) {
    Info.ID_ = Grid->ID();
  }
  MPI_Bcast(&Info.ID_, 1, MPI_INT, RootRank, Comm);

  int NameLength;
  if (IsRoot) NameLength = Grid->Name().length();
  MPI_Bcast(&NameLength, 1, MPI_INT, RootRank, Comm);
  array<char> NameChars({NameLength});
  if (IsRoot) NameChars.Fill(Grid->Name().begin());
  MPI_Bcast(NameChars.Data(), NameLength, MPI_CHAR, RootRank, Comm);
  Info.Name_.assign(NameChars.LinearBegin(), NameChars.LinearEnd());

  Info.RootRank_ = RootRank;

  int NumDims;
  if (IsRoot) {
    NumDims = Grid->Dimension();
  }
  MPI_Bcast(&NumDims, 1, MPI_INT, RootRank, Comm);

  Info.Cart_ = MakeEmptyCart(NumDims);

  if (IsRoot) {
    Info.Cart_.Range() = Grid->Cart().Range();
  }
  MPI_Bcast(Info.Cart_.Range().Begin().Data(), MAX_DIMS, MPI_INT, RootRank, Comm);
  MPI_Bcast(Info.Cart_.Range().End().Data(), MAX_DIMS, MPI_INT, RootRank, Comm);

  tuple<int> PeriodicInt;
  if (IsRoot) {
    PeriodicInt = tuple<int>(Grid->Cart().Periodic());
  }
  MPI_Bcast(PeriodicInt.Data(), MAX_DIMS, MPI_INT, RootRank, Comm);
  Info.Cart_.Periodic() = tuple<bool>(PeriodicInt);

  int PeriodicStorageInt;
  if (IsRoot) {
    PeriodicStorageInt = int(Grid->Cart().PeriodicStorage());
  }
  MPI_Bcast(&PeriodicStorageInt, 1, MPI_INT, RootRank, Comm);
  Info.Cart_.PeriodicStorage() = periodic_storage(PeriodicStorageInt);

  if (IsRoot) {
    Info.PeriodicLength_ = Grid->PeriodicLength();
  }
  MPI_Bcast(Info.PeriodicLength_.Data(), MAX_DIMS, MPI_DOUBLE, RootRank, Comm);

  int GeometryTypeInt;
  if (IsRoot) {
    GeometryTypeInt = int(Grid->GeometryType());
  }
  MPI_Bcast(&GeometryTypeInt, 1, MPI_INT, RootRank, Comm);
  Info.GeometryType_ = geometry_type(GeometryTypeInt);

  Info.IsLocal_ = IsLocal;

}

void DestroyGridInfo(grid_info &Info) {

  Info.Name_.clear();

}

}

void GetGridInfoID(const grid_info &Info, int &ID) {

  ID = Info.ID_;

}

void GetGridInfoName(const grid_info &Info, std::string &Name) {

  Name = Info.Name_;

}

void GetGridInfoRootRank(const grid_info &Info, int &RootRank) {

  RootRank = Info.RootRank_;

}

void GetGridInfoCart(const grid_info &Info, cart &Cart) {

  Cart = Info.Cart_;

}


void GetGridInfoPeriodicLength(const grid_info &Info, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = Info.PeriodicLength_[iDim];
  }

}

void GetGridInfoGeometryType(const grid_info &Info, geometry_type &GeometryType) {

  GeometryType = Info.GeometryType_;

}

void GetGridInfoGlobalRange(const grid_info &Info, range &GlobalRange) {

  GlobalRange = Info.Cart_.Range();

}

void GetGridInfoIsLocal(const grid_info &Info, bool &IsLocal) {

  IsLocal = Info.IsLocal_;

}

}
