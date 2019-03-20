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
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/TextProcessing.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <map>
#include <set>
#include <string>
#include <utility>

namespace ovk {

namespace {
void CreateNeighbors(grid &Grid);
void DestroyNeighbors(grid &Grid);
void PrintGridSummary(const grid &Grid);
void PrintGridDecomposition(const grid &Grid);
}

namespace core {

void CreateGrid(grid &Grid, int ID, const grid_params &Params, core::logger &Logger,
  core::error_handler &ErrorHandler) {

  int NumDims = Params.NumDims_;

  Grid.Comm_ = core::comm(Params.Comm_);

  MPI_Barrier(Grid.Comm_);

  Grid.Logger_ = &Logger;
  Grid.ErrorHandler_ = &ErrorHandler;

  Grid.ID_ = ID;

  if (Params.Name_.length() > 0) {
    Grid.Name_ = Params.Name_;
  } else {
    Grid.Name_ = "Grid" + std::to_string(ID);
  }

  Grid.Cart_ = cart(NumDims, {MakeUniformTuple<int>(0), Params.Size_}, Params.Periodic_,
    Params.PeriodicStorage_);

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Grid.PeriodicLength_[iDim] = Params.PeriodicLength_[iDim];
  }

  Grid.GeometryType_ = Params.GeometryType_;

  Grid.LocalRange_ = Params.LocalRange_;

  core::CreatePartitionHash(Grid.PartitionHash_, Grid.Cart_.Dimension(), Grid.Comm_,
    Grid.Cart_.Range(), Grid.LocalRange_);

  CreateNeighbors(Grid);

  if (Grid.Comm_.Rank() == 0) {
    PrintGridSummary(Grid);
  }

  if (OVK_DEBUG) {
    PrintGridDecomposition(Grid);
  }

  MPI_Barrier(Grid.Comm_);

}

void DestroyGrid(grid &Grid) {

  MPI_Barrier(Grid.Comm_);

  core::DestroyPartitionHash(Grid.PartitionHash_);

  DestroyNeighbors(Grid);

  MPI_Barrier(Grid.Comm_);

  LogStatus(*Grid.Logger_, Grid.Comm_.Rank() == 0, 0, "Destroyed grid %s.", Grid.Name_);

  Grid.Comm_.Reset();

}

}

void GetGridID(const grid &Grid, int &ID) {

  ID = Grid.ID_;

}

void GetGridName(const grid &Grid, std::string &Name) {

  Name = Grid.Name_;

}

void GetGridDimension(const grid &Grid, int &NumDims) {

  NumDims = Grid.Cart_.Dimension();

}

void GetGridComm(const grid &Grid, MPI_Comm &Comm) {

  Comm = Grid.Comm_.Get();

}

void GetGridCommSize(const grid &Grid, int &CommSize) {

  CommSize = Grid.Comm_.Size();

}

void GetGridCommRank(const grid &Grid, int &CommRank) {

  CommRank = Grid.Comm_.Rank();

}

void GetGridCart(const grid &Grid, cart &Cart) {

  Cart = Grid.Cart_;

}

void GetGridSize(const grid &Grid, int *Size) {

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Size[iDim] = Grid.Cart_.Range().Size(iDim);
  }

}

void GetGridPeriodic(const grid &Grid, bool *Periodic) {

  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Periodic[iDim] = Grid.Cart_.Periodic(iDim);
  }

}

void GetGridPeriodicStorage(const grid &Grid, periodic_storage &PeriodicStorage) {

  PeriodicStorage = Grid.Cart_.PeriodicStorage();

}

void GetGridPeriodicLength(const grid &Grid, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = Grid.PeriodicLength_[iDim];
  }

}

void GetGridGeometryType(const grid &Grid, geometry_type &GeometryType) {

  GeometryType = Grid.GeometryType_;

}

void GetGridGlobalRange(const grid &Grid, range &GlobalRange) {

  GlobalRange = Grid.Cart_.Range();

}

void GetGridLocalRange(const grid &Grid, range &LocalRange) {

  LocalRange = Grid.LocalRange_;

}

void GetGridGlobalCount(const grid &Grid, long long &NumGlobal) {

  NumGlobal = Grid.Cart_.Range().Count();

}

void GetGridLocalCount(const grid &Grid, long long &NumLocal) {

  NumLocal = Grid.LocalRange_.Count();

}

namespace {

void CreateNeighbors(grid &Grid) {

  int NumDims = Grid.Cart_.Dimension();

  const range &GlobalRange = Grid.Cart_.Range();
  const range &LocalRange = Grid.LocalRange_;
  const tuple<bool> &Periodic = Grid.Cart_.Periodic();

  tuple<bool> HasNeighborsBefore;
  tuple<bool> HasNeighborsAfter;
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    bool HasPartitionBefore = LocalRange.Begin(iDim) > GlobalRange.Begin(iDim);
    bool HasPartitionAfter = LocalRange.End(iDim) < GlobalRange.End(iDim);
    if (Periodic[iDim]) {
      HasNeighborsBefore[iDim] = HasPartitionBefore || HasPartitionAfter;
      HasNeighborsAfter[iDim] = HasPartitionBefore || HasPartitionAfter;
    } else {
      HasNeighborsBefore[iDim] = HasPartitionBefore;
      HasNeighborsAfter[iDim] = HasPartitionAfter;
    }
  }

  array<range> NeighborFaces;
  NeighborFaces.Reserve(2*MAX_DIMS);
//   static_array<range,2*MAX_DIMS> NeighborFaces;
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    if (HasNeighborsBefore[iDim]) {
      range &Face = NeighborFaces.Append(LocalRange);
      Face.Begin(iDim) -= 1;
      Face.End(iDim) = Face.Begin(iDim)+1;
      for (int jDim = iDim+1; jDim < NumDims; ++jDim) {
        Face.Begin(jDim) -= int(HasNeighborsBefore[jDim]);
        Face.End(jDim) += int(HasNeighborsAfter[jDim]);
      }
    }
    if (HasNeighborsAfter[iDim]) {
      range &Face = NeighborFaces.Append(LocalRange);
      Face.End(iDim) += 1;
      Face.Begin(iDim) = Face.End(iDim)-1;
      for (int jDim = iDim+1; jDim < NumDims; ++jDim) {
        Face.Begin(jDim) -= int(HasNeighborsBefore[jDim]);
        Face.End(jDim) += int(HasNeighborsAfter[jDim]);
      }
    }
  }

  int NumNeighborFaces = NeighborFaces.Count();

  long long NumNeighborPoints = 0;
  for (int iFace = 0; iFace < NumNeighborFaces; ++iFace) {
    NumNeighborPoints += NeighborFaces(iFace).Count();
  }

  array<int,2> NeighborPoints({{MAX_DIMS,NumNeighborPoints}});

  long long iNextPoint = 0;
  for (int iFace = 0; iFace < NumNeighborFaces; ++iFace) {
    range &Face = NeighborFaces(iFace);
    for (int k = Face.Begin(2); k < Face.End(2); ++k) {
      for (int j = Face.Begin(1); j < Face.End(1); ++j) {
        for (int i = Face.Begin(0); i < Face.End(0); ++i) {
          tuple<int> Point = Grid.Cart_.PeriodicAdjust({i,j,k});
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            NeighborPoints(iDim,iNextPoint) = Point[iDim];
          }
          ++iNextPoint;
        }
      }
    }
  }

  array<int> NeighborPointBinIndices({NumNeighborPoints});

  core::MapToPartitionBins(Grid.PartitionHash_, NeighborPoints, NeighborPointBinIndices);

  std::map<int, core::partition_bin> Bins;

  for (long long iPoint = 0; iPoint < NumNeighborPoints; ++iPoint) {
    Bins.emplace(NeighborPointBinIndices(iPoint), core::partition_bin());
  }

  core::RetrievePartitionBins(Grid.PartitionHash_, Bins);

  array<int> NeighborPointRanks({NumNeighborPoints});

  core::FindPartitions(Grid.PartitionHash_, Bins, NeighborPoints, NeighborPointBinIndices,
    NeighborPointRanks);

  Bins.clear();

  std::set<int> NeighborRanks;

  for (long long iPoint = 0; iPoint < NumNeighborPoints; ++iPoint) {
    NeighborRanks.insert(NeighborPointRanks(iPoint));
  }

  NeighborPoints.Clear();
  NeighborPointBinIndices.Clear();
  NeighborPointRanks.Clear();

  int NumNeighbors = NeighborRanks.size();

  Grid.Neighbors_.Resize({NumNeighbors});

  auto RankIter = NeighborRanks.begin();
  for (auto &Neighbor : Grid.Neighbors_) {
    Neighbor.Rank = *RankIter;
    ++RankIter;
  }

  array<int,2> LocalRangeValues({{2,MAX_DIMS}});
//   static_array<int,2*MAX_DIMS,2> LocalRangeValues({{2,MAX_DIMS}});
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    LocalRangeValues(0,iDim) = LocalRange.Begin(iDim);
    LocalRangeValues(1,iDim) = LocalRange.End(iDim);
  }

  array<int,3> NeighborRangeValues({{NumNeighbors,2,MAX_DIMS}});

  array<MPI_Request> Requests;
  Requests.Reserve(2*NumNeighbors);
  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    core::grid_neighbor &Neighbor = Grid.Neighbors_(iNeighbor);
    MPI_Request &RecvRequest = Requests.Append();
    MPI_Irecv(NeighborRangeValues.Data(iNeighbor,0,0), 2*MAX_DIMS, MPI_INT, Neighbor.Rank, 0,
      Grid.Comm_, &RecvRequest);
    MPI_Request &SendRequest = Requests.Append();
    MPI_Isend(LocalRangeValues.Data(), 2*MAX_DIMS, MPI_INT, Neighbor.Rank, 0, Grid.Comm_,
      &SendRequest);
  }
  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    core::grid_neighbor &Neighbor = Grid.Neighbors_(iNeighbor);
    Neighbor.LocalRange = range(NeighborRangeValues.Data(iNeighbor,0,0),
      NeighborRangeValues.Data(iNeighbor,1,0));
  }

}

void DestroyNeighbors(grid &Grid) {

  Grid.Neighbors_.Clear();

}

void PrintGridSummary(const grid &Grid) {

  int NumDims = Grid.Cart_.Dimension();

  std::string GlobalSizeString;
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    GlobalSizeString += core::FormatNumber(Grid.Cart_.Range().Size(iDim));
    if (iDim != NumDims-1) GlobalSizeString += " x ";
  }

  std::string TotalPointsString = core::FormatNumber(Grid.Cart_.Range().Count(), "points", "point");
  std::string ProcessesString = core::FormatNumber(Grid.Comm_.Size(), "processes", "process");

  LogStatus(*Grid.Logger_, true, 0, "Created grid %s (ID=%i): %s (%s) on %s.", Grid.Name_, Grid.ID_,
    GlobalSizeString, TotalPointsString, ProcessesString);

}

void PrintGridDecomposition(const grid &Grid) {

  int NumDims = Grid.Cart_.Dimension();

  const char *DimNames[3] = {"i", "j", "k"};
  std::string LocalRangeString;
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    std::string LocalBeginString = core::FormatNumber(Grid.LocalRange_.Begin(iDim));
    std::string LocalEndString = core::FormatNumber(Grid.LocalRange_.End(iDim));
    LocalRangeString += std::string(DimNames[iDim]) + '=' + LocalBeginString + ':' + LocalEndString;
    if (iDim != NumDims-1) LocalRangeString += ", ";
  }

  std::string TotalLocalPointsString = core::FormatNumber(Grid.LocalRange_.Count(), "points",
    "point");

  std::string NeighborRanksString;
  for (int iNeighbor = 0; iNeighbor < int(Grid.Neighbors_.Count()); ++iNeighbor) {
    // List separated by commas, so don't add thousands separators
    NeighborRanksString += std::to_string(Grid.Neighbors_(iNeighbor).Rank);
    if (iNeighbor != int(Grid.Neighbors_.Count())-1) NeighborRanksString += ", ";
  }

  LogStatus(*Grid.Logger_, Grid.Comm_.Rank() == 0, 0, "Grid %s decomposition info:", Grid.Name_);

  MPI_Barrier(Grid.Comm_);

  for (int OtherRank = 0; OtherRank < Grid.Comm_.Size(); ++OtherRank) {
    if (OtherRank == Grid.Comm_.Rank()) {
      std::string RankString = core::FormatNumber(Grid.Comm_.Rank());
      LogStatus(*Grid.Logger_, true, 1, "Rank %s (global rank @rank@) contains %s (%s).",
        RankString, LocalRangeString, TotalLocalPointsString);
      if (Grid.Neighbors_.Count() > 0) {
        LogStatus(*Grid.Logger_, true, 1, "Rank %s has neighbors: %s", RankString,
          NeighborRanksString);
      }
    }
    MPI_Barrier(Grid.Comm_);
  }

}

}

namespace core {

comm_view GetGridComm(const grid &Grid) {

  return Grid.Comm_;

}

const array<grid_neighbor> &GetGridNeighbors(const grid &Grid) {

  return Grid.Neighbors_;

}

const partition_hash &GetGridPartitionHash(const grid &Grid) {

  return Grid.PartitionHash_;

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
    IsRoot = Grid->Comm_.Rank() == 0;
  }

  int RootRank;
  if (IsRoot) RootRank = Comm.Rank();
  core::BroadcastAnySource(&RootRank, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) {
    Info.ID_ = Grid->ID_;
  }
  MPI_Bcast(&Info.ID_, 1, MPI_INT, RootRank, Comm);

  int NameLength;
  if (IsRoot) NameLength = Grid->Name_.length();
  MPI_Bcast(&NameLength, 1, MPI_INT, RootRank, Comm);
  array<char> NameChars({NameLength});
  if (IsRoot) NameChars.Fill(Grid->Name_.begin());
  MPI_Bcast(NameChars.Data(), NameLength, MPI_CHAR, RootRank, Comm);
  Info.Name_.assign(NameChars.LinearBegin(), NameChars.LinearEnd());

  Info.RootRank_ = RootRank;

  int NumDims;
  if (IsRoot) {
    NumDims = Grid->Cart_.Dimension();
  }
  MPI_Bcast(&NumDims, 1, MPI_INT, RootRank, Comm);

  Info.Cart_ = MakeEmptyCart(NumDims);

  if (IsRoot) {
    Info.Cart_.Range() = Grid->Cart_.Range();
  }
  MPI_Bcast(Info.Cart_.Range().Begin().Data(), MAX_DIMS, MPI_INT, RootRank, Comm);
  MPI_Bcast(Info.Cart_.Range().End().Data(), MAX_DIMS, MPI_INT, RootRank, Comm);

  tuple<int> PeriodicInt;
  if (IsRoot) {
    PeriodicInt = tuple<int>(Grid->Cart_.Periodic());
  }
  MPI_Bcast(PeriodicInt.Data(), MAX_DIMS, MPI_INT, RootRank, Comm);
  Info.Cart_.Periodic() = tuple<bool>(PeriodicInt);

  int PeriodicStorageInt;
  if (IsRoot) {
    PeriodicStorageInt = int(Grid->Cart_.PeriodicStorage());
  }
  MPI_Bcast(&PeriodicStorageInt, 1, MPI_INT, RootRank, Comm);
  Info.Cart_.PeriodicStorage() = periodic_storage(PeriodicStorageInt);

  if (IsRoot) {
    Info.PeriodicLength_ = Grid->PeriodicLength_;
  }
  MPI_Bcast(Info.PeriodicLength_.Data(), MAX_DIMS, MPI_DOUBLE, RootRank, Comm);

  int GeometryTypeInt;
  if (IsRoot) {
    GeometryTypeInt = int(Grid->GeometryType_);
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
