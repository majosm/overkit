// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Grid.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Partition.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/TextProcessing.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <memory>
#include <set>
#include <string>
#include <utility>

namespace ovk {

namespace grid_internal {

grid_base::grid_base(std::shared_ptr<context> &&Context, std::string &&Name, MPI_Comm Comm):
  Context_(std::move(Context)),
  Name_(std::move(Name)),
  Comm_(Comm)
{} 

grid_base::~grid_base() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Logger.LogDebug(Comm_.Rank() == 0, 0, "Destroyed grid %s.", *Name_);
  }

}

}

grid::grid(std::shared_ptr<context> &&Context, params &&Params):
  grid_base(std::move(Context), std::move(*Params.Name_), Params.Comm_),
  FloatingRefGenerator_(*this),
  NumDims_(Params.NumDims_),
  Cart_(NumDims_, {Params.Size_}, Params.Periodic_, Params.PeriodicStorage_),
  PeriodicLength_(Params.PeriodicLength_),
  GeometryType_(Params.GeometryType_),
  PartitionHash_(NumDims_, Comm_, Cart_.Range(), Params.LocalRange_),
  Partition_(std::make_shared<core::partition>(Context_, Cart_, Comm_, Params.LocalRange_, 1,
    1, core::DetectNeighbors(Cart_, Comm_, Params.LocalRange_, PartitionHash_)))
{

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingDebug()) {

    if (Comm_.Rank() == 0) {
      std::string ProcessesString = core::FormatNumber(Comm_.Size(), "processes", "process");
      Logger.LogDebug(true, 0, "Created grid %s on %s.", *Name_, ProcessesString);
    }

    Logger.LogDebug(Comm_.Rank() == 0, 0, "Grid %s info:", *Name_);

    if (Comm_.Rank() == 0) {
      std::string GlobalSizeString;
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        GlobalSizeString += core::FormatNumber(Cart_.Range().Size(iDim));
        if (iDim != NumDims_-1) GlobalSizeString += " x ";
      }
      std::string TotalPointsString = core::FormatNumber(Cart_.Range().Count(), "points", "point");
      Logger.LogDebug(true, 1, "Size: %s (%s)", GlobalSizeString, TotalPointsString);
    }

    const range &LocalRange = Partition_->LocalRange();
    const array<core::partition_info> &Neighbors = Partition_->Neighbors();

    const char *DimNames[3] = {"i", "j", "k"};
    std::string LocalRangeString;
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      std::string LocalBeginString = core::FormatNumber(LocalRange.Begin(iDim));
      std::string LocalEndString = core::FormatNumber(LocalRange.End(iDim));
      LocalRangeString += std::string(DimNames[iDim]) + '=' + LocalBeginString + ':' +
        LocalEndString;
      if (iDim != NumDims_-1) LocalRangeString += ", ";
    }

    std::string TotalLocalPointsString = core::FormatNumber(LocalRange.Count(), "points", "point");

    std::string NeighborRanksString;
    for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
      // List separated by commas, so don't add thousands separators
      NeighborRanksString += std::to_string(Neighbors(iNeighbor).Rank);
      if (iNeighbor != Neighbors.Count()-1) NeighborRanksString += ", ";
    }

    Logger.LogDebug(Comm_.Rank() == 0, 1, "Decomposition:");

    MPI_Barrier(Comm_);

    for (int OtherRank = 0; OtherRank < Comm_.Size(); ++OtherRank) {
      if (OtherRank == Comm_.Rank()) {
        std::string RankString = core::FormatNumber(Comm_.Rank());
        Logger.LogDebug(true, 2, "Rank %s (global rank @rank@) contains %s (%s).", RankString,
          LocalRangeString, TotalLocalPointsString);
        if (Neighbors.Count() > 0) {
          Logger.LogDebug(true, 2, "Rank %s has neighbors: %s", RankString, NeighborRanksString);
        }
      }
      MPI_Barrier(Comm_);
    }

  }

}

grid::~grid() noexcept {

  if (Context_) {
    // Barrier before cleaning up
    MPI_Barrier(Comm_);
  }

}

grid grid::internal_Create(std::shared_ptr<context> &&Context, params &&Params) {

  return {std::move(Context), std::move(Params)};

}

namespace core {

grid CreateGrid(std::shared_ptr<context> Context, grid::params Params) {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return grid::internal_Create(std::move(Context), std::move(Params));

}

}

grid::params &grid::params::SetName(std::string Name) {

  Name_ = std::move(Name);

  return *this;

}

grid::params &grid::params::SetDimension(int NumDims) {

  OVK_DEBUG_ASSERT(NumDims == 2 || NumDims == 3, "Invalid dimension.");

  NumDims_ = NumDims;

  for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
    Size_(iDim) = 1;
    Periodic_(iDim) = false;
    PeriodicLength_(iDim) = 0.;
    LocalRange_.Begin(iDim) = 0;
    LocalRange_.End(iDim) = 1;
  }

  return *this;

}

grid::params &grid::params::SetComm(MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Comm_ = Comm;

  return *this;

}

grid::params &grid::params::SetSize(const tuple<int> &Size) {

  if (OVK_DEBUG) {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      OVK_DEBUG_ASSERT(Size(iDim) > 0, "Size must be greater than 0 in each dimension.");
    }
  }

  Size_ = Size;

  for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
    Size_(iDim) = 1;
  }

  return *this;

}

grid::params &grid::params::SetPeriodic(const tuple<bool> &Periodic) {

  Periodic_ = Periodic;

  for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
    Periodic_(iDim) = false;
  }

  return *this;

}

grid::params &grid::params::SetPeriodicStorage(periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(PeriodicStorage == periodic_storage::UNIQUE, "Duplicated periodic storage "
    "is not currently supported.");

  PeriodicStorage_ = PeriodicStorage;

  return *this;

}

grid::params &grid::params::SetPeriodicLength(const tuple<double> &PeriodicLength) {

  if (OVK_DEBUG) {
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      OVK_DEBUG_ASSERT(PeriodicLength(iDim) >= 0., "Periodic length must be nonnegative.");
    }
  }

  PeriodicLength_ = PeriodicLength;

  for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength_(iDim) = 0.;
  }

  return *this;

}

grid::params &grid::params::SetGeometryType(geometry_type GeometryType) {

  OVK_DEBUG_ASSERT(ValidGeometryType(GeometryType), "Invalid geometry type.");

  GeometryType_ = GeometryType;

  return *this;

}

grid::params &grid::params::SetLocalRange(const range &LocalRange) {

  LocalRange_ = LocalRange;

  for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
    LocalRange_.Begin(iDim) = 0;
    LocalRange_.End(iDim) = 1;
  }

  return *this;

}

grid_info::grid_info(grid *MaybeGrid, core::comm_view Comm) {

  IsLocal_ = MaybeGrid != nullptr;
  bool IsRoot = false;
  if (IsLocal_) {
    IsRoot = MaybeGrid->CommRank() == 0;
  }

  if (IsRoot) RootRank_ = Comm.Rank();
  core::BroadcastAnySource(&RootRank_, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) Name_ = MaybeGrid->Name();
  core::BroadcastString(*Name_, RootRank_, Comm);

  int NumDims;
  if (IsRoot) {
    NumDims = MaybeGrid->Dimension();
  }
  MPI_Bcast(&NumDims, 1, MPI_INT, RootRank_, Comm);

  Cart_ = MakeEmptyCart(NumDims);

  if (IsRoot) {
    Cart_.Range() = MaybeGrid->Cart().Range();
  }
  MPI_Bcast(Cart_.Range().Begin().Data(), MAX_DIMS, MPI_INT, RootRank_, Comm);
  MPI_Bcast(Cart_.Range().End().Data(), MAX_DIMS, MPI_INT, RootRank_, Comm);

  tuple<int> PeriodicInt;
  if (IsRoot) {
    PeriodicInt = tuple<int>(MaybeGrid->Cart().Periodic());
  }
  MPI_Bcast(PeriodicInt.Data(), MAX_DIMS, MPI_INT, RootRank_, Comm);
  Cart_.Periodic() = tuple<bool>(PeriodicInt);

  int PeriodicStorageInt;
  if (IsRoot) {
    PeriodicStorageInt = int(MaybeGrid->Cart().PeriodicStorage());
  }
  MPI_Bcast(&PeriodicStorageInt, 1, MPI_INT, RootRank_, Comm);
  Cart_.PeriodicStorage() = periodic_storage(PeriodicStorageInt);

  if (IsRoot) {
    PeriodicLength_ = MaybeGrid->PeriodicLength();
  }
  MPI_Bcast(PeriodicLength_.Data(), MAX_DIMS, MPI_DOUBLE, RootRank_, Comm);

  int GeometryTypeInt;
  if (IsRoot) {
    GeometryTypeInt = int(MaybeGrid->GeometryType());
  }
  MPI_Bcast(&GeometryTypeInt, 1, MPI_INT, RootRank_, Comm);
  GeometryType_ = geometry_type(GeometryTypeInt);

}

grid_info grid_info::internal_Create(grid *MaybeGrid, core::comm_view Comm) {

  return {MaybeGrid, Comm};

}

namespace core {

grid_info CreateGridInfo(grid *MaybeGrid, comm_view Comm) {

  return grid_info::internal_Create(MaybeGrid, Comm);

}

}

}
