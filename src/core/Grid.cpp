// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Grid.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/CommunicationOps.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Decomp.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Partition.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/TextProcessing.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <memory>
#include <string>
#include <utility>

namespace ovk {

namespace grid_internal {

grid_base::grid_base(std::shared_ptr<context> &&Context, std::string &&Name, comm &&Comm):
  Context_(std::move(Context)),
  Name_(std::move(Name)),
  Comm_(std::move(Comm))
{} 

grid_base::~grid_base() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Logger.LogStatus(Comm_.Rank() == 0, "Destroyed grid %s.", *Name_);
  }

}

}

grid::grid(std::shared_ptr<context> &&Context, params &&Params):
  grid(std::move(Context), std::move(Params), Params.NumDims_, DuplicateComm(Params.Comm_),
    Params.Cart_, Params.LocalRange_)
{}

grid::grid(std::shared_ptr<context> &&Context, params &&Params, int NumDims, comm &&Comm, const cart
  &Cart, const range &LocalRange):
  grid(std::move(Context), std::move(Params), NumDims, std::move(Comm), Cart, LocalRange,
    core::CartPointToCell(Cart), core::RangePointToCell(Cart, LocalRange))
{}

grid::grid(std::shared_ptr<context> &&Context, params &&Params, int NumDims, comm &&Comm, const cart
  &Cart, const range &LocalRange, const cart &CellCart, const range &CellLocalRange):
  grid(std::move(Context), std::move(Params), NumDims, std::move(Comm), Cart, LocalRange, CellCart,
    CellLocalRange, core::ExtendLocalRange(CellCart, CellLocalRange, 2))
{}

grid::grid(std::shared_ptr<context> &&Context, params &&Params, int NumDims, comm &&Comm, const cart
  &Cart, const range &LocalRange, const cart &CellCart, const range &CellLocalRange, const range
  &CellExtendedRange):
  grid(std::move(Context), std::move(Params), NumDims, std::move(Comm), Cart, LocalRange,
    core::RangeCellToPointAll(Cart, CellExtendedRange), CellCart, CellLocalRange, CellExtendedRange)
{}

grid::grid(std::shared_ptr<context> &&Context, params &&Params, int NumDims, comm &&Comm, const cart
  &Cart, const range &LocalRange, const range &ExtendedRange, const cart &CellCart, const range
  &CellLocalRange, const range &CellExtendedRange):
  grid(std::move(Context), std::move(Params), NumDims, std::move(Comm), Cart, LocalRange,
    ExtendedRange, CellCart, CellLocalRange, CellExtendedRange, core::DetectNeighbors(Cart, Comm,
    LocalRange, core::CreateDecompHash(NumDims, Comm, LocalRange)))
{}

grid::grid(std::shared_ptr<context> &&Context, params &&Params, int NumDims, comm &&Comm, const cart
  &Cart, const range &LocalRange, const range &ExtendedRange, const cart &CellCart, const range
  &CellLocalRange, const range &CellExtendedRange, const array<int> &NeighborRanks):
  grid_base(std::move(Context), std::move(*Params.Name_), std::move(Comm)),
  NumDims_(NumDims),
  Partition_(std::make_shared<partition>(Context_, Cart, Comm_, LocalRange, ExtendedRange, 1,
    NeighborRanks)),
  CellPartition_(std::make_shared<partition>(Context_, CellCart, Comm_, CellLocalRange,
    CellExtendedRange, 1, NeighborRanks))
{

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {

    if (Comm_.Rank() == 0) {
      std::string ProcessesString = core::FormatNumber(Comm_.Size(), "processes", "process");
      Logger.LogStatus(true, "Created grid %s on %s.", *Name_, ProcessesString);
    }

    Logger.LogStatus(Comm_.Rank() == 0, "Grid %s info:", *Name_);
    auto Indent1 = Logger.IndentStatus();

    if (Comm_.Rank() == 0) {
      const cart &Cart = Partition_->Cart();
      std::string GlobalSizeString;
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        GlobalSizeString += core::FormatNumber(Cart.Range().Size(iDim));
        if (iDim != NumDims_-1) GlobalSizeString += " x ";
      }
      std::string TotalPointsString = core::FormatNumber(Cart.Range().Count(), "points", "point");
      Logger.LogStatus(true, "Size: %s (%s)", GlobalSizeString, TotalPointsString);
    }

    auto Level1 = Logger.IncreaseStatusLevel(10);

    if (Logger.LoggingStatus()) {

      const range &LocalRange = Partition_->LocalRange();
      const map<int,partition::neighbor_info> &Neighbors = Partition_->Neighbors();

      const char *DimNames[3] = {"i", "j", "k"};
      std::string LocalRangeString;
      for (int iDim = 0; iDim < NumDims_; ++iDim) {
        std::string LocalBeginString = core::FormatNumber(LocalRange.Begin(iDim));
        std::string LocalEndString = core::FormatNumber(LocalRange.End(iDim));
        LocalRangeString += std::string(DimNames[iDim]) + '=' + LocalBeginString + ':' +
          LocalEndString;
        if (iDim != NumDims_-1) LocalRangeString += ", ";
      }

      std::string TotalLocalPointsString = core::FormatNumber(LocalRange.Count(), "points",
        "point");

      std::string NeighborRanksString;
      for (int iNeighbor = 0; iNeighbor < Neighbors.Count(); ++iNeighbor) {
        // List separated by commas, so don't add thousands separators
        NeighborRanksString += std::to_string(Neighbors[iNeighbor].Key());
        if (iNeighbor != Neighbors.Count()-1) NeighborRanksString += ", ";
      }

      Logger.LogStatus(Comm_.Rank() == 0, "Decomposition:");
      auto Indent2 = Logger.IndentStatus();

      MPI_Barrier(Comm_);

      for (int OtherRank = 0; OtherRank < Comm_.Size(); ++OtherRank) {
        if (OtherRank == Comm_.Rank()) {
          std::string RankString = core::FormatNumber(Comm_.Rank());
          Logger.LogStatus(true, "Rank %s (global rank @rank@) contains %s (%s).", RankString,
            LocalRangeString, TotalLocalPointsString);
          if (Neighbors.Count() > 0) {
            Logger.LogStatus(true, "Rank %s has neighbors: %s", RankString, NeighborRanksString);
          }
        }
        MPI_Barrier(Comm_);
      }

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

  range GlobalRange = MakeEmptyRange(NumDims);
  tuple<bool> Periodic = MakeUniformTuple<bool>(NumDims, false);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    GlobalRange.Begin(iDim) = Cart_.Range().Begin(iDim);
    GlobalRange.End(iDim) = Cart_.Range().End(iDim);
    Periodic(iDim) = Cart_.Periodic(iDim);
  }
  Cart_ = cart(NumDims, GlobalRange, Periodic, Cart_.PeriodicStorage());

  for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
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

grid::params &grid::params::SetCart(const cart &Cart) {

  OVK_DEBUG_ASSERT(Cart.Dimension() == NumDims_, "Cart has incorrect dimension.");
  if (OVK_DEBUG) {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      OVK_DEBUG_ASSERT(Cart.Range().Begin(iDim) == 0, "Cart range begin must be 0.");
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      OVK_DEBUG_ASSERT(Cart.Range().Size(iDim) > 0, "Size must be greater than 0 in each "
        "dimension.");
    }
  }
  OVK_DEBUG_ASSERT(Cart.PeriodicStorage() == periodic_storage::UNIQUE, "Duplicated periodic "
    "storage is not currently supported.");

  Cart_ = Cart;

  return *this;

}

grid::params &grid::params::SetGlobalRange(const range &GlobalRange) {

  OVK_DEBUG_ASSERT(!GlobalRange.Empty(), "Global range cannot be empty.");
  if (OVK_DEBUG) {
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      OVK_DEBUG_ASSERT(GlobalRange.Begin(iDim) == 0, "Global range has incorrect dimension.");
      OVK_DEBUG_ASSERT(GlobalRange.End(iDim) == 1, "Global range has incorrect dimension.");
    }
  }

  Cart_.Range() = GlobalRange;

  return *this;

}

grid::params &grid::params::SetLocalRange(const range &LocalRange) {

  OVK_DEBUG_ASSERT(!LocalRange.Empty(), "Local range cannot be empty.");
  if (OVK_DEBUG) {
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      OVK_DEBUG_ASSERT(LocalRange.Begin(iDim) == 0, "Local range has incorrect dimension.");
      OVK_DEBUG_ASSERT(LocalRange.End(iDim) == 1, "Local range has incorrect dimension.");
    }
  }

  LocalRange_ = LocalRange;

  return *this;

}

grid::params &grid::params::SetPeriodic(const tuple<bool> &Periodic) {

  if (OVK_DEBUG) {
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      OVK_DEBUG_ASSERT(Periodic(iDim) == false, "Periodic has incorrect dimension.");
    }
  }

  Cart_.Periodic() = Periodic;

  return *this;

}

grid::params &grid::params::SetPeriodicStorage(periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(PeriodicStorage == periodic_storage::UNIQUE, "Duplicated periodic storage "
    "is not currently supported.");

  Cart_.PeriodicStorage() = PeriodicStorage;

  return *this;

}

grid_info::grid_info(grid *MaybeGrid, comm_view Comm) {

  IsLocal_ = MaybeGrid != nullptr;
  bool IsRoot = false;
  if (IsLocal_) {
    IsRoot = MaybeGrid->Comm().Rank() == 0;
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

  if (IsRoot) {
    Cart_.Periodic() = MaybeGrid->Cart().Periodic();
  }
  MPI_Bcast(Cart_.Periodic().Data(), MAX_DIMS, MPI_C_BOOL, RootRank_, Comm);

  if (IsRoot) {
    Cart_.PeriodicStorage() = MaybeGrid->Cart().PeriodicStorage();
  }
  MPI_Bcast(&Cart_.PeriodicStorage(), 1, core::GetMPIDataType<periodic_storage>(), RootRank_, Comm);

  CellCart_ = core::CartPointToCell(Cart_);

}

grid_info grid_info::internal_Create(grid *MaybeGrid, comm_view Comm) {

  return {MaybeGrid, Comm};

}

namespace core {

grid_info CreateGridInfo(grid *MaybeGrid, comm_view Comm) {

  return grid_info::internal_Create(MaybeGrid, Comm);

}

}

}
