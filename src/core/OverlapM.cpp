// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/OverlapM.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Event.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Logger.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

namespace overlap_m_internal {

overlap_m_base::overlap_m_base(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
  &&DestinationGridInfo):
  Context_(std::move(Context)),
  Grid_(&Grid),
  DestinationGridInfo_(std::move(DestinationGridInfo)),
  Comm_(Grid_->Comm())
{
  MPI_Barrier(Comm_);
}

overlap_m_base::~overlap_m_base() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Logger.LogDebug(Comm_.Rank() == 0, 0, "Destroyed M side of overlap (%s,%s).",
      Grid_->Name(), DestinationGridInfo_.Name());
  }

}

}

overlap_m::overlap_m(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
  &&DestinationGridInfo):
  overlap_m_base(std::move(Context), Grid, std::move(DestinationGridInfo)),
  NumDims_(Grid_->Dimension()),
  Count_(0),
  Cells_({{MAX_DIMS,0}}),
  Coords_({{MAX_DIMS,0}}),
  Destinations_({{MAX_DIMS,0}}),
  DestinationRanks_({0})
{

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();
  Logger.LogDebug(Comm_.Rank() == 0, 0, "Created M side of overlap (%s,%s).", Grid_->Name(),
    DestinationGridInfo_.Name());

}

overlap_m::~overlap_m() noexcept {

  if (Context_) {
    // Barrier before cleaning up
    MPI_Barrier(Comm_);
  }

}

overlap_m overlap_m::internal_Create(std::shared_ptr<context> &&Context, const grid &Grid,
  grid_info &&DestinationGridInfo) {

  return {std::move(Context), Grid, std::move(DestinationGridInfo)};

}

namespace core {

overlap_m CreateOverlapM(std::shared_ptr<context> Context, const grid &Grid, grid_info
  DestinationGridInfo) {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return overlap_m::internal_Create(std::move(Context), Grid, std::move(DestinationGridInfo));

}

}

void overlap_m::Resize(long long Count) {

  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");

  MPI_Barrier(Comm_);

  OVK_DEBUG_ASSERT(!CellsEditor_.Active(), "Cannot resize while editing cells.");
  OVK_DEBUG_ASSERT(!CoordsEditor_.Active(), "Cannot resize while editing coords.");
  OVK_DEBUG_ASSERT(!DestinationsEditor_.Active(), "Cannot resize while editing destinations.");
  OVK_DEBUG_ASSERT(!DestinationRanksEditor_.Active(), "Cannot resize while editing destination "
    "ranks.");

  Count_ = Count;

  Cells_.Resize({{MAX_DIMS,Count}});
  Coords_.Resize({{MAX_DIMS,Count}});
  Destinations_.Resize({{MAX_DIMS,Count}});
  DestinationRanks_.Resize({Count});

  for (long long iCell = 0; iCell < Count_; ++iCell) {
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      Cells_(iDim,iCell) = Grid_->CellGlobalRange().Begin(iDim)-1;
    }
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      Cells_(iDim,iCell) = 0;
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Coords_(iDim,iCell) = 0.;
    }
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      Destinations_(iDim,iCell) = DestinationGridInfo_.GlobalRange().Begin(iDim)-1;
    }
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      Destinations_(iDim,iCell) = 0;
    }
    DestinationRanks_(iCell) = -1;
  }

  MPI_Barrier(Comm_);

  ResizeEvent_.Trigger();
  CellsEvent_.Trigger();
  CoordsEvent_.Trigger();
  DestinationsEvent_.Trigger();
  DestinationRanksEvent_.Trigger();

  MPI_Barrier(Comm_);

}

bool overlap_m::EditingCells() const {

  return CellsEditor_.Active();

}

edit_handle<array<int,2>> overlap_m::EditCells() {

  if (!CellsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<overlap_m> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      overlap_m &OverlapM = *FloatingRef;
      MPI_Barrier(OverlapM.Comm_);
      OverlapM.CellsEvent_.Trigger();
      MPI_Barrier(OverlapM.Comm_);
    };
    CellsEditor_.Activate(std::move(DeactivateFunc));
  }

  return CellsEditor_.Edit(Cells_);

}

void overlap_m::RestoreCells() {

  OVK_DEBUG_ASSERT(CellsEditor_.Active(), "Unable to restore cells; not currently being "
    "edited.");

  CellsEditor_.Restore();

}

bool overlap_m::EditingCoords() const {

  return CoordsEditor_.Active();

}

edit_handle<array<double,2>> overlap_m::EditCoords() {

  if (!CoordsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<overlap_m> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      overlap_m &OverlapM = *FloatingRef;
      MPI_Barrier(OverlapM.Comm_);
      OverlapM.CoordsEvent_.Trigger();
      MPI_Barrier(OverlapM.Comm_);
    };
    CoordsEditor_.Activate(std::move(DeactivateFunc));
  }

  return CoordsEditor_.Edit(Coords_);

}

void overlap_m::RestoreCoords() {

  OVK_DEBUG_ASSERT(CoordsEditor_.Active(), "Unable to restore coords; not currently being edited.");

  CoordsEditor_.Restore();

}

bool overlap_m::EditingDestinations() const {

  return DestinationsEditor_.Active();

}

edit_handle<array<int,2>> overlap_m::EditDestinations() {

  if (!DestinationsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<overlap_m> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      overlap_m &OverlapM = *FloatingRef;
      MPI_Barrier(OverlapM.Comm_);
      OverlapM.DestinationsEvent_.Trigger();
      MPI_Barrier(OverlapM.Comm_);
    };
    DestinationsEditor_.Activate(std::move(DeactivateFunc));
  }

  return DestinationsEditor_.Edit(Destinations_);

}

void overlap_m::RestoreDestinations() {

  OVK_DEBUG_ASSERT(DestinationsEditor_.Active(), "Unable to restore destinations; not currently "
    "being edited.");

  DestinationsEditor_.Restore();

}

bool overlap_m::EditingDestinationRanks() const {

  return DestinationRanksEditor_.Active();

}

edit_handle<array<int>> overlap_m::EditDestinationRanks() {

  if (!DestinationRanksEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<overlap_m> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      overlap_m &OverlapM = *FloatingRef;
      MPI_Barrier(OverlapM.Comm_);
      OverlapM.DestinationRanksEvent_.Trigger();
      MPI_Barrier(OverlapM.Comm_);
    };
    DestinationRanksEditor_.Activate(std::move(DeactivateFunc));
  }

  return DestinationRanksEditor_.Edit(DestinationRanks_);

}

void overlap_m::RestoreDestinationRanks() {

  OVK_DEBUG_ASSERT(DestinationRanksEditor_.Active(), "Unable to restore destination ranks; not "
    "currently being edited.");

  DestinationRanksEditor_.Restore();

}

}
