// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ConnectivityM.hpp"

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

namespace connectivity_m_internal {

connectivity_m_base::connectivity_m_base(std::shared_ptr<context> &&Context, const grid &Grid,
  grid_info &&DestinationGridInfo):
  Context_(std::move(Context)),
  Grid_(&Grid),
  DestinationGridInfo_(std::move(DestinationGridInfo)),
  Comm_(Grid_->Comm())
{
  MPI_Barrier(Comm_);
}

connectivity_m_base::~connectivity_m_base() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Logger.LogDebug(Comm_.Rank() == 0, 0, "Destroyed M side of connectivity (%s,%s).",
      Grid_->Name(), DestinationGridInfo_.Name());
  }

}

}

connectivity_m::connectivity_m(std::shared_ptr<context> &&Context, const grid &Grid, grid_info
  &&DestinationGridInfo):
  connectivity_m_base(std::move(Context), Grid, std::move(DestinationGridInfo)),
  NumDims_(Grid_->Dimension()),
  Count_(0),
  MaxSize_(1),
  Extents_({{2,MAX_DIMS,0}}),
  Coords_({{MAX_DIMS,0}}),
  InterpCoefs_({{MAX_DIMS,0,0}}),
  Destinations_({{MAX_DIMS,0}}),
  DestinationRanks_({0})
{

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();
  Logger.LogDebug(Comm_.Rank() == 0, 0, "Created M side of connectivity (%s,%s).", Grid_->Name(),
    DestinationGridInfo_.Name());

}

connectivity_m::~connectivity_m() noexcept {

  if (Context_) {
    // Barrier before cleaning up
    MPI_Barrier(Comm_);
  }

}

connectivity_m connectivity_m::internal_Create(std::shared_ptr<context> &&Context, const grid &Grid,
  grid_info &&DestinationGridInfo) {

  return {std::move(Context), Grid, std::move(DestinationGridInfo)};

}

namespace core {

connectivity_m CreateConnectivityM(std::shared_ptr<context> Context, const grid &Grid, grid_info
  DestinationGridInfo) {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return connectivity_m::internal_Create(std::move(Context), Grid, std::move(DestinationGridInfo));

}

}

void connectivity_m::Resize(long long Count, int MaxSize) {

  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(MaxSize > 0, "Invalid max size.");

  MPI_Barrier(Comm_);

  OVK_DEBUG_ASSERT(!ExtentsEditor_.Active(), "Cannot resize while editing extents.");
  OVK_DEBUG_ASSERT(!CoordsEditor_.Active(), "Cannot resize while editing coords.");
  OVK_DEBUG_ASSERT(!InterpCoefsEditor_.Active(), "Cannot resize while editing interp coefs.");
  OVK_DEBUG_ASSERT(!DestinationsEditor_.Active(), "Cannot resize while editing destinations.");
  OVK_DEBUG_ASSERT(!DestinationRanksEditor_.Active(), "Cannot resize while editing destination "
    "ranks.");

  Count_ = Count;
  MaxSize_ = MaxSize;

  Extents_.Resize({{2,MAX_DIMS,Count}});
  for (long long iCell = 0; iCell < Count; ++iCell) {
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      Extents_(0,iDim,iCell) = 0;
      Extents_(1,iDim,iCell) = 0;
    }
    for (int iDim = NumDims_; iDim < MAX_DIMS; ++iDim) {
      Extents_(0,iDim,iCell) = 0;
      Extents_(1,iDim,iCell) = 1;
    }
  }

  Coords_.Resize({{MAX_DIMS,Count}}, 0.);

  InterpCoefs_.Resize({{MAX_DIMS,MaxSize,Count}});
  for (long long iCell = 0; iCell < Count; ++iCell) {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      InterpCoefs_(iDim,0,iCell) = 1.;
      for (int iPoint = 1; iPoint < MaxSize; ++iPoint) {
        InterpCoefs_(iDim,iPoint,iCell) = 0.;
      }
    }
  }

  Destinations_.Resize({{MAX_DIMS,Count}}, 0);

  DestinationRanks_.Resize({Count}, -1);

  MPI_Barrier(Comm_);

  ResizeEvent_.Trigger();
  ExtentsEvent_.Trigger();
  CoordsEvent_.Trigger();
  InterpCoefsEvent_.Trigger();
  DestinationsEvent_.Trigger();
  DestinationRanksEvent_.Trigger();

  MPI_Barrier(Comm_);

}

bool connectivity_m::EditingExtents() const {

  return ExtentsEditor_.Active();

}

edit_handle<array<int,3>> connectivity_m::EditExtents() {

  if (!ExtentsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<connectivity_m> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      connectivity_m &ConnectivityM = *FloatingRef;
      MPI_Barrier(ConnectivityM.Comm_);
      ConnectivityM.ExtentsEvent_.Trigger();
      MPI_Barrier(ConnectivityM.Comm_);
    };
    ExtentsEditor_.Activate(std::move(DeactivateFunc));
  }

  return ExtentsEditor_.Edit(Extents_);

}

void connectivity_m::RestoreExtents() {

  OVK_DEBUG_ASSERT(ExtentsEditor_.Active(), "Unable to restore extents; not currently being "
    "edited.");

  ExtentsEditor_.Restore();

}

bool connectivity_m::EditingCoords() const {

  return CoordsEditor_.Active();

}

edit_handle<array<double,2>> connectivity_m::EditCoords() {

  if (!CoordsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<connectivity_m> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      connectivity_m &ConnectivityM = *FloatingRef;
      MPI_Barrier(ConnectivityM.Comm_);
      ConnectivityM.CoordsEvent_.Trigger();
      MPI_Barrier(ConnectivityM.Comm_);
    };
    CoordsEditor_.Activate(std::move(DeactivateFunc));
  }

  return CoordsEditor_.Edit(Coords_);

}

void connectivity_m::RestoreCoords() {

  OVK_DEBUG_ASSERT(CoordsEditor_.Active(), "Unable to restore coords; not currently being edited.");

  CoordsEditor_.Restore();

}

bool connectivity_m::EditingInterpCoefs() const {

  return InterpCoefsEditor_.Active();

}

edit_handle<array<double,3>> connectivity_m::EditInterpCoefs() {

  if (!InterpCoefsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<connectivity_m> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      connectivity_m &ConnectivityM = *FloatingRef;
      MPI_Barrier(ConnectivityM.Comm_);
      ConnectivityM.InterpCoefsEvent_.Trigger();
      MPI_Barrier(ConnectivityM.Comm_);
    };
    InterpCoefsEditor_.Activate(std::move(DeactivateFunc));
  }

  return InterpCoefsEditor_.Edit(InterpCoefs_);

}

void connectivity_m::RestoreInterpCoefs() {

  OVK_DEBUG_ASSERT(InterpCoefsEditor_.Active(), "Unable to restore interp coefs; not currently "
    "being edited.");

  InterpCoefsEditor_.Restore();

}

bool connectivity_m::EditingDestinations() const {

  return DestinationsEditor_.Active();

}

edit_handle<array<int,2>> connectivity_m::EditDestinations() {

  if (!DestinationsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<connectivity_m> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      connectivity_m &ConnectivityM = *FloatingRef;
      MPI_Barrier(ConnectivityM.Comm_);
      ConnectivityM.DestinationsEvent_.Trigger();
      MPI_Barrier(ConnectivityM.Comm_);
    };
    DestinationsEditor_.Activate(std::move(DeactivateFunc));
  }

  return DestinationsEditor_.Edit(Destinations_);

}

void connectivity_m::RestoreDestinations() {

  OVK_DEBUG_ASSERT(DestinationsEditor_.Active(), "Unable to restore destinations; not currently "
    "being edited.");

  DestinationsEditor_.Restore();

}

bool connectivity_m::EditingDestinationRanks() const {

  return DestinationRanksEditor_.Active();

}

edit_handle<array<int>> connectivity_m::EditDestinationRanks() {

  if (!DestinationsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<connectivity_m> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      connectivity_m &ConnectivityM = *FloatingRef;
      MPI_Barrier(ConnectivityM.Comm_);
      ConnectivityM.DestinationRanksEvent_.Trigger();
      MPI_Barrier(ConnectivityM.Comm_);
    };
    DestinationRanksEditor_.Activate(std::move(DeactivateFunc));
  }

  return DestinationRanksEditor_.Edit(DestinationRanks_);

}

void connectivity_m::RestoreDestinationRanks() {

  OVK_DEBUG_ASSERT(DestinationRanksEditor_.Active(), "Unable to restore destination ranks; not "
    "currently being edited.");

  DestinationRanksEditor_.Restore();

}

}
