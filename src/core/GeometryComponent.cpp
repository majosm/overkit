// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/GeometryComponent.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/DomainBase.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Geometry.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Set.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <algorithm>
#include <string>
#include <utility>

namespace ovk {

namespace geometry_component_internal {

geometry_component_base::geometry_component_base(const core::domain_base &Domain, std::string
  &&Name):
  Context_(Domain.Context().GetFloatingRef()),
  Domain_(Domain.GetFloatingRef()),
  Name_(std::move(Name))
{
  MPI_Barrier(Domain.Comm());
}

geometry_component_base::~geometry_component_base() noexcept {

  if (Context_) {
    const core::domain_base &Domain = *Domain_;
    MPI_Barrier(Domain.Comm());
    core::logger &Logger = Context_->core_Logger();
    Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroyed geometry component %s.%s.",
      Domain.Name(), *Name_);
  }

}

}

geometry_component::geometry_component(const core::domain_base &Domain, params Params):
  geometry_component_base(Domain, std::move(*Params.Name_))
{

  floating_ref<geometry_component> FloatingRef = FloatingRefGenerator_.Generate(*this);

  GridEventListener_ = Domain.AddGridEventListener([FloatingRef](int GridID, grid_event_flags Flags,
    bool LastInSequence) {
    geometry_component &GeometryComponent = *FloatingRef;
    grid_event_flags &AccumulatedFlags = GeometryComponent.GridEventFlags_.Fetch(GridID,
      grid_event_flags::NONE);
    AccumulatedFlags |= Flags;
    if (LastInSequence) {
      GeometryComponent.OnGridEvent_();
    }
  });

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Created geometry component %s.%s.", Domain.Name(),
    *Name_);

}

void geometry_component::OnGridEvent_() {

  DestroyGeometriesForDyingGrids_();

  GridEventFlags_.Clear();

}

void geometry_component::DestroyGeometriesForDyingGrids_() {

  set<int> DyingGridIDs;

  for (auto &EventEntry : GridEventFlags_) {
    int GridID = EventEntry.Key();
    grid_event_flags EventFlags = EventEntry.Value();
    if ((EventFlags & grid_event_flags::DESTROY) != grid_event_flags::NONE) {
      DyingGridIDs.Insert(GridID);
    }
  }

  DestroyGeometries(DyingGridIDs);

}

void geometry_component::StartEdit() {

  // Nothing to do here

}

void geometry_component::EndEdit() {

  SyncEdits_();

}

void geometry_component::SyncEdits_() {

  const core::domain_base &Domain = *Domain_;

  int HasEdits = false;
  for (auto &LocalEntry : Locals_) {
    local &Local = LocalEntry.Value();
    if (Local.EventFlags != geometry_event_flags::NONE) {
      HasEdits = true;
      break;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &HasEdits, 1, MPI_INT, MPI_MAX, Domain.Comm());

  if (!HasEdits) return;

  int NumGeometries = GeometryRecords_.Count();

  map<int,int> GridIDsToIndex;

  int NextIndex = 0;
  for (int GridID : GeometryRecords_.Keys()) {
    GridIDsToIndex.Insert(GridID, NextIndex);
    ++NextIndex;
  }

  array<geometry_event_flags> AllGeometryEventFlags({NumGeometries}, geometry_event_flags::NONE);

  for (auto &LocalEntry : Locals_) {
    int GridID = LocalEntry.Key();
    local &Local = LocalEntry.Value();
    int iGeometry = GridIDsToIndex(GridID);
    AllGeometryEventFlags(iGeometry) |= Local.EventFlags;
  }

  MPI_Allreduce(MPI_IN_PLACE, AllGeometryEventFlags.Data(), NumGeometries, MPI_INT, MPI_BOR,
    Domain.Comm());

  int NumTriggers = 0;
  for (auto EventFlags : AllGeometryEventFlags) {
    if (EventFlags != geometry_event_flags::NONE) ++NumTriggers;
  }

  int iTrigger = 0;
  for (int GridID : GeometryRecords_.Keys()) {
    int iGeometry = GridIDsToIndex(GridID);
    geometry_event_flags EventFlags = AllGeometryEventFlags(iGeometry);
    if (EventFlags != geometry_event_flags::NONE) {
      GeometryEvent_.Trigger(GridID, EventFlags, iTrigger == NumTriggers-1);
      ++iTrigger;
    }
  }

  MPI_Barrier(Domain.Comm());

  for (auto &LocalEntry : Locals_) {
    local &Local = LocalEntry.Value();
    Local.EventFlags = geometry_event_flags::NONE;
  }

}

bool geometry_component::GeometryExists(int GridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);

  return GeometryRecords_.Contains(GridID);

}

void geometry_component::CreateGeometry(int GridID, optional<geometry::params> MaybeParams) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(!GeometryExists(GridID), "Geometry %i already exists.", GridID);

  SyncEdits_();

  core::logger &Logger = Context_->core_Logger();

  const grid_info &GridInfo = Domain.GridInfo(GridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating geometry %s.%s...", Domain.Name(),
    GridInfo.Name());

  const std::shared_ptr<context> &SharedContext = Domain.SharedContext();

  if (GridInfo.IsLocal()) {
    const grid &Grid = Domain.Grid(GridID);
    if (MaybeParams.Present()) {
      geometry Geometry = core::CreateGeometry(SharedContext, Grid, MaybeParams.Release());
      Locals_.Insert(GridID, std::move(Geometry));
    } else {
      geometry Geometry = core::CreateGeometry(SharedContext, Grid);
      Locals_.Insert(GridID, std::move(Geometry));
    }
  }

  GeometryRecords_.Insert(GridID);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating geometry %s.%s.", Domain.Name(),
    GridInfo.Name());

  GeometryEvent_.Trigger(GridID, geometry_event_flags::CREATE, true);

  MPI_Barrier(Domain.Comm());

}

void geometry_component::CreateGeometries(array_view<const int> GridIDs, array<optional<
  geometry::params>> MaybeParams) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int NumCreates = GridIDs.Count();

  OVK_DEBUG_ASSERT(MaybeParams.Count() == NumCreates || MaybeParams.Count() == 0, "Incorrect "
    "params array size.");

  if (OVK_DEBUG) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int GridID = GridIDs(iCreate);
      OVK_DEBUG_ASSERT(GridID >= 0, "Invalid M grid ID.");
      OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
      OVK_DEBUG_ASSERT(!GeometryExists(GridID), "Geometry %i already exists.", GridID);
    }
  }

  SyncEdits_();

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int GridID = GridIDs(iCreate);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating geometry %s.%s...", Domain.Name(),
        GridInfo.Name());
    }
  }

  const std::shared_ptr<context> &SharedContext = Domain.SharedContext();

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int GridID = GridIDs(iCreate);
    if (Domain.GridIsLocal(GridID)) {
      const grid &Grid = Domain.Grid(GridID);
      if (MaybeParams.Count() > 0 && MaybeParams(iCreate).Present()) {
        geometry Geometry = core::CreateGeometry(SharedContext, Grid, MaybeParams(iCreate)
          .Release());
        Locals_.Insert(GridID, std::move(Geometry));
      } else {
        geometry Geometry = core::CreateGeometry(SharedContext, Grid);
        Locals_.Insert(GridID, std::move(Geometry));
      }
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int GridID = GridIDs(iCreate);
    GeometryRecords_.Insert(GridID);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int GridID = GridIDs(iCreate);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating geometry %s.%s.",
        Domain.Name(), GridInfo.Name());
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int GridID = GridIDs(iCreate);
    GeometryEvent_.Trigger(GridID, geometry_event_flags::CREATE, iCreate == NumCreates-1);
  }

  MPI_Barrier(Domain.Comm());

}

void geometry_component::DestroyGeometry(int GridID) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(GeometryExists(GridID), "Geometry %i does not exist.", GridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    int Editing = 0;
    if (GridInfo.IsLocal()) {
      Editing = Editing || Locals_(GridID).Editor.Active();
    }
    MPI_Allreduce(MPI_IN_PLACE, &Editing, 1, MPI_INT, MPI_LOR, Domain.Comm());
    OVK_DEBUG_ASSERT(!Editing, "Cannot destroy geometry %s.%s; still being edited.", Domain.Name(),
      GridInfo.Name());
  }

  SyncEdits_();

  GeometryEvent_.Trigger(GridID, geometry_event_flags::DESTROY, true);

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  const grid_info &GridInfo = Domain.GridInfo(GridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying geometry %s.%s...", Domain.Name(),
    GridInfo.Name());

  Locals_.Erase(GridID);

  GeometryRecords_.Erase(GridID);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying geometry %s.%s.", Domain.Name(),
    GridInfo.Name());

}

void geometry_component::DestroyGeometries(array_view<const int> GridIDs) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int NumDestroys = GridIDs.Count();

  if (OVK_DEBUG) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
      OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
      OVK_DEBUG_ASSERT(GeometryExists(GridID), "Geometry %i does not exist.", GridID);
    }
  }

  if (OVK_DEBUG) {
    array<int> Editing({NumDestroys}, 0);
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      if (GridInfo.IsLocal()) {
        Editing(iDestroy) = Editing(iDestroy) || Locals_(GridID).Editor.Active();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, Editing.Data(), NumDestroys, MPI_INT, MPI_LOR, Domain.Comm());
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      OVK_DEBUG_ASSERT(!Editing(iDestroy), "Cannot destroy geometry %s.%s; still being edited.",
        Domain.Name(), GridInfo.Name());
    }
  }

  SyncEdits_();

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int GridID = GridIDs(iDestroy);
    GeometryEvent_.Trigger(GridID, geometry_event_flags::DESTROY, iDestroy == NumDestroys-1);
  }

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying geometry %s.%s...", Domain.Name(),
        GridInfo.Name());
    }
  }

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int GridID = GridIDs(iDestroy);
    Locals_.Erase(GridID);
  }

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int GridID = GridIDs(iDestroy);
    GeometryRecords_.Erase(GridID);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying geometry %s.%s.",
        Domain.Name(), GridInfo.Name());
    }
  }

}

void geometry_component::ClearGeometries() {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Clearing geometries...");

  int NumGeometries = GeometryCount();

  if (OVK_DEBUG) {
    array<int> Editing({NumGeometries}, 0);
    int iGeometry = 0;
    for (int GridID : GeometryRecords_.Keys()) {
      if (Domain.GridIsLocal(GridID)) {
        Editing(iGeometry) = Editing(iGeometry) || Locals_(GridID).Editor.Active();
      }
      ++iGeometry;
    }
    MPI_Allreduce(MPI_IN_PLACE, Editing.Data(), NumGeometries, MPI_INT, MPI_LOR, Domain.Comm());
    iGeometry = 0;
    for (int GridID : GeometryRecords_.Keys()) {
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      OVK_DEBUG_ASSERT(!Editing(iGeometry), "Cannot destroy geometry %s.%s; still being edited.",
        Domain.Name(), GridInfo.Name());
      ++iGeometry;
    }
  }

  SyncEdits_();

  int iGeometry = 0;
  for (int GridID : GeometryRecords_.Keys()) {
    GeometryEvent_.Trigger(GridID, geometry_event_flags::DESTROY, iGeometry == NumGeometries-1);
    ++iGeometry;
  }

  MPI_Barrier(Domain.Comm());

  Locals_.Clear();

  GeometryRecords_.Clear();

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done clearing geometries.");

}

const geometry &geometry_component::Geometry(int GridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "Grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  return Locals_(GridID).Geometry;

}

bool geometry_component::EditingGeometry(int GridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(GeometryExists(GridID), "Geometry %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "Grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  const local &Local = Locals_(GridID);
  const editor &Editor = Local.Editor;

  return Editor.Active();

}

edit_handle<geometry> geometry_component::EditGeometry(int GridID) {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(GeometryExists(GridID), "Geometry %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "Grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local &Local = Locals_(GridID);
  geometry &Geometry = Local.Geometry;
  editor &Editor = Local.Editor;

  if (!Editor.Active()) {
    floating_ref<const grid> GridRef = Domain.Grid(GridID).GetFloatingRef();
    MPI_Barrier(GridRef->Comm());
    auto DeactivateFunc = [GridRef] { MPI_Barrier(GridRef->Comm()); };
    Editor.Activate(std::move(DeactivateFunc));
  }

  return Editor.Edit(Geometry);

}

void geometry_component::RestoreGeometry(int GridID) {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(GeometryExists(GridID), "Geometry %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "Grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local &Local = Locals_(GridID);
  editor &Editor = Local.Editor;

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore geometry %s.%s; not currently being "
      "edited.", Domain.Name(), GridInfo.Name());
  }

  Editor.Restore();

}

geometry_component::local::local(geometry Geometry_):
  Geometry(std::move(Geometry_)),
  EventFlags(geometry_event_flags::NONE)
{

  floating_ref<geometry_event_flags> EventFlagsRef = FloatingRefGenerator.Generate(EventFlags);

  CoordsEventListener = Geometry.AddCoordsEventListener([EventFlagsRef] {
    *EventFlagsRef |= geometry_event_flags::EDIT_COORDS;
  });

}

geometry_component::params &geometry_component::params::SetName(std::string Name) {

  Name_ = std::move(Name);

  return *this;

}

}
