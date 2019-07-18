// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/StateComponent.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/DomainBase.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/IDMap.hpp"
#include "ovk/core/IDSet.hpp"
#include "ovk/core/State.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <algorithm>
#include <string>
#include <utility>

namespace ovk {

namespace state_component_internal {

state_component_base::state_component_base(const core::domain_base &Domain, std::string
  &&Name):
  Context_(Domain.Context().GetFloatingRef()),
  Domain_(Domain.GetFloatingRef()),
  Name_(std::move(Name))
{
  MPI_Barrier(Domain.Comm());
}

state_component_base::~state_component_base() noexcept {

  if (Context_) {
    const core::domain_base &Domain = *Domain_;
    MPI_Barrier(Domain.Comm());
    core::logger &Logger = Context_->core_Logger();
    Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroyed state component %s.%s.",
      Domain.Name(), *Name_);
  }

}

}

state_component::state_component(const core::domain_base &Domain, params Params):
  state_component_base(Domain, std::move(*Params.Name_))
{

  floating_ref<state_component> FloatingRef = FloatingRefGenerator_.Generate(*this);

  GridEventListener_ = Domain.AddGridEventListener([FloatingRef](int GridID, grid_event_flags Flags,
    bool LastInSequence) {
    state_component &StateComponent = *FloatingRef;
    grid_event_flags &AccumulatedFlags = StateComponent.GridEventFlags_.Fetch(GridID,
      grid_event_flags::NONE);
    AccumulatedFlags |= Flags;
    if (LastInSequence) {
      StateComponent.OnGridEvent_();
    }
  });

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Created state component %s.%s.", Domain.Name(),
    *Name_);

}

void state_component::OnGridEvent_() {

  DestroyStatesForDyingGrids_();

  GridEventFlags_.Clear();

}

void state_component::DestroyStatesForDyingGrids_() {

  id_set<1> DyingGridIDs;

  for (auto &EventEntry : GridEventFlags_) {
    int GridID = EventEntry.Key();
    grid_event_flags EventFlags = EventEntry.Value();
    if ((EventFlags & grid_event_flags::DESTROY) != grid_event_flags::NONE) {
      DyingGridIDs.Insert(GridID);
    }
  }

  if (DyingGridIDs.Count() > 0) {
    array<int> GridIDs({DyingGridIDs.Count()}, DyingGridIDs.Begin());
    DestroyStates(GridIDs);
  }

}

void state_component::StartEdit() {

  // Nothing to do here

}

void state_component::EndEdit() {

  SyncEdits_();

}

void state_component::SyncEdits_() {

  const core::domain_base &Domain = *Domain_;

  int HasEdits = false;
  for (auto &LocalEntry : Locals_) {
    local &Local = LocalEntry.Value();
    if (Local.EventFlags != state_event_flags::NONE) {
      HasEdits = true;
      break;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &HasEdits, 1, MPI_INT, MPI_MAX, Domain.Comm());

  if (!HasEdits) return;

  int NumStates = StateRecords_.Count();

  id_map<1,int> GridIDsToIndex;

  int NextIndex = 0;
  for (int GridID : StateRecords_.Keys()) {
    GridIDsToIndex.Insert(GridID, NextIndex);
    ++NextIndex;
  }

  array<state_event_flags> AllStateEventFlags({NumStates}, state_event_flags::NONE);

  for (auto &LocalEntry : Locals_) {
    int GridID = LocalEntry.Key();
    local &Local = LocalEntry.Value();
    int iState = GridIDsToIndex(GridID);
    AllStateEventFlags(iState) |= Local.EventFlags;
  }

  MPI_Allreduce(MPI_IN_PLACE, AllStateEventFlags.Data(), NumStates, MPI_INT, MPI_BOR,
    Domain.Comm());

  int NumTriggers = 0;
  for (auto EventFlags : AllStateEventFlags) {
    if (EventFlags != state_event_flags::NONE) ++NumTriggers;
  }

  int iTrigger = 0;
  for (int GridID : StateRecords_.Keys()) {
    int iState = GridIDsToIndex(GridID);
    state_event_flags EventFlags = AllStateEventFlags(iState);
    if (EventFlags != state_event_flags::NONE) {
      StateEvent_.Trigger(GridID, EventFlags, iTrigger == NumTriggers-1);
      ++iTrigger;
    }
  }

  MPI_Barrier(Domain.Comm());

  for (auto &LocalEntry : Locals_) {
    local &Local = LocalEntry.Value();
    Local.EventFlags = state_event_flags::NONE;
  }

}

int state_component::StateCount() const {

  return StateRecords_.Count();

}

const id_set<1> &state_component::StateIDs() const {

  return StateRecords_.Keys();

}

bool state_component::StateExists(int GridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);

  return StateRecords_.Contains(GridID);

}

void state_component::CreateState(int GridID, optional<state::params> MaybeParams) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(!StateExists(GridID), "State %i already exists.", GridID);

  SyncEdits_();

  core::logger &Logger = Context_->core_Logger();

  const grid_info &GridInfo = Domain.GridInfo(GridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating state %s.%s...", Domain.Name(),
    GridInfo.Name());

  const std::shared_ptr<context> &SharedContext = Domain.SharedContext();

  if (GridInfo.IsLocal()) {
    const grid &Grid = Domain.Grid(GridID);
    if (MaybeParams.Present()) {
      state State = core::CreateState(SharedContext, Grid, MaybeParams.Release());
      Locals_.Insert(GridID, std::move(State));
    } else {
      state State = core::CreateState(SharedContext, Grid);
      Locals_.Insert(GridID, std::move(State));
    }
  }

  StateRecords_.Insert(GridID);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating state %s.%s.", Domain.Name(),
    GridInfo.Name());

  StateEvent_.Trigger(GridID, state_event_flags::CREATE, true);

  MPI_Barrier(Domain.Comm());

}

void state_component::CreateStates(array_view<const int> GridIDs, array<optional<
  state::params>> MaybeParams) {

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
      OVK_DEBUG_ASSERT(!StateExists(GridID), "State %i already exists.", GridID);
    }
  }

  SyncEdits_();

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int GridID = GridIDs(iCreate);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating state %s.%s...", Domain.Name(),
        GridInfo.Name());
    }
  }

  const std::shared_ptr<context> &SharedContext = Domain.SharedContext();

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int GridID = GridIDs(iCreate);
    if (Domain.GridIsLocal(GridID)) {
      const grid &Grid = Domain.Grid(GridID);
      if (MaybeParams.Count() > 0 && MaybeParams(iCreate).Present()) {
        state State = core::CreateState(SharedContext, Grid, MaybeParams(iCreate)
          .Release());
        Locals_.Insert(GridID, std::move(State));
      } else {
        state State = core::CreateState(SharedContext, Grid);
        Locals_.Insert(GridID, std::move(State));
      }
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int GridID = GridIDs(iCreate);
    StateRecords_.Insert(GridID);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int GridID = GridIDs(iCreate);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating state %s.%s.",
        Domain.Name(), GridInfo.Name());
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int GridID = GridIDs(iCreate);
    StateEvent_.Trigger(GridID, state_event_flags::CREATE, iCreate == NumCreates-1);
  }

  MPI_Barrier(Domain.Comm());

}

void state_component::DestroyState(int GridID) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(StateExists(GridID), "State %i does not exist.", GridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    int Editing = 0;
    if (GridInfo.IsLocal()) {
      Editing = Editing || Locals_(GridID).Editor.Active();
    }
    MPI_Allreduce(MPI_IN_PLACE, &Editing, 1, MPI_INT, MPI_LOR, Domain.Comm());
    OVK_DEBUG_ASSERT(!Editing, "Cannot destroy state %s.%s; still being edited.", Domain.Name(),
      GridInfo.Name());
  }

  SyncEdits_();

  StateEvent_.Trigger(GridID, state_event_flags::DESTROY, true);

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  const grid_info &GridInfo = Domain.GridInfo(GridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying state %s.%s...", Domain.Name(),
    GridInfo.Name());

  Locals_.Erase(GridID);

  StateRecords_.Erase(GridID);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying state %s.%s.", Domain.Name(),
    GridInfo.Name());

}

void state_component::DestroyStates(array_view<const int> GridIDs) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int NumDestroys = GridIDs.Count();

  if (OVK_DEBUG) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
      OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
      OVK_DEBUG_ASSERT(StateExists(GridID), "State %i does not exist.", GridID);
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
      OVK_DEBUG_ASSERT(!Editing(iDestroy), "Cannot destroy state %s.%s; still being edited.",
        Domain.Name(), GridInfo.Name());
    }
  }

  SyncEdits_();

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int GridID = GridIDs(iDestroy);
    StateEvent_.Trigger(GridID, state_event_flags::DESTROY, iDestroy == NumDestroys-1);
  }

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying state %s.%s...", Domain.Name(),
        GridInfo.Name());
    }
  }

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int GridID = GridIDs(iDestroy);
    Locals_.Erase(GridID);
  }

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int GridID = GridIDs(iDestroy);
    StateRecords_.Erase(GridID);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int GridID = GridIDs(iDestroy);
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying state %s.%s.", Domain.Name(),
        GridInfo.Name());
    }
  }

}

void state_component::ClearStates() {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Clearing states...");

  int NumStates = StateCount();

  if (OVK_DEBUG) {
    array<int> Editing({NumStates}, 0);
    int iState = 0;
    for (int GridID : StateRecords_.Keys()) {
      if (Domain.GridIsLocal(GridID)) {
        Editing(iState) = Editing(iState) || Locals_(GridID).Editor.Active();
      }
      ++iState;
    }
    MPI_Allreduce(MPI_IN_PLACE, Editing.Data(), NumStates, MPI_INT, MPI_LOR, Domain.Comm());
    iState = 0;
    for (int GridID : StateRecords_.Keys()) {
      const grid_info &GridInfo = Domain.GridInfo(GridID);
      OVK_DEBUG_ASSERT(!Editing(iState), "Cannot destroy state %s.%s; still being edited.",
        Domain.Name(), GridInfo.Name());
      ++iState;
    }
  }

  SyncEdits_();

  int iState = 0;
  for (int GridID : StateRecords_.Keys()) {
    StateEvent_.Trigger(GridID, state_event_flags::DESTROY, iState == NumStates-1);
    ++iState;
  }

  MPI_Barrier(Domain.Comm());

  Locals_.Clear();

  StateRecords_.Clear();

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done clearing states.");

}

int state_component::LocalStateCount() const {

  return Locals_.Count();

}

const id_set<1> &state_component::LocalStateIDs() const {

  return Locals_.Keys();

}

const state &state_component::State(int GridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "Grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  return Locals_(GridID).State;

}

bool state_component::EditingState(int GridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(StateExists(GridID), "State %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "Grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  const local &Local = Locals_(GridID);
  const editor &Editor = Local.Editor;

  return Editor.Active();

}

edit_handle<state> state_component::EditState(int GridID) {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(StateExists(GridID), "State %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "Grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local &Local = Locals_(GridID);
  state &State = Local.State;
  editor &Editor = Local.Editor;

  if (!Editor.Active()) {
    floating_ref<const grid> GridRef = Domain.Grid(GridID).GetFloatingRef();
    MPI_Barrier(GridRef->Comm());
    auto DeactivateFunc = [GridRef] { MPI_Barrier(GridRef->Comm()); };
    Editor.Activate(std::move(DeactivateFunc));
  }

  return Editor.Edit(State);

}

void state_component::RestoreState(int GridID) {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(StateExists(GridID), "State %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "Grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local &Local = Locals_(GridID);
  editor &Editor = Local.Editor;

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore state %s.%s; not currently being "
      "edited.", Domain.Name(), GridInfo.Name());
  }

  Editor.Restore();

}

state_component::local::local(state State_):
  State(std::move(State_)),
  EventFlags(state_event_flags::NONE)
{

  floating_ref<state_event_flags> EventFlagsRef = FloatingRefGenerator.Generate(EventFlags);

  CoordsEventListener = State.AddFlagsEventListener([EventFlagsRef] {
    *EventFlagsRef |= state_event_flags::EDIT_FLAGS;
  });

}

state_component::params &state_component::params::SetName(std::string Name) {

  Name_ = std::move(Name);

  return *this;

}

}
