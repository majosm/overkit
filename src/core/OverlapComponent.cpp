// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/OverlapComponent.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/DomainBase.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/ElemMap.hpp"
#include "ovk/core/ElemSet.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/OverlapM.hpp"
#include "ovk/core/OverlapN.hpp"
#include "ovk/core/Set.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <algorithm>
#include <string>
#include <utility>

namespace ovk {

namespace overlap_component_internal {

overlap_component_base::overlap_component_base(const core::domain_base &Domain, std::string &&Name):
  Context_(Domain.Context().GetFloatingRef()),
  Domain_(Domain.GetFloatingRef()),
  Name_(std::move(Name))
{
  MPI_Barrier(Domain.Comm());
}

overlap_component_base::~overlap_component_base() noexcept {

  if (Context_) {
    const core::domain_base &Domain = *Domain_;
    MPI_Barrier(Domain.Comm());
    core::logger &Logger = Context_->core_Logger();
    Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroyed overlap component %s.%s.",
      Domain.Name(), *Name_);
  }

}

}

overlap_component::overlap_component(const core::domain_base &Domain, params Params):
  overlap_component_base(Domain, std::move(*Params.Name_))
{

  floating_ref<overlap_component> FloatingRef = FloatingRefGenerator_.Generate(*this);

  GridEventListener_ = Domain.AddGridEventListener([FloatingRef](int GridID, grid_event_flags Flags,
    bool LastInSequence) {
    overlap_component &OverlapComponent = *FloatingRef;
    grid_event_flags &AccumulatedFlags = OverlapComponent.GridEventFlags_.Fetch(GridID,
      grid_event_flags::NONE);
    AccumulatedFlags |= Flags;
    if (LastInSequence) {
      OverlapComponent.OnGridEvent_();
    }
  });

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Created overlap component %s.%s.", Domain.Name(),
    *Name_);

}

void overlap_component::OnGridEvent_() {

  DestroyOverlapsForDyingGrids_();

  GridEventFlags_.Clear();

}

void overlap_component::DestroyOverlapsForDyingGrids_() {

  set<int> DyingGridIDs;

  for (auto &EventEntry : GridEventFlags_) {
    int GridID = EventEntry.Key();
    grid_event_flags EventFlags = EventEntry.Value();
    if ((EventFlags & grid_event_flags::DESTROY) != grid_event_flags::NONE) {
      DyingGridIDs.Insert(GridID);
    }
  }

  elem_set<int,2> DyingOverlapIDs;

  for (auto &IDPair : OverlapRecords_.Keys()) {
    if (DyingGridIDs.Contains(IDPair(0)) || DyingGridIDs.Contains(IDPair(1))) {
      DyingOverlapIDs.Insert(IDPair);
    }
  }

  DestroyOverlaps(DyingOverlapIDs);

}

void overlap_component::StartEdit() {

  // Nothing to do here

}

void overlap_component::EndEdit() {

  SyncEdits_();

}

void overlap_component::SyncEdits_() {

  const core::domain_base &Domain = *Domain_;

  int HasEdits = false;
  for (auto &LocalMEntry : LocalMs_) {
    local_m &LocalM = LocalMEntry.Value();
    if (LocalM.EventFlags != overlap_event_flags::NONE) {
      HasEdits = true;
      goto done_looping;
    }
  }
  for (auto &LocalNEntry : LocalNs_) {
    local_n &LocalN = LocalNEntry.Value();
    if (LocalN.EventFlags != overlap_event_flags::NONE) {
      HasEdits = true;
      goto done_looping;
    }
  }
  done_looping:;
  MPI_Allreduce(MPI_IN_PLACE, &HasEdits, 1, MPI_INT, MPI_MAX, Domain.Comm());

  if (!HasEdits) return;

  int NumOverlaps = OverlapRecords_.Count();

  elem_map<int,2,int> GridIDsToIndex;

  int NextIndex = 0;
  for (auto &IDPair : OverlapRecords_.Keys()) {
    GridIDsToIndex.Insert({IDPair,NextIndex});
    ++NextIndex;
  }

  array<overlap_event_flags> AllOverlapEventFlags({NumOverlaps},
    overlap_event_flags::NONE);

  for (auto &LocalMEntry : LocalMs_) {
    local_m &LocalM = LocalMEntry.Value();
    int iOverlap = GridIDsToIndex(LocalMEntry.Key());
    AllOverlapEventFlags(iOverlap) |= LocalM.EventFlags;
  }

  for (auto &LocalNEntry : LocalNs_) {
    local_n &LocalN = LocalNEntry.Value();
    int iOverlap = GridIDsToIndex(LocalNEntry.Key());
    AllOverlapEventFlags(iOverlap) |= LocalN.EventFlags;
  }

  MPI_Allreduce(MPI_IN_PLACE, AllOverlapEventFlags.Data(), NumOverlaps, MPI_INT, MPI_BOR,
    Domain.Comm());

  int NumTriggers = 0;
  for (auto EventFlags : AllOverlapEventFlags) {
    if (EventFlags != overlap_event_flags::NONE) ++NumTriggers;
  }
  for (auto &IDPair : OverlapRecords_.Keys()) {
    int iOverlap = GridIDsToIndex(IDPair);
    overlap_event_flags EventFlags = AllOverlapEventFlags(iOverlap);
    if (EventFlags != overlap_event_flags::NONE) {
      OverlapEvent_.Trigger(IDPair, EventFlags, --NumTriggers == 0);
    }
  }

  MPI_Barrier(Domain.Comm());

  for (auto &LocalMEntry : LocalMs_) {
    local_m &LocalM = LocalMEntry.Value();
    LocalM.EventFlags = overlap_event_flags::NONE;
  }

  for (auto &LocalNEntry : LocalNs_) {
    local_n &LocalN = LocalNEntry.Value();
    LocalN.EventFlags = overlap_event_flags::NONE;
  }

}

int overlap_component::OverlapCount() const {

  return OverlapRecords_.Count();

}

const elem_set<int,2> &overlap_component::OverlapIDs() const {

  return OverlapRecords_.Keys();

}

bool overlap_component::OverlapExists(const elem<int,2> &GridIDPair) const {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = GridIDPair(0);
    int NGridID = GridIDPair(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  }

  return OverlapRecords_.Contains(GridIDPair);

}

void overlap_component::CreateOverlap(const elem<int,2> &GridIDPair) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int MGridID = GridIDPair(0);
  int NGridID = GridIDPair(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(!OverlapExists(GridIDPair), "Overlap (%i,%i) already exists.", MGridID, NGridID);

  SyncEdits_();

  core::logger &Logger = Context_->core_Logger();

  const grid_info &MGridInfo = Domain.GridInfo(MGridID);
  const grid_info &NGridInfo = Domain.GridInfo(NGridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating overlap %s.(%s,%s)...", Domain.Name(),
    MGridInfo.Name(), NGridInfo.Name());

  const std::shared_ptr<context> &SharedContext = Domain.SharedContext();

  if (MGridInfo.IsLocal()) {
    const grid &MGrid = Domain.Grid(MGridID);
    overlap_m OverlapM = core::CreateOverlapM(SharedContext, MGrid, NGridInfo);
    LocalMs_.Insert(GridIDPair, std::move(OverlapM));
  }

  if (NGridInfo.IsLocal()) {
    const grid &NGrid = Domain.Grid(NGridID);
    overlap_n OverlapN = core::CreateOverlapN(SharedContext, NGrid, MGridInfo);
    LocalNs_.Insert(GridIDPair, std::move(OverlapN));
  }

  OverlapRecords_.Insert(GridIDPair);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating overlap %s.(%s,%s).", Domain.Name(),
    MGridInfo.Name(), NGridInfo.Name());

  OverlapEvent_.Trigger(GridIDPair, overlap_event_flags::CREATE, true);

  MPI_Barrier(Domain.Comm());

}

void overlap_component::CreateOverlaps(array_view<const elem<int,2>> GridIDPairs) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  if (OVK_DEBUG) {
    for (auto &IDPair : GridIDPairs) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
      OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
      OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
      OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
      OVK_DEBUG_ASSERT(!OverlapExists(IDPair), "Overlap (%i,%i) already exists.", MGridID, NGridID);
    }
  }

  SyncEdits_();

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (auto &IDPair : GridIDPairs) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating overlap %s.(%s,%s)...",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  const std::shared_ptr<context> &SharedContext = Domain.SharedContext();

  for (auto &IDPair : GridIDPairs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    if (MGridInfo.IsLocal()) {
      const grid &MGrid = Domain.Grid(MGridID);
      overlap_m OverlapM = core::CreateOverlapM(SharedContext, MGrid, NGridInfo);
      LocalMs_.Insert(IDPair, std::move(OverlapM));
    }
    if (NGridInfo.IsLocal()) {
      const grid &NGrid = Domain.Grid(NGridID);
      overlap_n OverlapN = core::CreateOverlapN(SharedContext, NGrid, MGridInfo);
      LocalNs_.Insert(IDPair, std::move(OverlapN));
    }
  }

  for (auto &IDPair : GridIDPairs) {
    OverlapRecords_.Insert(IDPair);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (auto &IDPair : GridIDPairs) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating overlap %s.(%s,%s).",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  int NumRemaining = GridIDPairs.Count();
  for (auto &IDPair : GridIDPairs) {
    OverlapEvent_.Trigger(IDPair, overlap_event_flags::CREATE, --NumRemaining == 0);
  }

  MPI_Barrier(Domain.Comm());

}

void overlap_component::DestroyOverlap(const elem<int,2> &GridIDPair) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int MGridID = GridIDPair(0);
  int NGridID = GridIDPair(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(GridIDPair), "Overlap (%i,%i) does not exist.", MGridID, NGridID);

  if (OVK_DEBUG) {
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    int Editing = 0;
    if (MGridInfo.IsLocal()) {
      Editing = Editing || LocalMs_(GridIDPair).Editor.Active();
    }
    if (NGridInfo.IsLocal()) {
      Editing = Editing || LocalNs_(GridIDPair).Editor.Active();
    }
    MPI_Allreduce(MPI_IN_PLACE, &Editing, 1, MPI_INT, MPI_LOR, Domain.Comm());
    OVK_DEBUG_ASSERT(!Editing, "Cannot destroy overlap %s.(%s,%s); still being edited.",
      Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
  }

  SyncEdits_();

  OverlapEvent_.Trigger(GridIDPair, overlap_event_flags::DESTROY, true);

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  const grid_info &MGridInfo = Domain.GridInfo(MGridID);
  const grid_info &NGridInfo = Domain.GridInfo(NGridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying overlap %s.(%s,%s)...", Domain.Name(),
    MGridInfo.Name(), NGridInfo.Name());

  LocalMs_.Erase(GridIDPair);
  LocalNs_.Erase(GridIDPair);

  OverlapRecords_.Erase(GridIDPair);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying overlap %s.(%s,%s).",
    Domain.Name(), MGridInfo.Name(), NGridInfo.Name());

}

void overlap_component::DestroyOverlaps(array_view<const elem<int,2>> GridIDPairs) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int NumDestroys = GridIDPairs.Count();

  if (OVK_DEBUG) {
    for (auto &IDPair : GridIDPairs) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
      OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
      OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
      OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
      OVK_DEBUG_ASSERT(OverlapExists(IDPair), "Overlap (%i,%i) does not exist.", MGridID, NGridID);
    }
  }

  if (OVK_DEBUG) {
    array<int> Editing({NumDestroys}, 0);
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      auto &IDPair = GridIDPairs(iDestroy);
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      if (Domain.GridIsLocal(MGridID)) {
        Editing(iDestroy) = Editing(iDestroy) || LocalMs_(IDPair).Editor.Active();
      }
      if (Domain.GridIsLocal(NGridID)) {
        Editing(iDestroy) = Editing(iDestroy) || LocalNs_(IDPair).Editor.Active();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, Editing.Data(), NumDestroys, MPI_INT, MPI_LOR, Domain.Comm());
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      auto &IDPair = GridIDPairs(iDestroy);
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      OVK_DEBUG_ASSERT(!Editing(iDestroy), "Cannot destroy overlap %s.(%s,%s); still being edited.",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  SyncEdits_();

  int NumRemaining = NumDestroys;
  for (auto &IDPair : GridIDPairs) {
    OverlapEvent_.Trigger(IDPair, overlap_event_flags::DESTROY, --NumRemaining == 0);
  }

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (auto &IDPair : GridIDPairs) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying overlap %s.(%s,%s)...",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  for (auto &IDPair : GridIDPairs) {
    LocalMs_.Erase(IDPair);
    LocalNs_.Erase(IDPair);
  }

  for (auto &IDPair : GridIDPairs) {
    OverlapRecords_.Erase(IDPair);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (auto &IDPair : GridIDPairs) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying overlap %s.(%s,%s).",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

}

void overlap_component::ClearOverlaps() {

  DestroyOverlaps(OverlapRecords_.Keys());

}

int overlap_component::LocalOverlapMCount() const {

  return LocalMs_.Count();

}

int overlap_component::LocalOverlapMCountForGrid(int MGridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);

  int NumLocalMsForGrid = 0;

  for (auto &LocalMEntry : LocalMs_) {
    if (LocalMEntry.Key()(0) == MGridID) ++NumLocalMsForGrid;
  }

  return NumLocalMsForGrid;

}

const elem_set<int,2> &overlap_component::LocalOverlapMIDs() const {

  return LocalMs_.Keys();

}

const overlap_m &overlap_component::OverlapM(const elem<int,2> &GridIDPair) const {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = GridIDPair(0);
    int NGridID = GridIDPair(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
    OVK_DEBUG_ASSERT(OverlapExists(GridIDPair), "Overlap (%i,%i) does not exist.", MGridID,
      NGridID);
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  return LocalMs_(GridIDPair).Overlap;

}

bool overlap_component::EditingOverlapM(const elem<int,2> &GridIDPair) const {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = GridIDPair(0);
    int NGridID = GridIDPair(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
    OVK_DEBUG_ASSERT(OverlapExists(GridIDPair), "Overlap (%i,%i) does not exist.", MGridID,
      NGridID);
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  const local_m &LocalM = LocalMs_(GridIDPair);
  const editor &Editor = LocalM.Editor;

  return Editor.Active();

}

edit_handle<overlap_m> overlap_component::EditOverlapM(const elem<int,2> &GridIDPair) {

  const core::domain_base &Domain = *Domain_;

  int MGridID = GridIDPair(0);
  int NGridID = GridIDPair(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(GridIDPair), "Overlap (%i,%i) does not exist.", MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_m &LocalM = LocalMs_(GridIDPair);
  overlap_m &OverlapM = LocalM.Overlap;
  editor &Editor = LocalM.Editor;

  if (!Editor.Active()) {
    floating_ref<const grid> GridRef = Domain.Grid(MGridID).GetFloatingRef();
    MPI_Barrier(GridRef->Comm());
    auto DeactivateFunc = [GridRef] { MPI_Barrier(GridRef->Comm()); };
    Editor.Activate(std::move(DeactivateFunc));
  }

  return Editor.Edit(OverlapM);

}

void overlap_component::RestoreOverlapM(const elem<int,2> &GridIDPair) {

  const core::domain_base &Domain = *Domain_;

  int MGridID = GridIDPair(0);
  int NGridID = GridIDPair(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(GridIDPair), "Overlap (%i,%i) does not exist.", MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_m &LocalM = LocalMs_(GridIDPair);
  editor &Editor = LocalM.Editor;

  if (OVK_DEBUG) {
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore overlap M %s.(%s,%s); not currently "
      "being edited.", Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
  }

  Editor.Restore();

}

int overlap_component::LocalOverlapNCount() const {

  return LocalNs_.Count();

}

int overlap_component::LocalOverlapNCountForGrid(int NGridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);

  int NumLocalNsForGrid = 0;

  for (auto &LocalNEntry : LocalNs_) {
    if (LocalNEntry.Key()(1) == NGridID) ++NumLocalNsForGrid;
  }

  return NumLocalNsForGrid;

}

const elem_set<int,2> &overlap_component::LocalOverlapNIDs() const {

  return LocalNs_.Keys();

}

const overlap_n &overlap_component::OverlapN(const elem<int,2> &GridIDPair) const {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = GridIDPair(0);
    int NGridID = GridIDPair(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
    OVK_DEBUG_ASSERT(OverlapExists(GridIDPair), "Overlap (%i,%i) does not exist.", MGridID,
      NGridID);
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  return LocalNs_(GridIDPair).Overlap;

}

bool overlap_component::EditingOverlapN(const elem<int,2> &GridIDPair) const {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = GridIDPair(0);
    int NGridID = GridIDPair(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
    OVK_DEBUG_ASSERT(OverlapExists(GridIDPair), "Overlap (%i,%i) does not exist.", MGridID,
      NGridID);
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  const local_n &LocalN = LocalNs_(GridIDPair);
  const editor &Editor = LocalN.Editor;

  return Editor.Active();

}

edit_handle<overlap_n> overlap_component::EditOverlapN(const elem<int,2> &GridIDPair) {

  const core::domain_base &Domain = *Domain_;

  int MGridID = GridIDPair(0);
  int NGridID = GridIDPair(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(GridIDPair), "Overlap (%i,%i) does not exist.", MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_n &LocalN = LocalNs_(GridIDPair);
  overlap_n &OverlapN = LocalN.Overlap;
  editor &Editor = LocalN.Editor;

  if (!Editor.Active()) {
    floating_ref<const grid> GridRef = Domain.Grid(NGridID).GetFloatingRef();
    MPI_Barrier(GridRef->Comm());
    auto DeactivateFunc = [GridRef] { MPI_Barrier(GridRef->Comm()); };
    Editor.Activate(std::move(DeactivateFunc));
  }

  return Editor.Edit(OverlapN);

}

void overlap_component::RestoreOverlapN(const elem<int,2> &GridIDPair) {

  const core::domain_base &Domain = *Domain_;

  int MGridID = GridIDPair(0);
  int NGridID = GridIDPair(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(GridIDPair), "Overlap (%i,%i) does not exist.", MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_n &LocalN = LocalNs_(GridIDPair);
  editor &Editor = LocalN.Editor;

  if (OVK_DEBUG) {
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore overlap N %s.(%s,%s); not currently "
      "being edited.", Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
  }

  Editor.Restore();

}

overlap_component::local_m::local_m(overlap_m Overlap_):
  Overlap(std::move(Overlap_)),
  EventFlags(overlap_event_flags::NONE)
{

  floating_ref<overlap_event_flags> EventFlagsRef = FloatingRefGenerator.Generate(EventFlags);

  ResizeEventListener = Overlap.AddResizeEventListener([EventFlagsRef] {
    *EventFlagsRef |= overlap_event_flags::RESIZE_M;
  });

  CellsEventListener = Overlap.AddCellsEventListener([EventFlagsRef] {
    *EventFlagsRef |= overlap_event_flags::EDIT_M_CELLS;
  });

  CoordsEventListener = Overlap.AddCoordsEventListener([EventFlagsRef] {
    *EventFlagsRef |= overlap_event_flags::EDIT_M_COORDS;
  });

  DestinationsEventListener = Overlap.AddDestinationsEventListener([EventFlagsRef] {
    *EventFlagsRef |= overlap_event_flags::EDIT_M_DESTINATIONS;
  });

  DestinationRanksEventListener = Overlap.AddDestinationRanksEventListener([EventFlagsRef] {
    *EventFlagsRef |= overlap_event_flags::EDIT_M_DESTINATIONS;
  });

}

overlap_component::local_n::local_n(overlap_n Overlap_):
  Overlap(std::move(Overlap_)),
  EventFlags(overlap_event_flags::NONE)
{

  floating_ref<overlap_event_flags> EventFlagsRef = FloatingRefGenerator.Generate(EventFlags);

  ResizeEventListener = Overlap.AddResizeEventListener([EventFlagsRef] {
    *EventFlagsRef |= overlap_event_flags::RESIZE_N;
  });

  PointsEventListener = Overlap.AddPointsEventListener([EventFlagsRef] {
    *EventFlagsRef |= overlap_event_flags::EDIT_N_POINTS;
  });

  SourcesEventListener = Overlap.AddSourcesEventListener([EventFlagsRef] {
    *EventFlagsRef |= overlap_event_flags::EDIT_N_SOURCES;
  });

  SourceRanksEventListener = Overlap.AddSourceRanksEventListener([EventFlagsRef] {
    *EventFlagsRef |= overlap_event_flags::EDIT_N_SOURCES;
  });

}

overlap_component::params &overlap_component::params::SetName(std::string Name) {

  Name_ = std::move(Name);

  return *this;

}

}
