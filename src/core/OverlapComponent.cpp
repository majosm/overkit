// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/OverlapComponent.hpp"

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
#include "ovk/core/OverlapM.hpp"
#include "ovk/core/OverlapN.hpp"
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

  id_set<1> DyingGridIDs;

  for (auto &EventEntry : GridEventFlags_) {
    int GridID = EventEntry.Key(0);
    grid_event_flags EventFlags = EventEntry.Value();
    if ((EventFlags & grid_event_flags::DESTROY) != grid_event_flags::NONE) {
      DyingGridIDs.Insert(GridID);
    }
  }

  id_set<2> DyingOverlapIDs;

  for (auto &IDPair : OverlapRecords_.Keys()) {
    if (DyingGridIDs.Contains(IDPair(0)) || DyingGridIDs.Contains(IDPair(1))) {
      DyingOverlapIDs.Insert(IDPair);
    }
  }

  if (DyingOverlapIDs.Count() > 0) {

    array<int> MGridIDs, NGridIDs;

    MGridIDs.Reserve(DyingOverlapIDs.Count());
    NGridIDs.Reserve(DyingOverlapIDs.Count());

    for (auto &IDPair : DyingOverlapIDs) {
      MGridIDs.Append(IDPair(0));
      NGridIDs.Append(IDPair(1));
    }

    DestroyOverlaps(MGridIDs, NGridIDs);

  }

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

  id_map<2,int> GridIDsToIndex;

  int NextIndex = 0;
  for (auto &IDPair : OverlapRecords_.Keys()) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    GridIDsToIndex.Insert({MGridID,NGridID}, NextIndex);
    ++NextIndex;
  }

  array<overlap_event_flags> AllOverlapEventFlags({NumOverlaps},
    overlap_event_flags::NONE);

  for (auto &LocalMEntry : LocalMs_) {
    int MGridID = LocalMEntry.Key(0);
    int NGridID = LocalMEntry.Key(1);
    local_m &LocalM = LocalMEntry.Value();
    int iOverlap = GridIDsToIndex(MGridID,NGridID);
    AllOverlapEventFlags(iOverlap) |= LocalM.EventFlags;
  }

  for (auto &LocalNEntry : LocalNs_) {
    int MGridID = LocalNEntry.Key(0);
    int NGridID = LocalNEntry.Key(1);
    local_n &LocalN = LocalNEntry.Value();
    int iOverlap = GridIDsToIndex(MGridID,NGridID);
    AllOverlapEventFlags(iOverlap) |= LocalN.EventFlags;
  }

  MPI_Allreduce(MPI_IN_PLACE, AllOverlapEventFlags.Data(), NumOverlaps, MPI_INT, MPI_BOR,
    Domain.Comm());

  int NumTriggers = 0;
  for (auto EventFlags : AllOverlapEventFlags) {
    if (EventFlags != overlap_event_flags::NONE) ++NumTriggers;
  }

  int iTrigger = 0;
  for (auto &IDPair : OverlapRecords_.Keys()) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    int iOverlap = GridIDsToIndex(MGridID,NGridID);
    overlap_event_flags EventFlags = AllOverlapEventFlags(iOverlap);
    if (EventFlags != overlap_event_flags::NONE) {
      OverlapEvent_.Trigger(MGridID, NGridID, EventFlags, iTrigger == NumTriggers-1);
      ++iTrigger;
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

const id_set<2> &overlap_component::OverlapIDs() const {

  return OverlapRecords_.Keys();

}

bool overlap_component::OverlapExists(int MGridID, int NGridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);

  return OverlapRecords_.Contains(MGridID,NGridID);

}

void overlap_component::CreateOverlap(int MGridID, int NGridID) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(!OverlapExists(MGridID, NGridID), "Overlap (%i,%i) already exists.", MGridID,
    NGridID);

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
    LocalMs_.Insert({MGridID,NGridID}, std::move(OverlapM));
  }

  if (NGridInfo.IsLocal()) {
    const grid &NGrid = Domain.Grid(NGridID);
    overlap_n OverlapN = core::CreateOverlapN(SharedContext, NGrid, MGridInfo);
    LocalNs_.Insert({MGridID,NGridID}, std::move(OverlapN));
  }

  OverlapRecords_.Insert({MGridID,NGridID});

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating overlap %s.(%s,%s).", Domain.Name(),
    MGridInfo.Name(), NGridInfo.Name());

  OverlapEvent_.Trigger(MGridID, NGridID, overlap_event_flags::CREATE, true);

  MPI_Barrier(Domain.Comm());

}

void overlap_component::CreateOverlaps(array_view<const int> MGridIDs, array_view<const
  int> NGridIDs) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int NumCreates = MGridIDs.Count();

  OVK_DEBUG_ASSERT(NGridIDs.Count() == NumCreates, "Incorrect N grid IDs array size.");

  if (OVK_DEBUG) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int MGridID = MGridIDs(iCreate);
      int NGridID = NGridIDs(iCreate);
      OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
      OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
      OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
      OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
      OVK_DEBUG_ASSERT(!OverlapExists(MGridID, NGridID), "Overlap (%i,%i) already exists.", MGridID,
        NGridID);
    }
  }

  SyncEdits_();

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int MGridID = MGridIDs(iCreate);
      int NGridID = NGridIDs(iCreate);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating overlap %s.(%s,%s)...",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  const std::shared_ptr<context> &SharedContext = Domain.SharedContext();

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int MGridID = MGridIDs(iCreate);
    int NGridID = NGridIDs(iCreate);
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    if (MGridInfo.IsLocal()) {
      const grid &MGrid = Domain.Grid(MGridID);
      overlap_m OverlapM = core::CreateOverlapM(SharedContext, MGrid, NGridInfo);
      LocalMs_.Insert({MGridID,NGridID}, std::move(OverlapM));
    }
    if (NGridInfo.IsLocal()) {
      const grid &NGrid = Domain.Grid(NGridID);
      overlap_n OverlapN = core::CreateOverlapN(SharedContext, NGrid, MGridInfo);
      LocalNs_.Insert({MGridID,NGridID}, std::move(OverlapN));
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int MGridID = MGridIDs(iCreate);
    int NGridID = NGridIDs(iCreate);
    OverlapRecords_.Insert({MGridID,NGridID});
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int MGridID = MGridIDs(iCreate);
      int NGridID = NGridIDs(iCreate);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating overlap %s.(%s,%s).",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int MGridID = MGridIDs(iCreate);
    int NGridID = NGridIDs(iCreate);
    OverlapEvent_.Trigger(MGridID, NGridID, overlap_event_flags::CREATE, iCreate == NumCreates-1);
  }

  MPI_Barrier(Domain.Comm());

}

void overlap_component::DestroyOverlap(int MGridID, int NGridID) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
    NGridID);

  if (OVK_DEBUG) {
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    int Editing = 0;
    if (MGridInfo.IsLocal()) {
      Editing = Editing || LocalMs_(MGridID,NGridID).Editor.Active();
    }
    if (NGridInfo.IsLocal()) {
      Editing = Editing || LocalNs_(MGridID,NGridID).Editor.Active();
    }
    MPI_Allreduce(MPI_IN_PLACE, &Editing, 1, MPI_INT, MPI_LOR, Domain.Comm());
    OVK_DEBUG_ASSERT(!Editing, "Cannot destroy overlap %s.(%s,%s); still being edited.",
      Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
  }

  SyncEdits_();

  OverlapEvent_.Trigger(MGridID, NGridID, overlap_event_flags::DESTROY, true);

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  const grid_info &MGridInfo = Domain.GridInfo(MGridID);
  const grid_info &NGridInfo = Domain.GridInfo(NGridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying overlap %s.(%s,%s)...", Domain.Name(),
    MGridInfo.Name(), NGridInfo.Name());

  LocalMs_.Erase(MGridID,NGridID);
  LocalNs_.Erase(MGridID,NGridID);

  OverlapRecords_.Erase(MGridID,NGridID);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying overlap %s.(%s,%s).",
    Domain.Name(), MGridInfo.Name(), NGridInfo.Name());

}

void overlap_component::DestroyOverlaps(array_view<const int> MGridIDs, array_view<const int>
  NGridIDs) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int NumDestroys = MGridIDs.Count();

  OVK_DEBUG_ASSERT(NGridIDs.Count() == NumDestroys, "Incorrect N grid IDs array size.");

  if (OVK_DEBUG) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int MGridID = MGridIDs(iDestroy);
      int NGridID = NGridIDs(iDestroy);
      OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
      OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
      OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
      OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
      OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
        NGridID);
    }
  }

  if (OVK_DEBUG) {
    array<int> Editing({NumDestroys}, 0);
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int MGridID = MGridIDs(iDestroy);
      int NGridID = NGridIDs(iDestroy);
      if (Domain.GridIsLocal(MGridID)) {
        Editing(iDestroy) = Editing(iDestroy) || LocalMs_(MGridID,NGridID).Editor.Active();
      }
      if (Domain.GridIsLocal(NGridID)) {
        Editing(iDestroy) = Editing(iDestroy) || LocalNs_(MGridID,NGridID).Editor.Active();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, Editing.Data(), NumDestroys, MPI_INT, MPI_LOR, Domain.Comm());
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int MGridID = MGridIDs(iDestroy);
      int NGridID = NGridIDs(iDestroy);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      OVK_DEBUG_ASSERT(!Editing(iDestroy), "Cannot destroy overlap %s.(%s,%s); still being edited.",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  SyncEdits_();

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int MGridID = MGridIDs(iDestroy);
    int NGridID = NGridIDs(iDestroy);
    OverlapEvent_.Trigger(MGridID, NGridID, overlap_event_flags::DESTROY, iDestroy ==
      NumDestroys-1);
  }

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int MGridID = MGridIDs(iDestroy);
      int NGridID = NGridIDs(iDestroy);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying overlap %s.(%s,%s)...",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int MGridID = MGridIDs(iDestroy);
    int NGridID = NGridIDs(iDestroy);
    LocalMs_.Erase(MGridID,NGridID);
    LocalNs_.Erase(MGridID,NGridID);
  }

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int MGridID = MGridIDs(iDestroy);
    int NGridID = NGridIDs(iDestroy);
    OverlapRecords_.Erase(MGridID,NGridID);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int MGridID = MGridIDs(iDestroy);
      int NGridID = NGridIDs(iDestroy);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying overlap %s.(%s,%s).",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

}

void overlap_component::ClearOverlaps() {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Clearing overlaps...");

  int NumOverlaps = OverlapCount();

  if (OVK_DEBUG) {
    array<int> Editing({NumOverlaps}, 0);
    int iOverlap = 0;
    for (auto &IDPair : OverlapRecords_.Keys()) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      if (Domain.GridIsLocal(MGridID)) {
        Editing(iOverlap) = Editing(iOverlap) || LocalMs_(MGridID,NGridID).Editor.
          Active();
      }
      if (Domain.GridIsLocal(NGridID)) {
        Editing(iOverlap) = Editing(iOverlap) || LocalNs_(MGridID,NGridID).Editor.
          Active();
      }
      ++iOverlap;
    }
    MPI_Allreduce(MPI_IN_PLACE, Editing.Data(), NumOverlaps, MPI_INT, MPI_LOR, Domain.Comm());
    iOverlap = 0;
    for (auto &IDPair : OverlapRecords_.Keys()) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      OVK_DEBUG_ASSERT(!Editing(iOverlap), "Cannot destroy overlap %s.(%s,%s); still being edited.",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
      ++iOverlap;
    }
  }

  SyncEdits_();

  int iOverlap = 0;
  for (auto &IDPair : OverlapRecords_.Keys()) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    OverlapEvent_.Trigger(MGridID, NGridID, overlap_event_flags::DESTROY, iOverlap
      == NumOverlaps-1);
    ++iOverlap;
  }

  MPI_Barrier(Domain.Comm());

  LocalMs_.Clear();
  LocalNs_.Clear();

  OverlapRecords_.Clear();

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done clearing overlaps.");

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
    if (LocalMEntry.Key(0) == MGridID) ++NumLocalMsForGrid;
  }

  return NumLocalMsForGrid;

}

const id_set<2> &overlap_component::LocalOverlapMIDs() const {

  return LocalMs_.Keys();

}

const overlap_m &overlap_component::OverlapM(int MGridID, int NGridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
    NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  return LocalMs_(MGridID,NGridID).Overlap;

}

bool overlap_component::EditingOverlapM(int MGridID, int NGridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
    NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  const local_m &LocalM = LocalMs_(MGridID,NGridID);
  const editor &Editor = LocalM.Editor;

  return Editor.Active();

}

edit_handle<overlap_m> overlap_component::EditOverlapM(int MGridID, int NGridID) {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
    NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_m &LocalM = LocalMs_(MGridID,NGridID);
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

void overlap_component::RestoreOverlapM(int MGridID, int NGridID) {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
    NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_m &LocalM = LocalMs_(MGridID,NGridID);
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
    if (LocalNEntry.Key(1) == NGridID) ++NumLocalNsForGrid;
  }

  return NumLocalNsForGrid;

}

const id_set<2> &overlap_component::LocalOverlapNIDs() const {

  return LocalNs_.Keys();

}

const overlap_n &overlap_component::OverlapN(int MGridID, int NGridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
    NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  return LocalNs_(MGridID,NGridID).Overlap;

}

bool overlap_component::EditingOverlapN(int MGridID, int NGridID) const {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
    NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  const local_n &LocalN = LocalNs_(MGridID,NGridID);
  const editor &Editor = LocalN.Editor;

  return Editor.Active();

}

edit_handle<overlap_n> overlap_component::EditOverlapN(int MGridID, int NGridID) {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
    NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_n &LocalN = LocalNs_(MGridID,NGridID);
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

void overlap_component::RestoreOverlapN(int MGridID, int NGridID) {

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(OverlapExists(MGridID, NGridID), "Overlap (%i,%i) does not exist.", MGridID,
    NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_n &LocalN = LocalNs_(MGridID,NGridID);
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
