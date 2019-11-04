// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ConnectivityComponent.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/ConnectivityM.hpp"
#include "ovk/core/ConnectivityN.hpp"
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
#include "ovk/core/Set.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <algorithm>
#include <string>
#include <utility>

namespace ovk {

namespace connectivity_component_internal {

connectivity_component_base::connectivity_component_base(const core::domain_base &Domain,
  std::string &&Name):
  Context_(Domain.Context().GetFloatingRef()),
  Domain_(Domain.GetFloatingRef()),
  Name_(std::move(Name))
{
  MPI_Barrier(Domain.Comm());
}

connectivity_component_base::~connectivity_component_base() noexcept {

  if (Context_) {
    const core::domain_base &Domain = *Domain_;
    MPI_Barrier(Domain.Comm());
    core::logger &Logger = Context_->core_Logger();
    Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroyed connectivity component %s.%s.",
      Domain.Name(), *Name_);
  }

}

}

connectivity_component::connectivity_component(const core::domain_base &Domain, params Params):
  connectivity_component_base(Domain, std::move(*Params.Name_))
{

  floating_ref<connectivity_component> FloatingRef = FloatingRefGenerator_.Generate(*this);

  GridEventListener_ = Domain.AddGridEventListener([FloatingRef](int GridID, grid_event_flags Flags,
    bool LastInSequence) {
    connectivity_component &ConnectivityComponent = *FloatingRef;
    grid_event_flags &AccumulatedFlags = ConnectivityComponent.GridEventFlags_.Fetch(GridID,
      grid_event_flags::NONE);
    AccumulatedFlags |= Flags;
    if (LastInSequence) {
      ConnectivityComponent.OnGridEvent_();
    }
  });

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Created connectivity component %s.%s.",
    Domain.Name(), *Name_);

}

void connectivity_component::OnGridEvent_() {

  DestroyConnectivitiesForDyingGrids_();

  GridEventFlags_.Clear();

}

void connectivity_component::DestroyConnectivitiesForDyingGrids_() {

  set<int> DyingGridIDs;

  for (auto &EventEntry : GridEventFlags_) {
    int GridID = EventEntry.Key();
    grid_event_flags EventFlags = EventEntry.Value();
    if ((EventFlags & grid_event_flags::DESTROY) != grid_event_flags::NONE) {
      DyingGridIDs.Insert(GridID);
    }
  }

  elem_set<int,2> DyingConnectivityIDs;

  for (auto &ConnectivityID : ConnectivityRecords_.Keys()) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    if (DyingGridIDs.Contains(MGridID) || DyingGridIDs.Contains(NGridID)) {
      DyingConnectivityIDs.Insert(ConnectivityID);
    }
  }

  DestroyConnectivities(DyingConnectivityIDs);

}

void connectivity_component::StartEdit() {

  // Nothing to do here

}

void connectivity_component::EndEdit() {

  SyncEdits_();

}

void connectivity_component::SyncEdits_() {

  const core::domain_base &Domain = *Domain_;

  int HasEdits = false;
  for (auto &LocalMEntry : LocalMs_) {
    local_m &LocalM = LocalMEntry.Value();
    if (LocalM.EventFlags != connectivity_event_flags::NONE) {
      HasEdits = true;
      goto done_looping;
    }
  }
  for (auto &LocalNEntry : LocalNs_) {
    local_n &LocalN = LocalNEntry.Value();
    if (LocalN.EventFlags != connectivity_event_flags::NONE) {
      HasEdits = true;
      goto done_looping;
    }
  }
  done_looping:;
  MPI_Allreduce(MPI_IN_PLACE, &HasEdits, 1, MPI_INT, MPI_MAX, Domain.Comm());

  if (!HasEdits) return;

  int NumConnectivities = ConnectivityRecords_.Count();

  elem_map<int,2,int> GridIDsToIndex;

  int NextIndex = 0;
  for (auto &ConnectivityID : ConnectivityRecords_.Keys()) {
    GridIDsToIndex.Insert(ConnectivityID, NextIndex);
    ++NextIndex;
  }

  array<connectivity_event_flags> AllConnectivityEventFlags({NumConnectivities},
    connectivity_event_flags::NONE);

  for (auto &LocalMEntry : LocalMs_) {
    local_m &LocalM = LocalMEntry.Value();
    int iConnectivity = GridIDsToIndex(LocalMEntry.Key());
    AllConnectivityEventFlags(iConnectivity) |= LocalM.EventFlags;
  }

  for (auto &LocalNEntry : LocalNs_) {
    local_n &LocalN = LocalNEntry.Value();
    int iConnectivity = GridIDsToIndex(LocalNEntry.Key());
    AllConnectivityEventFlags(iConnectivity) |= LocalN.EventFlags;
  }

  MPI_Allreduce(MPI_IN_PLACE, AllConnectivityEventFlags.Data(), NumConnectivities, MPI_INT, MPI_BOR,
    Domain.Comm());

  int NumTriggers = 0;
  for (auto EventFlags : AllConnectivityEventFlags) {
    if (EventFlags != connectivity_event_flags::NONE) ++NumTriggers;
  }
  for (auto &ConnectivityID : ConnectivityRecords_.Keys()) {
    int iConnectivity = GridIDsToIndex(ConnectivityID);
    connectivity_event_flags EventFlags = AllConnectivityEventFlags(iConnectivity);
    if (EventFlags != connectivity_event_flags::NONE) {
      ConnectivityEvent_.Trigger(ConnectivityID, EventFlags, --NumTriggers == 0);
    }
  }

  MPI_Barrier(Domain.Comm());

  for (auto &LocalMEntry : LocalMs_) {
    local_m &LocalM = LocalMEntry.Value();
    LocalM.EventFlags = connectivity_event_flags::NONE;
  }

  for (auto &LocalNEntry : LocalNs_) {
    local_n &LocalN = LocalNEntry.Value();
    LocalN.EventFlags = connectivity_event_flags::NONE;
  }

}

bool connectivity_component::ConnectivityExists(const elem<int,2> &ConnectivityID) const {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  }

  return ConnectivityRecords_.Contains(ConnectivityID);

}

void connectivity_component::CreateConnectivity(const elem<int,2> &ConnectivityID) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int MGridID = ConnectivityID(0);
  int NGridID = ConnectivityID(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(!ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) already exists.",
    MGridID, NGridID);

  SyncEdits_();

  core::logger &Logger = Context_->core_Logger();

  const grid_info &MGridInfo = Domain.GridInfo(MGridID);
  const grid_info &NGridInfo = Domain.GridInfo(NGridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating connectivity %s.(%s,%s)...",
    Domain.Name(), MGridInfo.Name(), NGridInfo.Name());

  const std::shared_ptr<context> &SharedContext = Domain.SharedContext();

  if (MGridInfo.IsLocal()) {
    const grid &MGrid = Domain.Grid(MGridID);
    connectivity_m ConnectivityM = core::CreateConnectivityM(SharedContext, MGrid, NGridInfo);
    LocalMs_.Insert(ConnectivityID, std::move(ConnectivityM));
  }

  if (NGridInfo.IsLocal()) {
    const grid &NGrid = Domain.Grid(NGridID);
    connectivity_n ConnectivityN = core::CreateConnectivityN(SharedContext, NGrid, MGridInfo);
    LocalNs_.Insert(ConnectivityID, std::move(ConnectivityN));
  }

  ConnectivityRecords_.Insert(ConnectivityID);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating connectivity %s.(%s,%s).",
    Domain.Name(), MGridInfo.Name(), NGridInfo.Name());

  ConnectivityEvent_.Trigger(ConnectivityID, connectivity_event_flags::CREATE, true);

  MPI_Barrier(Domain.Comm());

}

void connectivity_component::CreateConnectivities(array_view<const elem<int,2>> ConnectivityIDs) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  if (OVK_DEBUG) {
    for (auto &ConnectivityID : ConnectivityIDs) {
      int MGridID = ConnectivityID(0);
      int NGridID = ConnectivityID(1);
      OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
      OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
      OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
      OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
      OVK_DEBUG_ASSERT(!ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) already exists.",
        MGridID, NGridID);
    }
  }

  SyncEdits_();

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (auto &ConnectivityID : ConnectivityIDs) {
      int MGridID = ConnectivityID(0);
      int NGridID = ConnectivityID(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating connectivity %s.(%s,%s)...",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  const std::shared_ptr<context> &SharedContext = Domain.SharedContext();

  for (auto &ConnectivityID : ConnectivityIDs) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    if (MGridInfo.IsLocal()) {
      const grid &MGrid = Domain.Grid(MGridID);
      connectivity_m ConnectivityM = core::CreateConnectivityM(SharedContext, MGrid, NGridInfo);
      LocalMs_.Insert(ConnectivityID, std::move(ConnectivityM));
    }
    if (NGridInfo.IsLocal()) {
      const grid &NGrid = Domain.Grid(NGridID);
      connectivity_n ConnectivityN = core::CreateConnectivityN(SharedContext, NGrid, MGridInfo);
      LocalNs_.Insert(ConnectivityID, std::move(ConnectivityN));
    }
  }

  for (auto &ConnectivityID : ConnectivityIDs) {
    ConnectivityRecords_.Insert(ConnectivityID);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (auto &ConnectivityID : ConnectivityIDs) {
      int MGridID = ConnectivityID(0);
      int NGridID = ConnectivityID(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating connectivity %s.(%s,%s).",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  int NumRemaining = ConnectivityIDs.Count();
  for (auto &ConnectivityID : ConnectivityIDs) {
    ConnectivityEvent_.Trigger(ConnectivityID, connectivity_event_flags::CREATE,
      --NumRemaining == 0);
  }

  MPI_Barrier(Domain.Comm());

}

void connectivity_component::DestroyConnectivity(const elem<int,2> &ConnectivityID) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int MGridID = ConnectivityID(0);
  int NGridID = ConnectivityID(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);

  if (OVK_DEBUG) {
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    int Editing = 0;
    if (MGridInfo.IsLocal()) {
      Editing = Editing || LocalMs_(ConnectivityID).Editor.Active();
    }
    if (NGridInfo.IsLocal()) {
      Editing = Editing || LocalNs_(ConnectivityID).Editor.Active();
    }
    MPI_Allreduce(MPI_IN_PLACE, &Editing, 1, MPI_INT, MPI_LOR, Domain.Comm());
    OVK_DEBUG_ASSERT(!Editing, "Cannot destroy connectivity %s.(%s,%s); still being edited.",
      Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
  }

  SyncEdits_();

  ConnectivityEvent_.Trigger(ConnectivityID, connectivity_event_flags::DESTROY, true);

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  const grid_info &MGridInfo = Domain.GridInfo(MGridID);
  const grid_info &NGridInfo = Domain.GridInfo(NGridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying connectivity %s.(%s,%s)...",
    Domain.Name(), MGridInfo.Name(), NGridInfo.Name());

  LocalMs_.Erase(ConnectivityID);
  LocalNs_.Erase(ConnectivityID);

  ConnectivityRecords_.Erase(ConnectivityID);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying connectivity %s.(%s,%s).",
    Domain.Name(), MGridInfo.Name(), NGridInfo.Name());

}

void connectivity_component::DestroyConnectivities(array_view<const elem<int,2>> ConnectivityIDs) {

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  int NumDestroys = ConnectivityIDs.Count();

  if (OVK_DEBUG) {
    for (auto &ConnectivityID : ConnectivityIDs) {
      int MGridID = ConnectivityID(0);
      int NGridID = ConnectivityID(1);
      OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
      OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
      OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
      OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
      OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
        MGridID, NGridID);
    }
  }

  if (OVK_DEBUG) {
    array<int> Editing({NumDestroys}, 0);
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      auto &ConnectivityID = ConnectivityIDs(iDestroy);
      int MGridID = ConnectivityID(0);
      int NGridID = ConnectivityID(1);
      if (Domain.GridIsLocal(MGridID)) {
        Editing(iDestroy) = Editing(iDestroy) || LocalMs_(ConnectivityID).Editor.Active();
      }
      if (Domain.GridIsLocal(NGridID)) {
        Editing(iDestroy) = Editing(iDestroy) || LocalNs_(ConnectivityID).Editor.Active();
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, Editing.Data(), NumDestroys, MPI_INT, MPI_LOR, Domain.Comm());
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      auto &ConnectivityID = ConnectivityIDs(iDestroy);
      int MGridID = ConnectivityID(0);
      int NGridID = ConnectivityID(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      OVK_DEBUG_ASSERT(!Editing(iDestroy), "Cannot destroy connectivity %s.(%s,%s); still being "
        "edited.", Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  SyncEdits_();

  int NumRemaining = NumDestroys;
  for (auto &ConnectivityID : ConnectivityIDs) {
    ConnectivityEvent_.Trigger(ConnectivityID, connectivity_event_flags::DESTROY,
      --NumRemaining == 0);
  }

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  if (Logger.LoggingStatus()) {
    for (auto &ConnectivityID : ConnectivityIDs) {
      int MGridID = ConnectivityID(0);
      int NGridID = ConnectivityID(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying connectivity %s.(%s,%s)...",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  for (auto &ConnectivityID : ConnectivityIDs) {
    LocalMs_.Erase(ConnectivityID);
    LocalNs_.Erase(ConnectivityID);
  }

  for (auto &ConnectivityID : ConnectivityIDs) {
    ConnectivityRecords_.Erase(ConnectivityID);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (auto &ConnectivityID : ConnectivityIDs) {
      int MGridID = ConnectivityID(0);
      int NGridID = ConnectivityID(1);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying connectivity %s.(%s,%s).",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

}

void connectivity_component::ClearConnectivities() {

  DestroyConnectivities(ConnectivityRecords_.Keys());

}

const connectivity_m &connectivity_component::ConnectivityM(const elem<int,2> &ConnectivityID) const
  {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
    OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
      MGridID, NGridID);
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  return LocalMs_(ConnectivityID).Connectivity;

}

bool connectivity_component::EditingConnectivityM(const elem<int,2> &ConnectivityID) const {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
    OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
      MGridID, NGridID);
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  const local_m &LocalM = LocalMs_(ConnectivityID);
  const editor &Editor = LocalM.Editor;

  return Editor.Active();

}

edit_handle<connectivity_m> connectivity_component::EditConnectivityM(const elem<int,2>
  &ConnectivityID) {

  const core::domain_base &Domain = *Domain_;

  int MGridID = ConnectivityID(0);
  int NGridID = ConnectivityID(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_m &LocalM = LocalMs_(ConnectivityID);
  connectivity_m &ConnectivityM = LocalM.Connectivity;
  editor &Editor = LocalM.Editor;

  if (!Editor.Active()) {
    floating_ref<const grid> GridRef = Domain.Grid(MGridID).GetFloatingRef();
    MPI_Barrier(GridRef->Comm());
    auto DeactivateFunc = [GridRef] { MPI_Barrier(GridRef->Comm()); };
    Editor.Activate(std::move(DeactivateFunc));
  }

  return Editor.Edit(ConnectivityM);

}

void connectivity_component::RestoreConnectivityM(const elem<int,2> &ConnectivityID) {

  const core::domain_base &Domain = *Domain_;

  int MGridID = ConnectivityID(0);
  int NGridID = ConnectivityID(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "M grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_m &LocalM = LocalMs_(ConnectivityID);
  editor &Editor = LocalM.Editor;

  if (OVK_DEBUG) {
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore connectivity M %s.(%s,%s); not currently "
      "being edited.", Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
  }

  Editor.Restore();

}

const connectivity_n &connectivity_component::ConnectivityN(const elem<int,2> &ConnectivityID) const
  {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
    OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
      MGridID, NGridID);
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  return LocalNs_(ConnectivityID).Connectivity;

}

bool connectivity_component::EditingConnectivityN(const elem<int,2> &ConnectivityID) const {

  const core::domain_base &Domain = *Domain_;

  if (OVK_DEBUG) {
    int MGridID = ConnectivityID(0);
    int NGridID = ConnectivityID(1);
    OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
    OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
    OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
    OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
    OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
      MGridID, NGridID);
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  const local_n &LocalN = LocalNs_(ConnectivityID);
  const editor &Editor = LocalN.Editor;

  return Editor.Active();

}

edit_handle<connectivity_n> connectivity_component::EditConnectivityN(const elem<int,2>
  &ConnectivityID) {

  const core::domain_base &Domain = *Domain_;

  int MGridID = ConnectivityID(0);
  int NGridID = ConnectivityID(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_n &LocalN = LocalNs_(ConnectivityID);
  connectivity_n &ConnectivityN = LocalN.Connectivity;
  editor &Editor = LocalN.Editor;

  if (!Editor.Active()) {
    floating_ref<const grid> GridRef = Domain.Grid(NGridID).GetFloatingRef();
    MPI_Barrier(GridRef->Comm());
    auto DeactivateFunc = [GridRef] { MPI_Barrier(GridRef->Comm()); };
    Editor.Activate(std::move(DeactivateFunc));
  }

  return Editor.Edit(ConnectivityN);

}

void connectivity_component::RestoreConnectivityN(const elem<int,2> &ConnectivityID) {

  const core::domain_base &Domain = *Domain_;

  int MGridID = ConnectivityID(0);
  int NGridID = ConnectivityID(1);

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(ConnectivityID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(GridInfo.IsLocal(), "N grid %s is not local to rank @rank@.", GridInfo.Name());
  }

  local_n &LocalN = LocalNs_(ConnectivityID);
  editor &Editor = LocalN.Editor;

  if (OVK_DEBUG) {
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore connectivity N %s.(%s,%s); not currently "
      "being edited.", Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
  }

  Editor.Restore();

}

connectivity_component::local_m::local_m(connectivity_m Connectivity_):
  Connectivity(std::move(Connectivity_)),
  EventFlags(connectivity_event_flags::NONE)
{

  floating_ref<connectivity_event_flags> EventFlagsRef = FloatingRefGenerator.Generate(EventFlags);

  ResizeEventListener = Connectivity.AddResizeEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::RESIZE_M;
  });

  ExtentsEventListener = Connectivity.AddExtentsEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::EDIT_M_EXTENTS;
  });

  CoordsEventListener = Connectivity.AddCoordsEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::EDIT_M_COORDS;
  });

  InterpCoefsEventListener = Connectivity.AddInterpCoefsEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::EDIT_M_INTERP_COEFS;
  });

  DestinationsEventListener = Connectivity.AddDestinationsEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::EDIT_M_DESTINATIONS;
  });

  DestinationRanksEventListener = Connectivity.AddDestinationRanksEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::EDIT_M_DESTINATIONS;
  });

}

connectivity_component::local_n::local_n(connectivity_n Connectivity_):
  Connectivity(std::move(Connectivity_)),
  EventFlags(connectivity_event_flags::NONE)
{

  floating_ref<connectivity_event_flags> EventFlagsRef = FloatingRefGenerator.Generate(EventFlags);

  ResizeEventListener = Connectivity.AddResizeEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::RESIZE_N;
  });

  PointsEventListener = Connectivity.AddPointsEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::EDIT_N_POINTS;
  });

  SourcesEventListener = Connectivity.AddSourcesEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::EDIT_N_SOURCES;
  });

  SourceRanksEventListener = Connectivity.AddSourceRanksEventListener([EventFlagsRef] {
    *EventFlagsRef |= connectivity_event_flags::EDIT_N_SOURCES;
  });

}

connectivity_component::params &connectivity_component::params::SetName(std::string Name) {

  Name_ = std::move(Name);

  return *this;

}

}
