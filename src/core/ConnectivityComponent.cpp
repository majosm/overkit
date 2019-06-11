// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
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
#include "ovk/core/Editor.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/IDMap.hpp"
#include "ovk/core/IDSet.hpp"
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
  connectivity_component_base(Domain, std::move(*Params.Name_)),
  FloatingRefGenerator_(*this)
{

  floating_ref<connectivity_component> FloatingRef = FloatingRefGenerator_.Generate();

  GridEventListener_ = Domain.AddGridEventListener([FloatingRef](int GridID, grid_event_flags Flags,
    bool LastInSequence) {
    connectivity_component &ConnectivityComponent = *FloatingRef;
    grid_event_flags &AccumulatedFlags = ConnectivityComponent.GridEventFlags_.Get(GridID,
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

  id_set<1> DyingGridIDs;

  for (auto &EventEntry : GridEventFlags_) {
    int GridID = EventEntry.Key(0);
    grid_event_flags EventFlags = EventEntry.Value();
    if ((EventFlags & grid_event_flags::DESTROY) != grid_event_flags::NONE) {
      DyingGridIDs.Insert(GridID);
    }
  }

  id_set<2> DyingConnectivityIDs;

  for (auto &IDPair : ConnectivityRecords_.Keys()) {
    if (DyingGridIDs.Contains(IDPair(0)) || DyingGridIDs.Contains(IDPair(1))) {
      DyingConnectivityIDs.Insert(IDPair);
    }
  }

  if (DyingConnectivityIDs.Count() > 0) {

    array<int> MGridIDs, NGridIDs;

    MGridIDs.Reserve(DyingConnectivityIDs.Count());
    NGridIDs.Reserve(DyingConnectivityIDs.Count());

    for (auto &IDPair : DyingConnectivityIDs) {
      MGridIDs.Append(IDPair(0));
      NGridIDs.Append(IDPair(1));
    }

    DestroyConnectivities(MGridIDs, NGridIDs);

  }

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

  id_map<2,int> GridIDsToIndex;

  int NextIndex = 0;
  for (auto &IDPair : ConnectivityRecords_.Keys()) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    GridIDsToIndex.Insert({MGridID,NGridID}, NextIndex);
    ++NextIndex;
  }

  array<connectivity_event_flags> AllConnectivityEventFlags({NumConnectivities},
    connectivity_event_flags::NONE);

  for (auto &LocalMEntry : LocalMs_) {
    int MGridID = LocalMEntry.Key(0);
    int NGridID = LocalMEntry.Key(1);
    local_m &LocalM = LocalMEntry.Value();
    int iConnectivity = GridIDsToIndex(MGridID,NGridID);
    AllConnectivityEventFlags(iConnectivity) |= LocalM.EventFlags;
  }

  for (auto &LocalNEntry : LocalNs_) {
    int MGridID = LocalNEntry.Key(0);
    int NGridID = LocalNEntry.Key(1);
    local_n &LocalN = LocalNEntry.Value();
    int iConnectivity = GridIDsToIndex(MGridID,NGridID);
    AllConnectivityEventFlags(iConnectivity) |= LocalN.EventFlags;
  }

  MPI_Allreduce(MPI_IN_PLACE, AllConnectivityEventFlags.Data(), NumConnectivities, MPI_INT, MPI_BOR,
    Domain.Comm());

  int NumTriggers = 0;
  for (auto EventFlags : AllConnectivityEventFlags) {
    if (EventFlags != connectivity_event_flags::NONE) ++NumTriggers;
  }

  int iTrigger = 0;
  for (auto &IDPair : ConnectivityRecords_.Keys()) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    int iConnectivity = GridIDsToIndex(MGridID,NGridID);
    connectivity_event_flags EventFlags = AllConnectivityEventFlags(iConnectivity);
    if (EventFlags != connectivity_event_flags::NONE) {
      ConnectivityEvent_.Trigger(MGridID, NGridID, EventFlags, iTrigger == NumTriggers-1);
      ++iTrigger;
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

int connectivity_component::ConnectivityCount() const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  return ConnectivityRecords_.Count();

}

const id_set<2> &connectivity_component::ConnectivityIDs() const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  return ConnectivityRecords_.Keys();

}

bool connectivity_component::ConnectivityExists(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);

  return ConnectivityRecords_.Contains(MGridID,NGridID);

}

void connectivity_component::CreateConnectivity(int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(!ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) already exists.",
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
    connectivity_m ConnectivityM = core::CreateConnectivityM(SharedContext, MGridID, MGrid, NGridID,
      NGridInfo);
    LocalMs_.Insert({MGridID,NGridID}, std::move(ConnectivityM));
  }

  if (NGridInfo.IsLocal()) {
    const grid &NGrid = Domain.Grid(NGridID);
    connectivity_n ConnectivityN = core::CreateConnectivityN(SharedContext, NGridID, NGrid, MGridID,
      MGridInfo);
    LocalNs_.Insert({MGridID,NGridID}, std::move(ConnectivityN));
  }

  ConnectivityRecords_.Insert({MGridID,NGridID});

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating connectivity %s.(%s,%s).",
    Domain.Name(), MGridInfo.Name(), NGridInfo.Name());

  ConnectivityEvent_.Trigger(MGridID, NGridID, connectivity_event_flags::CREATE, true);

  MPI_Barrier(Domain.Comm());

}

void connectivity_component::CreateConnectivities(array_view<const int> MGridIDs, array_view<const
  int> NGridIDs) {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

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
      OVK_DEBUG_ASSERT(!ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) already "
        "exists.", MGridID, NGridID);
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
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Creating connectivity %s.(%s,%s)...",
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
      connectivity_m ConnectivityM = core::CreateConnectivityM(SharedContext, MGridID, MGrid,
        NGridID, NGridInfo);
      LocalMs_.Insert({MGridID,NGridID}, std::move(ConnectivityM));
    }
    if (NGridInfo.IsLocal()) {
      const grid &NGrid = Domain.Grid(NGridID);
      connectivity_n ConnectivityN = core::CreateConnectivityN(SharedContext, NGridID, NGrid,
        MGridID, MGridInfo);
      LocalNs_.Insert({MGridID,NGridID}, std::move(ConnectivityN));
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int MGridID = MGridIDs(iCreate);
    int NGridID = NGridIDs(iCreate);
    ConnectivityRecords_.Insert({MGridID,NGridID});
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
      int MGridID = MGridIDs(iCreate);
      int NGridID = NGridIDs(iCreate);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done creating connectivity %s.(%s,%s).",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

  for (int iCreate = 0; iCreate < NumCreates; ++iCreate) {
    int MGridID = MGridIDs(iCreate);
    int NGridID = NGridIDs(iCreate);
    ConnectivityEvent_.Trigger(MGridID, NGridID, connectivity_event_flags::CREATE, iCreate ==
      NumCreates-1);
  }

  MPI_Barrier(Domain.Comm());

}

void connectivity_component::DestroyConnectivity(int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);

  if (OVK_DEBUG) {
    int Editing = 0;
    if (Domain.GridIsLocal(MGridID)) {
      Editing = Editing || LocalMs_(MGridID,NGridID).Editor.Active();
    }
    if (Domain.GridIsLocal(NGridID)) {
      Editing = Editing || LocalNs_(MGridID,NGridID).Editor.Active();
    }
    MPI_Allreduce(MPI_IN_PLACE, &Editing, 1, MPI_INT, MPI_LOR, Domain.Comm());
    OVK_DEBUG_ASSERT(!Editing, "Cannot destroy connectivity (%i,%i); still being edited.", MGridID,
      NGridID);
  }

  SyncEdits_();

  ConnectivityEvent_.Trigger(MGridID, NGridID, connectivity_event_flags::DESTROY, true);

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  const grid_info &MGridInfo = Domain.GridInfo(MGridID);
  const grid_info &NGridInfo = Domain.GridInfo(NGridID);

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying connectivity %s.(%s,%s)...",
    Domain.Name(), MGridInfo.Name(), NGridInfo.Name());

  LocalMs_.Erase(MGridID,NGridID);
  LocalNs_.Erase(MGridID,NGridID);

  ConnectivityRecords_.Erase(MGridID,NGridID);

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying connectivity %s.(%s,%s).",
    Domain.Name(), MGridInfo.Name(), NGridInfo.Name());

}

void connectivity_component::DestroyConnectivities(array_view<const int> MGridIDs, array_view<const
  int> NGridIDs) {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

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
      OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
        MGridID, NGridID);
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
      OVK_DEBUG_ASSERT(!Editing(iDestroy), "Cannot destroy connectivity (%i,%i); still being "
        "edited.", MGridID, NGridID);
    }
  }

  SyncEdits_();

  for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
    int MGridID = MGridIDs(iDestroy);
    int NGridID = NGridIDs(iDestroy);
    ConnectivityEvent_.Trigger(MGridID, NGridID, connectivity_event_flags::DESTROY, iDestroy ==
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
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Destroying connectivity %s.(%s,%s)...",
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
    ConnectivityRecords_.Erase(MGridID,NGridID);
  }

  MPI_Barrier(Domain.Comm());

  if (Logger.LoggingStatus()) {
    for (int iDestroy = 0; iDestroy < NumDestroys; ++iDestroy) {
      int MGridID = MGridIDs(iDestroy);
      int NGridID = NGridIDs(iDestroy);
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done destroying connectivity %s.(%s,%s).",
        Domain.Name(), MGridInfo.Name(), NGridInfo.Name());
    }
  }

}

void connectivity_component::ClearConnectivities() {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Clearing connectivities...");

  int NumConnectivities = ConnectivityCount();

  if (OVK_DEBUG) {
    array<int> Editing({NumConnectivities}, 0);
    int iConnectivity = 0;
    for (auto &IDPair : ConnectivityRecords_.Keys()) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      if (Domain.GridIsLocal(MGridID)) {
        Editing(iConnectivity) = Editing(iConnectivity) || LocalMs_(MGridID,NGridID).Editor.
          Active();
      }
      if (Domain.GridIsLocal(NGridID)) {
        Editing(iConnectivity) = Editing(iConnectivity) || LocalNs_(MGridID,NGridID).Editor.
          Active();
      }
      ++iConnectivity;
    }
    MPI_Allreduce(MPI_IN_PLACE, Editing.Data(), NumConnectivities, MPI_INT, MPI_LOR, Domain.Comm());
    iConnectivity = 0;
    for (auto &IDPair : ConnectivityRecords_.Keys()) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      OVK_DEBUG_ASSERT(!Editing(iConnectivity), "Cannot destroy connectivity (%i,%i); still being "
        "edited.", MGridID, NGridID);
      ++iConnectivity;
    }
  }

  SyncEdits_();

  int iConnectivity = 0;
  for (auto &IDPair : ConnectivityRecords_.Keys()) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    ConnectivityEvent_.Trigger(MGridID, NGridID, connectivity_event_flags::DESTROY, iConnectivity
      == NumConnectivities-1);
    ++iConnectivity;
  }

  MPI_Barrier(Domain.Comm());

  LocalMs_.Clear();
  LocalNs_.Clear();

  ConnectivityRecords_.Clear();

  MPI_Barrier(Domain.Comm());

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done clearing connectivities.");

}

int connectivity_component::LocalConnectivityMCount() const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  return LocalMs_.Count();

}

int connectivity_component::LocalConnectivityMCountForGrid(int MGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);

  int NumLocalMsForGrid = 0;

  for (auto &LocalMEntry : LocalMs_) {
    if (LocalMEntry.Key(0) == MGridID) ++NumLocalMsForGrid;
  }

  return NumLocalMsForGrid;

}

const id_set<2> &connectivity_component::LocalConnectivityMIDs() const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  return LocalMs_.Keys();

}

const connectivity_m &connectivity_component::ConnectivityM(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &Info = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(Info.IsLocal(), "M grid %s is not local to rank @rank@.", Info.Name());
  }

  return LocalMs_(MGridID,NGridID).Connectivity;

}

bool connectivity_component::EditingConnectivityM(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &Info = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(Info.IsLocal(), "M grid %s is not local to rank @rank@.", Info.Name());
  }

  const local_m &LocalM = LocalMs_(MGridID,NGridID);
  const editor &Editor = LocalM.Editor;

  return Editor.Active();

}

edit_handle<connectivity_m> connectivity_component::EditConnectivityM(int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &Info = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(Info.IsLocal(), "M grid %s is not local to rank @rank@.", Info.Name());
  }

  local_m &LocalM = LocalMs_(MGridID,NGridID);
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

void connectivity_component::RestoreConnectivityM(int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &Info = Domain.GridInfo(MGridID);
    OVK_DEBUG_ASSERT(Info.IsLocal(), "M grid %s is not local to rank @rank@.", Info.Name());
  }

  local_m &LocalM = LocalMs_(MGridID,NGridID);
  editor &Editor = LocalM.Editor;

  OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore connectivity M (%i,%i); not currently "
    "being edited.", MGridID, NGridID);

  Editor.Restore();

}

int connectivity_component::LocalConnectivityNCount() const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  return LocalNs_.Count();

}

int connectivity_component::LocalConnectivityNCountForGrid(int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);

  int NumLocalNsForGrid = 0;

  for (auto &LocalNEntry : LocalNs_) {
    if (LocalNEntry.Key(1) == NGridID) ++NumLocalNsForGrid;
  }

  return NumLocalNsForGrid;

}

const id_set<2> &connectivity_component::LocalConnectivityNIDs() const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  return LocalNs_.Keys();

}

const connectivity_n &connectivity_component::ConnectivityN(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &Info = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(Info.IsLocal(), "N grid %s is not local to rank @rank@.", Info.Name());
  }

  return LocalNs_(MGridID,NGridID).Connectivity;

}

bool connectivity_component::EditingConnectivityN(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &Info = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(Info.IsLocal(), "N grid %s is not local to rank @rank@.", Info.Name());
  }

  const local_n &LocalN = LocalNs_(MGridID,NGridID);
  const editor &Editor = LocalN.Editor;

  return Editor.Active();

}

edit_handle<connectivity_n> connectivity_component::EditConnectivityN(int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &Info = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(Info.IsLocal(), "N grid %s is not local to rank @rank@.", Info.Name());
  }

  local_n &LocalN = LocalNs_(MGridID,NGridID);
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

void connectivity_component::RestoreConnectivityN(int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(Domain_, "Connectivity component is not bound to a domain.");

  const core::domain_base &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Domain.GridExists(MGridID), "Grid %i does not exist.", MGridID);
  OVK_DEBUG_ASSERT(Domain.GridExists(NGridID), "Grid %i does not exist.", NGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(MGridID, NGridID), "Connectivity (%i,%i) does not exist.",
    MGridID, NGridID);
  if (OVK_DEBUG) {
    const grid_info &Info = Domain.GridInfo(NGridID);
    OVK_DEBUG_ASSERT(Info.IsLocal(), "N grid %s is not local to rank @rank@.", Info.Name());
  }

  local_n &LocalN = LocalNs_(MGridID,NGridID);
  editor &Editor = LocalN.Editor;

  OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore connectivity N (%i,%i); not currently "
    "being edited.", MGridID, NGridID);

  Editor.Restore();

}

connectivity_component::local_m::local_m(connectivity_m Connectivity_):
  Connectivity(std::move(Connectivity_)),
  EventFlags(connectivity_event_flags::NONE),
  EventFlagsFloatingRefGenerator(EventFlags)
{

  floating_ref<connectivity_event_flags> EventFlagsRef = EventFlagsFloatingRefGenerator.Generate();

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
  EventFlags(connectivity_event_flags::NONE),
  EventFlagsFloatingRefGenerator(EventFlags)
{

  floating_ref<connectivity_event_flags> EventFlagsRef = EventFlagsFloatingRefGenerator.Generate();

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
