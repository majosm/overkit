// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Exchanger.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Collect.hpp"
#include "ovk/core/CollectMap.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/ConnectivityComponent.hpp"
#include "ovk/core/ConnectivityM.hpp"
#include "ovk/core/ConnectivityN.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Disperse.hpp"
#include "ovk/core/DisperseMap.hpp"
#include "ovk/core/ElemMap.hpp"
#include "ovk/core/ElemSet.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Recv.hpp"
#include "ovk/core/RecvMap.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/Send.hpp"
#include "ovk/core/SendMap.hpp"
#include "ovk/core/Set.hpp"

#include <mpi.h>

#include <memory>
#include <string>
#include <utility>

namespace ovk {

namespace {

long long BinDivide(long long N, int NumBins);

array<long long> GetSendRecvOrder(const array<int,2> &ReceiverPoints, const range
  &ReceiverGridGlobalRange);

}

exchanger::exchanger(std::shared_ptr<context> &&Context, params &&Params):
  Context_(std::move(Context)),
  Name_(std::move(*Params.Name_))
{}

exchanger::~exchanger() noexcept {

  if (Domain_) {
    Unbind();
  }

}

exchanger exchanger::internal_Create(std::shared_ptr<context> &&Context, params &&Params) {

  return {std::move(Context), std::move(Params)};

}

exchanger CreateExchanger(std::shared_ptr<context> Context, exchanger::params Params) {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return exchanger::internal_Create(std::move(Context), std::move(Params));

}

const context &exchanger::Context() const {

  return *Context_;

}

context &exchanger::Context() {

  return *Context_;

}

const std::shared_ptr<context> &exchanger::SharedContext() const {

  return Context_;

}

bool exchanger::Bound() const {

  return static_cast<bool>(Domain_);

}

void exchanger::Bind(const domain &Domain, bindings Bindings) {

  OVK_DEBUG_ASSERT(!Domain_, "Exchanger is already bound to a domain.");

  Domain_ = Domain.GetFloatingRef();

  floating_ref<exchanger> FloatingRef = FloatingRefGenerator_.Generate(*this);

  ComponentEventListener_ = Domain.AddComponentEventListener([FloatingRef](
    int ComponentID, component_event_flags Flags) {
    exchanger &Exchanger = *FloatingRef;
    Exchanger.OnComponentEvent_(ComponentID, Flags);
  });

  ConnectivityComponentID_ = Bindings.ConnectivityComponentID_;
  OVK_DEBUG_ASSERT(ConnectivityComponentID_ >= 0, "Invalid connectivity component ID.");
  OVK_DEBUG_ASSERT(Domain.ComponentExists(ConnectivityComponentID_), "Component %i does not exist.",
    ConnectivityComponentID_);

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);
  ConnectivityEventListener_ = ConnectivityComponent.AddConnectivityEventListener([FloatingRef](int
    MGridID, int NGridID, connectivity_event_flags Flags, bool LastInSequence) {
    exchanger &Exchanger = *FloatingRef;
    Exchanger.OnConnectivityEvent_(MGridID, NGridID, Flags, LastInSequence);
  });

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Bound exchanger %s to domain %s.", *Name_,
    Domain.Name());

  for (auto &IDPair : ConnectivityComponent.ConnectivityIDs()) {
    UpdateManifest_.CreateLocal.Insert(IDPair);
    UpdateManifest_.UpdateSourceDestRanks.Insert(IDPair);
    UpdateManifest_.ResetExchanges.Insert(IDPair);
  }

  Update_();

}

void exchanger::Unbind() {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  LocalMs_.Clear();
  LocalNs_.Clear();

  ConnectivityEventListener_.Reset();
  ConnectivityComponentID_ = -1;

  ComponentEventListener_.Reset();
  Domain_.Reset();

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Unbound exchanger %s.", *Name_);

}

const domain &exchanger::Domain() const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  return *Domain_;

}

void exchanger::OnComponentEvent_(int ComponentID, component_event_flags Flags) {

  auto FlagsMatchAny = [&Flags](component_event_flags Mask) -> bool {
    return (Flags & Mask) != component_event_flags::NONE;
  };

  bool DestroyConnectivityComponent = ComponentID == ConnectivityComponentID_ &&
    FlagsMatchAny(component_event_flags::DESTROY);

  if (DestroyConnectivityComponent) {
    Unbind();
  }

}

void exchanger::OnConnectivityEvent_(int MGridID, int NGridID, connectivity_event_flags Flags,
  bool LastInSequence) {

  elem<int,2> IDPair = {MGridID,NGridID};

  auto FlagsMatchAny = [&Flags](connectivity_event_flags Mask) -> bool {
    return (Flags & Mask) != connectivity_event_flags::NONE;
  };

  bool Create = FlagsMatchAny(connectivity_event_flags::CREATE);
  bool Destroy = FlagsMatchAny(connectivity_event_flags::DESTROY);
  bool EditAny = FlagsMatchAny(connectivity_event_flags::ALL_EDITS);
  bool EditDests = FlagsMatchAny(connectivity_event_flags::EDIT_M_DESTINATIONS);
  bool EditSources = FlagsMatchAny(connectivity_event_flags::EDIT_N_SOURCES);

  if (Create) {
    UpdateManifest_.CreateLocal.Insert(IDPair);
  }

  if (Destroy) {
    UpdateManifest_.DestroyLocal.Insert(IDPair);
  }

  if (Create || EditDests || EditSources) {
    UpdateManifest_.UpdateSourceDestRanks.Insert(IDPair);
  }

  if (Create || EditAny) {
    UpdateManifest_.ResetExchanges.Insert(IDPair);
  }

  if (LastInSequence) {
    Update_();
  }

}

void exchanger::Update_() {

  const domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Updating exchanger %s...", *Name_);

  CreateLocals_();
  DestroyLocals_();
  UpdateSourceDestRanks_();
  ResetExchanges_();

  UpdateManifest_.CreateLocal.Clear();
  UpdateManifest_.DestroyLocal.Clear();
  UpdateManifest_.UpdateSourceDestRanks.Clear();
  UpdateManifest_.ResetExchanges.Clear();

  Logger.LogStatus(Domain.Comm().Rank() == 0, 0, "Done updating exchanger %s.", *Name_);

}

void exchanger::CreateLocals_() {

  if (UpdateManifest_.CreateLocal.Empty()) return;

  const domain &Domain = *Domain_;
  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  for (auto &IDPair : UpdateManifest_.CreateLocal) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    if (Domain.GridIsLocal(MGridID)) {
      local_m &LocalM = LocalMs_.Insert(IDPair);
      LocalM.Connectivity = &ConnectivityComponent.ConnectivityM(IDPair);
    }
    if (Domain.GridIsLocal(NGridID)) {
      local_n &LocalN = LocalNs_.Insert(IDPair);
      LocalN.Connectivity = &ConnectivityComponent.ConnectivityN(IDPair);
    }
  }

}

void exchanger::DestroyLocals_() {

  if (UpdateManifest_.DestroyLocal.Empty()) return;

  for (auto &IDPair : UpdateManifest_.DestroyLocal) {
    LocalMs_.Erase(IDPair);
    LocalNs_.Erase(IDPair);
  }

}

void exchanger::UpdateSourceDestRanks_() {

  using range_indexer = indexer<long long, int, 3, ovk::array_layout::COLUMN_MAJOR>;

  if (UpdateManifest_.UpdateSourceDestRanks.Empty()) return;

  const domain &Domain = *Domain_;
  const comm &Comm = Domain.Comm();

  MPI_Barrier(Comm);

  elem_set<int,2> ConnectivityMGridIDs;
  elem_set<int,2> ConnectivityNGridIDs;

  for (auto &IDPair : UpdateManifest_.UpdateSourceDestRanks) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    if (Domain.GridIsLocal(MGridID)) {
      ConnectivityMGridIDs.Insert(IDPair);
    }
    if (Domain.GridIsLocal(NGridID)) {
      ConnectivityNGridIDs.Insert(IDPair);
    }
  }

  for (auto &IDPair : ConnectivityMGridIDs) {
    int MGridID = IDPair(0);
    const grid &MGrid = Domain.Grid(MGridID);
    local_m &LocalM = LocalMs_(IDPair);
    const connectivity_m &ConnectivityM = *LocalM.Connectivity;
    LocalM.DestinationRanks = ConnectivityM.DestinationRanks();
    // Ensure only process containing lower corner of donor cell communicates
    const array<int,3> &Extents = ConnectivityM.Extents();
    for (long long iDonor = 0; iDonor < ConnectivityM.Count(); ++iDonor) {
      tuple<int> CellLower = {
        Extents(0,0,iDonor),
        Extents(0,1,iDonor),
        Extents(0,2,iDonor)
      };
      if (!MGrid.LocalRange().Contains(CellLower)) {
        LocalM.DestinationRanks(iDonor) = -1;
      }
    }
  }

  for (auto &IDPair : ConnectivityNGridIDs) {
    local_n &LocalN = LocalNs_(IDPair);
    const connectivity_n &ConnectivityN = *LocalN.Connectivity;
    LocalN.SourceRanks = ConnectivityN.SourceRanks();
  }

  long long TotalPoints = 0;
  for (int GridID : Domain.GridIDs()) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    TotalPoints += GridInfo.Cart().Range().Count();
  }

  struct linear_partition {
    interval<long long> Interval;
    array<int> MRanks;
    array<int> NRanks;
  };

  linear_partition LinearPartition;

  long long LinearPartitionSize = BinDivide(TotalPoints, Comm.Size());
  LinearPartition.Interval.Begin(0) = LinearPartitionSize*Comm.Rank();
  LinearPartition.Interval.End(0) = Min(LinearPartitionSize*(Comm.Rank()+1), TotalPoints);
  LinearPartition.MRanks.Resize(LinearPartition.Interval, -1);
  LinearPartition.NRanks.Resize(LinearPartition.Interval, -1);

  map<int,long long> NumPointsBeforeGrid;
  long long NumPointsPartial = 0;
  for (int GridID : Domain.GridIDs()) {
    const grid_info &GridInfo = Domain.GridInfo(GridID);
    NumPointsBeforeGrid.Insert(GridID, NumPointsPartial);
    NumPointsPartial += GridInfo.Cart().Range().Count();
  }

  struct send_recv {
    long long Count;
    array<long long> PointIndices;
    array<int> Ranks;
    send_recv():
      Count(0)
    {}
  };

  map<int,send_recv> MSends, NSends;

  for (auto &IDPair : ConnectivityMGridIDs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid &MGrid = Domain.Grid(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    range_indexer NGridGlobalIndexer(NGridInfo.Cart().Range());
    const local_m &LocalM = LocalMs_(IDPair);
    const connectivity_m &ConnectivityM = *LocalM.Connectivity;
    const array<int,3> &Extents = ConnectivityM.Extents();
    const array<int,2> &Destinations = ConnectivityM.Destinations();
    const array<int> &DestinationRanks = LocalM.DestinationRanks;
    for (long long iDonor = 0; iDonor < ConnectivityM.Count(); ++iDonor) {
      tuple<int> CellLower = {
        Extents(0,0,iDonor),
        Extents(0,1,iDonor),
        Extents(0,2,iDonor)
      };
      if (MGrid.LocalRange().Contains(CellLower) && DestinationRanks(iDonor) < 0) {
        tuple<int> Point = {
          Destinations(0,iDonor),
          Destinations(1,iDonor),
          Destinations(2,iDonor)
        };
        long long iLinearPoint = NumPointsBeforeGrid(NGridID) + NGridGlobalIndexer.ToIndex(Point);
        int iLinearPartition = int(iLinearPoint/LinearPartitionSize);
        send_recv &Send = MSends.Fetch(iLinearPartition);
        ++Send.Count;
      }
    }
  }

  for (auto &IDPair : ConnectivityNGridIDs) {
    int NGridID = IDPair(1);
    const grid &NGrid = Domain.Grid(NGridID);
    range_indexer NGridGlobalIndexer(NGrid.GlobalRange());
    const local_n &LocalN = LocalNs_(IDPair);
    const connectivity_n &ConnectivityN = *LocalN.Connectivity;
    const array<int,2> &Points = ConnectivityN.Points();
    const array<int> &SourceRanks = LocalN.SourceRanks;
    for (long long iReceiver = 0; iReceiver < ConnectivityN.Count(); ++iReceiver) {
      if (SourceRanks(iReceiver) < 0) {
        tuple<int> Point = {
          Points(0,iReceiver),
          Points(1,iReceiver),
          Points(2,iReceiver)
        };
        long long iLinearPoint = NumPointsBeforeGrid(NGridID) + NGridGlobalIndexer.ToIndex(Point);
        int iLinearPartition = int(iLinearPoint/LinearPartitionSize);
        send_recv &Send = NSends.Fetch(iLinearPartition);
        ++Send.Count;
      }
    }
  }

  for (auto &Entry : MSends) {
    send_recv &Send = Entry.Value();
    Send.PointIndices.Reserve(Send.Count);
  }

  for (auto &Entry : NSends) {
    send_recv &Send = Entry.Value();
    Send.PointIndices.Reserve(Send.Count);
  }

  for (auto &IDPair : ConnectivityMGridIDs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid &MGrid = Domain.Grid(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    range_indexer NGridGlobalIndexer(NGridInfo.Cart().Range());
    const local_m &LocalM = LocalMs_(IDPair);
    const connectivity_m &ConnectivityM = *LocalM.Connectivity;
    const array<int,3> &Extents = ConnectivityM.Extents();
    const array<int,2> &Destinations = ConnectivityM.Destinations();
    const array<int> &DestinationRanks = LocalM.DestinationRanks;
    for (long long iDonor = 0; iDonor < ConnectivityM.Count(); ++iDonor) {
      tuple<int> CellLower = {
        Extents(0,0,iDonor),
        Extents(0,1,iDonor),
        Extents(0,2,iDonor)
      };
      if (MGrid.LocalRange().Contains(CellLower) && DestinationRanks(iDonor) < 0) {
        tuple<int> Point = {
          Destinations(0,iDonor),
          Destinations(1,iDonor),
          Destinations(2,iDonor)
        };
        long long iLinearPoint = NumPointsBeforeGrid(NGridID) + NGridGlobalIndexer.ToIndex(Point);
        int iLinearPartition = int(iLinearPoint/LinearPartitionSize);
        send_recv &Send = MSends(iLinearPartition);
        Send.PointIndices.Append(iLinearPoint);
      }
    }
  }

  for (auto &IDPair : ConnectivityNGridIDs) {
    int NGridID = IDPair(1);
    const grid &NGrid = Domain.Grid(NGridID);
    range_indexer NGridGlobalIndexer(NGrid.GlobalRange());
    const local_n &LocalN = LocalNs_(IDPair);
    const connectivity_n &ConnectivityN = *LocalN.Connectivity;
    const array<int,2> &Points = ConnectivityN.Points();
    const array<int> &SourceRanks = LocalN.SourceRanks;
    for (long long iReceiver = 0; iReceiver < ConnectivityN.Count(); ++iReceiver) {
      if (SourceRanks(iReceiver) < 0) {
        tuple<int> Point = {
          Points(0,iReceiver),
          Points(1,iReceiver),
          Points(2,iReceiver)
        };
        long long iLinearPoint = NumPointsBeforeGrid(NGridID) + NGridGlobalIndexer.ToIndex(Point);
        int iLinearPartition = int(iLinearPoint/LinearPartitionSize);
        send_recv &Send = NSends(iLinearPartition);
        Send.PointIndices.Append(iLinearPoint);
      }
    }
  }

  for (auto &Entry : MSends) {
    send_recv &Send = Entry.Value();
    Send.Ranks.Resize({Send.Count});
  }

  for (auto &Entry : NSends) {
    send_recv &Send = Entry.Value();
    Send.Ranks.Resize({Send.Count});
  }

  array<int> MRecvFromRanks = core::DynamicHandshake(Comm, MSends.Keys());
  array<int> NRecvFromRanks = core::DynamicHandshake(Comm, NSends.Keys());

  map<int,send_recv> MRecvs, NRecvs;

  for (int Rank : MRecvFromRanks) {
    MRecvs.Insert(Rank);
  }

  for (int Rank : NRecvFromRanks) {
    NRecvs.Insert(Rank);
  }

  MRecvFromRanks.Clear();
  NRecvFromRanks.Clear();

  int NumMSends = MSends.Count();
  int NumNSends = NSends.Count();

  int NumMRecvs = MRecvs.Count();
  int NumNRecvs = NRecvs.Count();

  array<MPI_Request> Requests;

  Requests.Reserve(NumMSends+NumNSends+NumMRecvs+NumNRecvs);

  auto Isend = [&Requests](void *Buffer, long long Count, MPI_Datatype DataType, int DestRank,
    int Tag, MPI_Comm SendComm) {
    OVK_DEBUG_ASSERT(Count <= std::numeric_limits<int>::max(), "Send count too large.");
    MPI_Request &Request = Requests.Append();
    MPI_Isend(Buffer, int(Count), DataType, DestRank, Tag, SendComm, &Request);
  };

  auto Irecv = [&Requests](void *Buffer, long long Count, MPI_Datatype DataType, int SourceRank,
    int Tag, MPI_Comm RecvComm) {
    OVK_DEBUG_ASSERT(Count <= std::numeric_limits<int>::max(), "Receive count too large.");
    MPI_Request &Request = Requests.Append();
    MPI_Irecv(Buffer, int(Count), DataType, SourceRank, Tag, RecvComm, &Request);
  };

  for (auto &Entry : MRecvs) {
    int Rank = Entry.Key();
    send_recv &Recv = Entry.Value();
    Irecv(&Recv.Count, 1, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Entry : NRecvs) {
    int Rank = Entry.Key();
    send_recv &Recv = Entry.Value();
    Irecv(&Recv.Count, 1, MPI_LONG_LONG, Rank, 1, Comm);
  }

  for (auto &Entry : MSends) {
    int Rank = Entry.Key();
    send_recv &Send = Entry.Value();
    Isend(&Send.Count, 1, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Entry : NSends) {
    int Rank = Entry.Key();
    send_recv &Send = Entry.Value();
    Isend(&Send.Count, 1, MPI_LONG_LONG, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  for (auto &Entry : MRecvs) {
    send_recv &Recv = Entry.Value();
    Recv.PointIndices.Resize({Recv.Count});
    Recv.Ranks.Resize({Recv.Count});
  }

  for (auto &Entry : NRecvs) {
    send_recv &Recv = Entry.Value();
    Recv.PointIndices.Resize({Recv.Count});
    Recv.Ranks.Resize({Recv.Count});
  }

  for (auto &Entry : MRecvs) {
    int Rank = Entry.Key();
    send_recv &Recv = Entry.Value();
    Irecv(Recv.PointIndices.Data(), Recv.Count, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Entry : NRecvs) {
    int Rank = Entry.Key();
    send_recv &Recv = Entry.Value();
    Irecv(Recv.PointIndices.Data(), Recv.Count, MPI_LONG_LONG, Rank, 1, Comm);
  }

  for (auto &Entry : MSends) {
    int Rank = Entry.Key();
    send_recv &Send = Entry.Value();
    Isend(Send.PointIndices.Data(), Send.Count, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Entry : NSends) {
    int Rank = Entry.Key();
    send_recv &Send = Entry.Value();
    Isend(Send.PointIndices.Data(), Send.Count, MPI_LONG_LONG, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  for (auto &Entry : MRecvs) {
    int Rank = Entry.Key();
    send_recv &Recv = Entry.Value();
    for (long long iPoint = 0; iPoint < Recv.Count; ++iPoint) {
      long long iLinearPoint = Recv.PointIndices(iPoint);
      LinearPartition.MRanks(iLinearPoint) = Rank;
    }
  }

  for (auto &Entry : NRecvs) {
    int Rank = Entry.Key();
    send_recv &Recv = Entry.Value();
    for (long long iPoint = 0; iPoint < Recv.Count; ++iPoint) {
      long long iLinearPoint = Recv.PointIndices(iPoint);
      LinearPartition.NRanks(iLinearPoint) = Rank;
    }
  }

  for (auto &Entry : MRecvs) {
    send_recv &Recv = Entry.Value();
    for (long long iPoint = 0; iPoint < Recv.Count; ++iPoint) {
      long long iLinearPoint = Recv.PointIndices(iPoint);
      Recv.Ranks(iPoint) = LinearPartition.NRanks(iLinearPoint);
    }
  }

  for (auto &Entry : NRecvs) {
    send_recv &Recv = Entry.Value();
    for (long long iPoint = 0; iPoint < Recv.Count; ++iPoint) {
      long long iLinearPoint = Recv.PointIndices(iPoint);
      Recv.Ranks(iPoint) = LinearPartition.MRanks(iLinearPoint);
    }
  }

  for (auto &Entry : MSends) {
    int Rank = Entry.Key();
    send_recv &Send = Entry.Value();
    Irecv(Send.Ranks.Data(), Send.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Entry : NSends) {
    int Rank = Entry.Key();
    send_recv &Send = Entry.Value();
    Irecv(Send.Ranks.Data(), Send.Count, MPI_INT, Rank, 1, Comm);
  }

  for (auto &Entry : MRecvs) {
    int Rank = Entry.Key();
    send_recv &Recv = Entry.Value();
    Isend(Recv.Ranks.Data(), Recv.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Entry : NRecvs) {
    int Rank = Entry.Key();
    send_recv &Recv = Entry.Value();
    Isend(Recv.Ranks.Data(), Recv.Count, MPI_INT, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  MRecvs.Clear();
  NRecvs.Clear();

  LinearPartition = linear_partition();

  for (auto &Entry : MSends) {
    send_recv &Send = Entry.Value();
    // Reuse count for unpacking
    Send.Count = 0;
  }

  for (auto &Entry : NSends) {
    send_recv &Recv = Entry.Value();
    // Reuse count for unpacking
    Recv.Count = 0;
  }

  for (auto &IDPair : ConnectivityMGridIDs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid &MGrid = Domain.Grid(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    range_indexer NGridGlobalIndexer(NGridInfo.Cart().Range());
    local_m &LocalM = LocalMs_(IDPair);
    const connectivity_m &ConnectivityM = *LocalM.Connectivity;
    const array<int,3> &Extents = ConnectivityM.Extents();
    const array<int,2> &Destinations = ConnectivityM.Destinations();
    array<int> &DestinationRanks = LocalM.DestinationRanks;
    for (long long iDonor = 0; iDonor < ConnectivityM.Count(); ++iDonor) {
      tuple<int> CellLower = {
        Extents(0,0,iDonor),
        Extents(0,1,iDonor),
        Extents(0,2,iDonor)
      };
      if (MGrid.LocalRange().Contains(CellLower) && DestinationRanks(iDonor) < 0) {
        tuple<int> Point = {
          Destinations(0,iDonor),
          Destinations(1,iDonor),
          Destinations(2,iDonor)
        };
        long long iLinearPoint = NumPointsBeforeGrid(NGridID) + NGridGlobalIndexer.ToIndex(Point);
        int iLinearPartition = int(iLinearPoint/LinearPartitionSize);
        send_recv &Send = MSends(iLinearPartition);
        DestinationRanks(iDonor) = Send.Ranks(Send.Count);
        ++Send.Count;
      }
    }
  }

  for (auto &IDPair : ConnectivityNGridIDs) {
    int NGridID = IDPair(1);
    const grid &NGrid = Domain.Grid(NGridID);
    range_indexer NGridGlobalIndexer(NGrid.GlobalRange());
    local_n &LocalN = LocalNs_(IDPair);
    const connectivity_n &ConnectivityN = *LocalN.Connectivity;
    const array<int,2> &Points = ConnectivityN.Points();
    array<int> &SourceRanks = LocalN.SourceRanks;
    for (long long iReceiver = 0; iReceiver < ConnectivityN.Count(); ++iReceiver) {
      if (SourceRanks(iReceiver) < 0) {
        tuple<int> Point = {
          Points(0,iReceiver),
          Points(1,iReceiver),
          Points(2,iReceiver)
        };
        long long iLinearPoint = NumPointsBeforeGrid(NGridID) + NGridGlobalIndexer.ToIndex(Point);
        int iLinearPartition = int(iLinearPoint/LinearPartitionSize);
        send_recv &Send = NSends(iLinearPartition);
        SourceRanks(iReceiver) = Send.Ranks(Send.Count);
        ++Send.Count;
      }
    }
  }

  if (OVK_DEBUG) {
    for (auto &IDPair : ConnectivityMGridIDs) {
      int MGridID = IDPair(0);
      const grid &MGrid = Domain.Grid(MGridID);
      const local_m &LocalM = LocalMs_(IDPair);
      const connectivity_m &ConnectivityM = *LocalM.Connectivity;
      const array<int,3> &Extents = ConnectivityM.Extents();
      const array<int> &DestinationRanks = LocalM.DestinationRanks;
      for (long long iDonor = 0; iDonor < ConnectivityM.Count(); ++iDonor) {
        tuple<int> CellLower = {
          Extents(0,0,iDonor),
          Extents(0,1,iDonor),
          Extents(0,2,iDonor)
        };
        OVK_DEBUG_ASSERT(DestinationRanks(iDonor) >= 0 || !MGrid.LocalRange().Contains(CellLower),
          "Failed to connect donor cell (%i,%i,%i) of grid %s to receiver point.", CellLower(0),
          CellLower(1), CellLower(2), MGrid.Name());
      }
    }
    for (auto &IDPair : ConnectivityNGridIDs) {
      int NGridID = IDPair(1);
      const grid &NGrid = Domain.Grid(NGridID);
      const local_n &LocalN = LocalNs_(IDPair);
      const connectivity_n &ConnectivityN = *LocalN.Connectivity;
      const array<int,2> &Points = ConnectivityN.Points();
      const array<int> &SourceRanks = LocalN.SourceRanks;
      for (long long iReceiver = 0; iReceiver < ConnectivityN.Count(); ++iReceiver) {
        OVK_DEBUG_ASSERT(SourceRanks(iReceiver) >= 0, "Failed to connect receiver point (%i,%i,%i) "
          "of grid %s to donor cell.", Points(0,iReceiver), Points(1,iReceiver), Points(2,
          iReceiver), NGrid.Name());
      }
    }
  }

  MPI_Barrier(Comm);

}

void exchanger::ResetExchanges_() {

  if (UpdateManifest_.ResetExchanges.Empty()) return;

  const domain &Domain = *Domain_;

  for (auto &IDPair : UpdateManifest_.ResetExchanges) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid_info &MGridInfo = Domain.GridInfo(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    if (MGridInfo.IsLocal()) {
      const grid &MGrid = Domain.Grid(MGridID);
      local_m &LocalM = LocalMs_(IDPair);
      const connectivity_m &ConnectivityM = *LocalM.Connectivity;
      LocalM.Collects.Clear();
      LocalM.Sends.Clear();
      LocalM.CollectMap = core::collect_map(MGrid.Cart(), MGrid.core_Partition(),
        ConnectivityM.Extents());
      array<long long> Order = GetSendRecvOrder(ConnectivityM.Destinations(),
        NGridInfo.Cart().Range());
      LocalM.SendMap = core::send_map(ConnectivityM.Count(), std::move(Order),
        LocalM.DestinationRanks);
    }
    if (NGridInfo.IsLocal()) {
      const grid &NGrid = Domain.Grid(NGridID);
      local_n &LocalN = LocalNs_(IDPair);
      const connectivity_n &ConnectivityN = *LocalN.Connectivity;
      LocalN.Recvs.Clear();
      LocalN.Disperses.Clear();
      array<long long> Order = GetSendRecvOrder(ConnectivityN.Points(), NGrid.GlobalRange());
      LocalN.RecvMap = core::recv_map(ConnectivityN.Count(), std::move(Order),
        LocalN.SourceRanks);
      LocalN.DisperseMap = core::disperse_map(ConnectivityN.Points());
    }
  }

}

const set<int> &exchanger::CollectIDs(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const local_m &LocalM = LocalMs_(IDPair);
  const map<int,core::collect> &Collects = LocalM.Collects;

  return Collects.Keys();

}

bool exchanger::CollectExists(int MGridID, int NGridID, int CollectID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const local_m &LocalM = LocalMs_(IDPair);
  const map<int,core::collect> &Collects = LocalM.Collects;

  return Collects.Contains(CollectID);

}

void exchanger::CreateCollect(int MGridID, int NGridID, int CollectID, collect_op CollectOp,
  data_type ValueType, int Count, const range &GridValuesRange, array_layout GridValuesLayout) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");
  OVK_DEBUG_ASSERT(ValidCollectOp(CollectOp), "Invalid collect operation.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(ValidArrayLayout(GridValuesLayout), "Invalid grid values layout.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const grid &MGrid = Domain.Grid(MGridID);
  const comm &GridComm = MGrid.Comm();

  OVK_DEBUG_ASSERT(GridValuesRange.Includes(MGrid.LocalRange()), "Invalid grid values range.");

  MPI_Barrier(GridComm);

  local_m &LocalM = LocalMs_(IDPair);
  const connectivity_m &ConnectivityM = *LocalM.Connectivity;
  const core::collect_map &CollectMap = LocalM.CollectMap;
  map<int,core::collect> &Collects = LocalM.Collects;

  const cart &Cart = MGrid.Cart();
  const range &LocalRange = MGrid.LocalRange();

  OVK_DEBUG_ASSERT(!Collects.Contains(CollectID), "Collect %i already exists.", CollectID);

  core::collect Collect;

  switch (CollectOp) {
  case collect_op::NONE:
    Collect = core::CreateCollectNone(Domain.SharedContext(), GridComm, Cart, LocalRange,
      CollectMap, ValueType, Count, GridValuesRange, GridValuesLayout);
    break;
  case collect_op::ANY:
    Collect = core::CreateCollectAny(Domain.SharedContext(), GridComm, Cart, LocalRange,
      CollectMap, ValueType, Count, GridValuesRange, GridValuesLayout);
    break;
  case collect_op::NOT_ALL:
    Collect = core::CreateCollectNotAll(Domain.SharedContext(), GridComm, Cart, LocalRange,
      CollectMap, ValueType, Count, GridValuesRange, GridValuesLayout);
    break;
  case collect_op::ALL:
    Collect = core::CreateCollectAll(Domain.SharedContext(), GridComm, Cart, LocalRange,
      CollectMap, ValueType, Count, GridValuesRange, GridValuesLayout);
    break;
  case collect_op::INTERPOLATE:
    {
    floating_ref<const array<double,3>> InterpCoefs = FloatingRefRebind(ConnectivityM.
      GetFloatingRef(), ConnectivityM.InterpCoefs());
    Collect = core::CreateCollectInterp(Domain.SharedContext(), GridComm, Cart, LocalRange,
      CollectMap, ValueType, Count, GridValuesRange, GridValuesLayout, InterpCoefs);
    }
    break;
  }

  Collects.Insert(CollectID, std::move(Collect));

  MPI_Barrier(GridComm);

}

void exchanger::DestroyCollect(int MGridID, int NGridID, int CollectID) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);

  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const grid &MGrid = Domain.Grid(MGridID);
  const comm &GridComm = MGrid.Comm();

  MPI_Barrier(GridComm);

  local_m &LocalM = LocalMs_(IDPair);
  map<int,core::collect> &Collects = LocalM.Collects;

  OVK_DEBUG_ASSERT(Collects.Contains(CollectID), "Collect %i does not exist.", CollectID);

  Collects.Erase(CollectID);

  MPI_Barrier(GridComm);

}

void exchanger::Collect(int MGridID, int NGridID, int CollectID, const void *GridValues, void
  *DonorValues) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const grid &MGrid = Domain.Grid(MGridID);
  const comm &GridComm = MGrid.Comm();

  MPI_Barrier(GridComm);

  local_m &LocalM = LocalMs_(IDPair);
  map<int,core::collect> &Collects = LocalM.Collects;

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.Start(COLLECT_TIME);

  OVK_DEBUG_ASSERT(Collects.Contains(CollectID), "Collect %i does not exist.", CollectID);

  core::collect &Collect = Collects(CollectID);

  Collect.Collect(GridValues, DonorValues);

  Profiler.Stop(COLLECT_TIME);

  MPI_Barrier(GridComm);

}

const set<int> &exchanger::SendIDs(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const local_m &LocalM = LocalMs_(IDPair);
  const map<int,core::send> &Sends = LocalM.Sends;

  return Sends.Keys();

}

bool exchanger::SendExists(int MGridID, int NGridID, int SendID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const local_m &LocalM = LocalMs_(IDPair);
  const map<int,core::send> &Sends = LocalM.Sends;

  return Sends.Contains(SendID);

}

void exchanger::CreateSend(int MGridID, int NGridID, int SendID, data_type ValueType, int Count, int
  Tag) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(Tag >= 0, "Invalid tag.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  local_m &LocalM = LocalMs_(IDPair);
  const core::send_map &SendMap = LocalM.SendMap;
  map<int,core::send> &Sends = LocalM.Sends;

  OVK_DEBUG_ASSERT(!Sends.Contains(SendID), "Send %i already exists.", SendID);

  const elem_set<int,2> &ConnectivityIDs = ConnectivityComponent.ConnectivityIDs();
  int GlobalTagMultiplier = ConnectivityIDs.Count();
  int GlobalTagOffset = ConnectivityIDs.Find(IDPair) - ConnectivityIDs.Begin();
  int GlobalTag = GlobalTagMultiplier*Tag + GlobalTagOffset;

  core::send Send = core::CreateSend(Context_, Domain.Comm(), SendMap, ValueType, Count, GlobalTag);

  Sends.Insert(SendID, std::move(Send));

}

void exchanger::DestroySend(int MGridID, int NGridID, int SendID) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  local_m &LocalM = LocalMs_(IDPair);
  map<int,core::send> &Sends = LocalM.Sends;

  OVK_DEBUG_ASSERT(Sends.Contains(SendID), "Send %i does not exist.", SendID);

  Sends.Erase(SendID);

}

request exchanger::Send(int MGridID, int NGridID, int SendID, const void *DonorValues) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  local_m &LocalM = LocalMs_(IDPair);
  map<int,core::send> &Sends = LocalM.Sends;

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.Start(SEND_RECV_TIME);

  OVK_DEBUG_ASSERT(Sends.Contains(SendID), "Send %i does not exist.", SendID);

  core::send &Send = Sends(SendID);

  request Request = Send.Send(DonorValues);

  Profiler.Stop(SEND_RECV_TIME);

  return Request;

}

const set<int> &exchanger::ReceiveIDs(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const local_n &LocalN = LocalNs_(IDPair);
  const map<int,core::recv> &Recvs = LocalN.Recvs;

  return Recvs.Keys();

}

bool exchanger::ReceiveExists(int MGridID, int NGridID, int RecvID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const local_n &LocalN = LocalNs_(IDPair);
  const map<int,core::recv> &Recvs = LocalN.Recvs;

  return Recvs.Contains(RecvID);

}

void exchanger::CreateReceive(int MGridID, int NGridID, int RecvID, data_type ValueType, int Count,
  int Tag) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(Tag >= 0, "Invalid tag.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_(IDPair);
  const core::recv_map &RecvMap = LocalN.RecvMap;
  map<int,core::recv> &Recvs = LocalN.Recvs;

  OVK_DEBUG_ASSERT(!Recvs.Contains(RecvID), "Receive %i already exists.", RecvID);

  const elem_set<int,2> &ConnectivityIDs = ConnectivityComponent.ConnectivityIDs();
  int GlobalTagMultiplier = ConnectivityIDs.Count();
  int GlobalTagOffset = ConnectivityIDs.Find(IDPair) - ConnectivityIDs.Begin();
  int GlobalTag = GlobalTagMultiplier*Tag + GlobalTagOffset;

  core::recv Recv = core::CreateRecv(Context_, Domain.Comm(), RecvMap, ValueType, Count, GlobalTag);

  Recvs.Insert(RecvID, std::move(Recv));

}

void exchanger::DestroyReceive(int MGridID, int NGridID, int RecvID) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_({MGridID,NGridID});
  map<int,core::recv> &Recvs = LocalN.Recvs;

  OVK_DEBUG_ASSERT(Recvs.Contains(RecvID), "Receive %i does not exist.", RecvID);

  Recvs.Erase(RecvID);

}

request exchanger::Receive(int MGridID, int NGridID, int RecvID, void *ReceiverValues) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_(IDPair);
  map<int,core::recv> &Recvs = LocalN.Recvs;

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.Start(SEND_RECV_TIME);

  OVK_DEBUG_ASSERT(Recvs.Contains(RecvID), "Receive %i does not exist.", RecvID);

  core::recv &Recv = Recvs(RecvID);

  request Request = Recv.Recv(ReceiverValues);

  Profiler.Stop(SEND_RECV_TIME);

  return Request;

}

const set<int> &exchanger::DisperseIDs(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const local_n &LocalN = LocalNs_(IDPair);
  const map<int,core::disperse> &Disperses = LocalN.Disperses;

  return Disperses.Keys();

}

bool exchanger::DisperseExists(int MGridID, int NGridID, int DisperseID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const local_n &LocalN = LocalNs_(IDPair);
  const map<int,core::disperse> &Disperses = LocalN.Disperses;

  return Disperses.Contains(DisperseID);

}

void exchanger::CreateDisperse(int MGridID, int NGridID, int DisperseID, disperse_op DisperseOp,
  data_type ValueType, int Count, const range &GridValuesRange, array_layout GridValuesLayout) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");
  OVK_DEBUG_ASSERT(ValidDisperseOp(DisperseOp), "Invalid disperse operation.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(ValidArrayLayout(GridValuesLayout), "Invalid grid values layout.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const grid &NGrid = Domain.Grid(NGridID);

  OVK_DEBUG_ASSERT(GridValuesRange.Includes(NGrid.LocalRange()), "Invalid grid values range.");

  local_n &LocalN = LocalNs_(IDPair);
  const core::disperse_map &DisperseMap = LocalN.DisperseMap;
  map<int,core::disperse> &Disperses = LocalN.Disperses;

  OVK_DEBUG_ASSERT(!Disperses.Contains(DisperseID), "Disperse %i already exists.", DisperseID);

  core::disperse Disperse;

  switch (DisperseOp) {
  case disperse_op::OVERWRITE:
    Disperse = core::CreateDisperseOverwrite(Context_, DisperseMap, ValueType, Count,
      GridValuesRange, GridValuesLayout);
    break;
  case disperse_op::APPEND:
    Disperse = core::CreateDisperseAppend(Context_, DisperseMap, ValueType, Count,
      GridValuesRange, GridValuesLayout);
    break;
  }

  Disperses.Insert(DisperseID, std::move(Disperse));

}

void exchanger::DestroyDisperse(int MGridID, int NGridID, int DisperseID) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_(IDPair);
  map<int,core::disperse> &Disperses = LocalN.Disperses;

  OVK_DEBUG_ASSERT(Disperses.Contains(DisperseID), "Disperse %i does not exist.", DisperseID);

  Disperses.Erase(DisperseID);

}

void exchanger::Disperse(int MGridID, int NGridID, int DisperseID, const void *ReceiverValues, void
  *GridValues) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");

  elem<int,2> IDPair = {MGridID,NGridID};

  auto &ConnectivityComponent = Domain.Component<connectivity_component>(ConnectivityComponentID_);

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(IDPair), "Connectivity (%i,%i) does "
    "not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_(IDPair);
  map<int,core::disperse> &Disperses = LocalN.Disperses;

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.Start(DISPERSE_TIME);

  OVK_DEBUG_ASSERT(Disperses.Contains(DisperseID), "Disperse %i does not exist.", DisperseID);

  core::disperse &Disperse = Disperses(DisperseID);

  Disperse.Disperse(ReceiverValues, GridValues);

  Profiler.Stop(DISPERSE_TIME);

}

exchanger::params &exchanger::params::SetName(std::string Name) {

  Name_ = std::move(Name);

  return *this;

}

exchanger::bindings &exchanger::bindings::SetConnectivityComponentID(int ConnectivityComponentID) {

  OVK_DEBUG_ASSERT(ConnectivityComponentID >= 0, "Invalid connectivity component ID.");

  ConnectivityComponentID_ = ConnectivityComponentID;

  return *this;

}

namespace {

long long BinDivide(long long N, int NumBins) {

  long long NumBinsLongLong = (long long)(NumBins);

  return (N + NumBinsLongLong - 1)/NumBinsLongLong;

}

array<long long> GetSendRecvOrder(const array<int,2> &ReceiverPoints, const range
  &ReceiverGridGlobalRange) {

  long long NumReceivers = ReceiverPoints.Size(1);

  array<long long> Order({NumReceivers});

  using range_indexer = indexer<long long, int, MAX_DIMS, array_layout::COLUMN_MAJOR>;
  range_indexer ReceiverGridGlobalIndexer(ReceiverGridGlobalRange);

  array<long long> ReceiverIndices({NumReceivers});

  for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
    tuple<int> Point = {
      ReceiverPoints(0,iReceiver),
      ReceiverPoints(1,iReceiver),
      ReceiverPoints(2,iReceiver)
    };
    ReceiverIndices(iReceiver) = ReceiverGridGlobalIndexer.ToIndex(Point);
  }

  bool Sorted = true;

  long long PrevIndex = 0;
  for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
    if (ReceiverIndices(iReceiver) < PrevIndex) {
      Sorted = false;
      break;
    }
    PrevIndex = ReceiverIndices(iReceiver);
  }

  if (Sorted) {
    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      Order(iReceiver) = iReceiver;
    }
  } else {
    core::SortPermutation(ReceiverIndices, Order);
  }

  return Order;

}

}

}
