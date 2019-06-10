// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
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
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/IDMap.hpp"
#include "ovk/core/IDSet.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Recv.hpp"
#include "ovk/core/RecvMap.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/Send.hpp"
#include "ovk/core/SendMap.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <map>
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
  FloatingRefGenerator_(*this),
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

  int ConnectivityComponentID = Bindings.ConnectivityComponentID_;

  OVK_DEBUG_ASSERT(ConnectivityComponentID >= 0, "Invalid connectivity component ID.");
  OVK_DEBUG_ASSERT(Domain.ComponentExists(ConnectivityComponentID), "Component %i does not "
    "exist.", ConnectivityComponentID);

  floating_ref<exchanger> FloatingRef = FloatingRefGenerator_.Generate();

  const connectivity_component &ConnectivityComponent = Domain.Component<connectivity_component>(
    ConnectivityComponentID);
  ConnectivityComponent_ = ConnectivityComponent.GetFloatingRef();

  ComponentEventListener_ = Domain.AddComponentEventListener([FloatingRef, ConnectivityComponentID](
    int ComponentID, component_event_flags Flags) {
    exchanger &Exchanger = *FloatingRef;
    if (ComponentID == ConnectivityComponentID && (Flags & component_event_flags::DESTROY) !=
      component_event_flags::NONE) {
      // Slightly questionable, since Unbind deletes this lambda, but I think we're ok as long as
      // no captured data is used afterwards
      Exchanger.Unbind();
    }
  });

  ConnectivityEventListener_ = ConnectivityComponent.AddConnectivityEventListener([FloatingRef](int
    MGridID, int NGridID, connectivity_event_flags Flags, bool LastInSequence) {
    exchanger &Exchanger = *FloatingRef;
    auto &AccumulatedFlags = Exchanger.ConnectivityEventFlags_.Get({MGridID,NGridID},
      connectivity_event_flags::NONE);
    AccumulatedFlags |= Flags; 
    if (LastInSequence) {
      Exchanger.OnConnectivityEvent_();
    }
  });

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.CommRank() == 0, 0, "Bound exchanger %s to domain %s.", *Name_,
    Domain.Name());

  for (auto &IDPair : ConnectivityComponent.ConnectivityIDs()) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    ConnectivityEventFlags_.Insert({MGridID,NGridID}, connectivity_event_flags::CREATE |
      connectivity_event_flags::ALL_EDITS);
  }

  OnConnectivityEvent_();

}

void exchanger::Unbind() {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  MPI_Barrier(Domain.Comm());

  LocalMs_.Clear();
  LocalNs_.Clear();

  ConnectivityEventListener_.Reset();

  ComponentEventListener_.Reset();

  ConnectivityComponent_.Reset();
  Domain_.Reset();

  MPI_Barrier(Domain.Comm());

  core::logger &Logger = Context_->core_Logger();
  Logger.LogStatus(Domain.CommRank() == 0, 0, "Unbound exchanger %s.", *Name_);

}

const domain &exchanger::Domain() const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  return *Domain_;

}

void exchanger::OnConnectivityEvent_() {

  const domain &Domain = *Domain_;
  core::logger &Logger = Context_->core_Logger();

  Logger.LogStatus(Domain.CommRank() == 0, 0, "Updating exchanger %s...", *Name_);

  DestroyExchangesForDyingConnectivities_();
  CreateExchangesForNewConnectivities_();
  UpdateExchangesForModifiedConnectivities_();

  ConnectivityEventFlags_.Clear();

  Logger.LogStatus(Domain.CommRank() == 0, 0, "Done updating exchanger %s.", *Name_);

}

void exchanger::DestroyExchangesForDyingConnectivities_() {

  for (auto &EventEntry : ConnectivityEventFlags_) {
    int MGridID = EventEntry.Key(0);
    int NGridID = EventEntry.Key(1);
    connectivity_event_flags Flags = EventEntry.Value();
    if ((Flags & connectivity_event_flags::DESTROY) != connectivity_event_flags::NONE) {
      LocalMs_.Erase(MGridID,NGridID);
      LocalNs_.Erase(MGridID,NGridID);
    }
  }

}

void exchanger::CreateExchangesForNewConnectivities_() {

  const domain &Domain = *Domain_;
  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  for (auto &EventEntry : ConnectivityEventFlags_) {
    int MGridID = EventEntry.Key(0);
    int NGridID = EventEntry.Key(1);
    connectivity_event_flags Flags = EventEntry.Value();
    if ((Flags & connectivity_event_flags::CREATE) != connectivity_event_flags::NONE) {
      if (Domain.GridIsLocal(MGridID)) {
        local_m &LocalM = LocalMs_.Insert(MGridID,NGridID);
        LocalM.Connectivity = &ConnectivityComponent.ConnectivityM(MGridID, NGridID);
      }
      if (Domain.GridIsLocal(NGridID)) {
        local_n &LocalN = LocalNs_.Insert(MGridID,NGridID);
        LocalN.Connectivity = &ConnectivityComponent.ConnectivityN(MGridID, NGridID);
      }
    }
  }

}

void exchanger::UpdateExchangesForModifiedConnectivities_() {

  UpdateSourceDestRanks_();
  PurgeExchanges_();

}

void exchanger::UpdateSourceDestRanks_() {

  using range_indexer = indexer<long long, int, 3, ovk::array_layout::GRID>;

  const domain &Domain = *Domain_;
  const core::comm &Comm = Domain.core_Comm();

  MPI_Barrier(Comm);
  
  for (auto &EventEntry : ConnectivityEventFlags_) {
    int MGridID = EventEntry.Key(0);
    int NGridID = EventEntry.Key(1);
    connectivity_event_flags Flags = EventEntry.Value();
    connectivity_event_flags MEditMask = connectivity_event_flags::RESIZE_M |
      connectivity_event_flags::EDIT_M_DESTINATIONS;
    connectivity_event_flags NEditMask = connectivity_event_flags::RESIZE_N |
      connectivity_event_flags::EDIT_N_SOURCES;
    if ((Flags & MEditMask) != connectivity_event_flags::NONE) {
      if (Domain.GridIsLocal(MGridID)) {
        const grid &MGrid = Domain.Grid(MGridID);
        local_m &LocalM = LocalMs_(MGridID,NGridID);
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
    }
    if ((Flags & NEditMask) != connectivity_event_flags::NONE) {
      if (Domain.GridIsLocal(NGridID)) {
        local_n &LocalN = LocalNs_(MGridID,NGridID);
        const connectivity_n &ConnectivityN = *LocalN.Connectivity;
        LocalN.SourceRanks = ConnectivityN.SourceRanks();
      }
    }
  }

  id_set<2> ConnectivityMGridIDs;
  id_set<2> ConnectivityNGridIDs;

  for (auto &EventEntry : ConnectivityEventFlags_) {
    int MGridID = EventEntry.Key(0);
    int NGridID = EventEntry.Key(1);
    connectivity_event_flags Flags = EventEntry.Value();
    connectivity_event_flags Mask = connectivity_event_flags::EDIT_M_DESTINATIONS |
      connectivity_event_flags::EDIT_N_SOURCES;
    if ((Flags & Mask) != connectivity_event_flags::NONE) {
      if (Domain.GridIsLocal(MGridID)) {
        ConnectivityMGridIDs.Insert(MGridID,NGridID);
      }
      if (Domain.GridIsLocal(NGridID)) {
        ConnectivityNGridIDs.Insert(MGridID,NGridID);
      }
    }
  }

  long long TotalPoints = 0;
  for (int GridID : Domain.GridIDs()) {
    const grid_info &Info = Domain.GridInfo(GridID);
    TotalPoints += Info.Cart().Range().Count();
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

  id_map<1,long long> NumPointsBeforeGrid;
  long long NumPointsPartial = 0;
  for (int GridID : Domain.GridIDs()) {
    const grid_info &Info = Domain.GridInfo(GridID);
    NumPointsBeforeGrid.Insert(GridID, NumPointsPartial);
    NumPointsPartial += Info.Cart().Range().Count();
  }

  struct send_recv {
    long long Count;
    array<long long> PointIndices;
    array<int> Ranks;
    send_recv():
      Count(0)
    {}
  };

  std::map<int, send_recv> MSends, NSends;

  for (auto &IDPair : ConnectivityMGridIDs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid &MGrid = Domain.Grid(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    range_indexer NGridGlobalIndexer(NGridInfo.Cart().Range());
    const local_m &LocalM = LocalMs_(MGridID,NGridID);
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
        send_recv &Send = MSends[iLinearPartition];
        ++Send.Count;
      }
    }
  }

  for (auto &IDPair : ConnectivityNGridIDs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid &NGrid = Domain.Grid(NGridID);
    range_indexer NGridGlobalIndexer(NGrid.GlobalRange());
    const local_n &LocalN = LocalNs_(MGridID,NGridID);
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
        send_recv &Send = NSends[iLinearPartition];
        ++Send.Count;
      }
    }
  }

  for (auto &Pair : MSends) {
    send_recv &Send = Pair.second;
    Send.PointIndices.Reserve(Send.Count);
  }

  for (auto &Pair : NSends) {
    send_recv &Send = Pair.second;
    Send.PointIndices.Reserve(Send.Count);
  }

  for (auto &IDPair : ConnectivityMGridIDs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid &MGrid = Domain.Grid(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    range_indexer NGridGlobalIndexer(NGridInfo.Cart().Range());
    const local_m &LocalM = LocalMs_(MGridID,NGridID);
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
        send_recv &Send = MSends.at(iLinearPartition);
        Send.PointIndices.Append(iLinearPoint);
      }
    }
  }

  for (auto &IDPair : ConnectivityNGridIDs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid &NGrid = Domain.Grid(NGridID);
    range_indexer NGridGlobalIndexer(NGrid.GlobalRange());
    const local_n &LocalN = LocalNs_(MGridID,NGridID);
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
        send_recv &Send = NSends.at(iLinearPartition);
        Send.PointIndices.Append(iLinearPoint);
      }
    }
  }

  for (auto &Pair : MSends) {
    send_recv &Send = Pair.second;
    Send.Ranks.Resize({Send.Count});
  }

  for (auto &Pair : NSends) {
    send_recv &Send = Pair.second;
    Send.Ranks.Resize({Send.Count});
  }

  int NumMSends = MSends.size();
  int NumNSends = NSends.size();

  array<int> MSendToRanks, NSendToRanks;

  MSendToRanks.Reserve(NumMSends);
  for (auto &Pair : MSends) {
    int Rank = Pair.first;
    MSendToRanks.Append(Rank);
  }

  NSendToRanks.Reserve(NumNSends);
  for (auto &Pair : NSends) {
    int Rank = Pair.first;
    NSendToRanks.Append(Rank);
  }

  array<int> MRecvFromRanks = core::DynamicHandshake(Comm, MSendToRanks);
  array<int> NRecvFromRanks = core::DynamicHandshake(Comm, NSendToRanks);

  MSendToRanks.Clear();
  NSendToRanks.Clear();

  std::map<int, send_recv> MRecvs, NRecvs;

  for (int Rank : MRecvFromRanks) {
    MRecvs.emplace(Rank, send_recv());
  }

  for (int Rank : NRecvFromRanks) {
    NRecvs.emplace(Rank, send_recv());
  }

  MRecvFromRanks.Clear();
  NRecvFromRanks.Clear();

  int NumMRecvs = MRecvs.size();
  int NumNRecvs = NRecvs.size();

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

  for (auto &Pair : MRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    Irecv(&Recv.Count, 1, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Pair : NRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    Irecv(&Recv.Count, 1, MPI_LONG_LONG, Rank, 1, Comm);
  }

  for (auto &Pair : MSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Isend(&Send.Count, 1, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Pair : NSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Isend(&Send.Count, 1, MPI_LONG_LONG, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  for (auto &Pair : MRecvs) {
    send_recv &Recv = Pair.second;
    Recv.PointIndices.Resize({Recv.Count});
    Recv.Ranks.Resize({Recv.Count});
  }

  for (auto &Pair : NRecvs) {
    send_recv &Recv = Pair.second;
    Recv.PointIndices.Resize({Recv.Count});
    Recv.Ranks.Resize({Recv.Count});
  }

  for (auto &Pair : MRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    Irecv(Recv.PointIndices.Data(), Recv.Count, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Pair : NRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    Irecv(Recv.PointIndices.Data(), Recv.Count, MPI_LONG_LONG, Rank, 1, Comm);
  }

  for (auto &Pair : MSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Isend(Send.PointIndices.Data(), Send.Count, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Pair : NSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Isend(Send.PointIndices.Data(), Send.Count, MPI_LONG_LONG, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  for (auto &Pair : MRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    for (long long iPoint = 0; iPoint < Recv.Count; ++iPoint) {
      long long iLinearPoint = Recv.PointIndices(iPoint);
      LinearPartition.MRanks(iLinearPoint) = Rank;
    }
  }

  for (auto &Pair : NRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    for (long long iPoint = 0; iPoint < Recv.Count; ++iPoint) {
      long long iLinearPoint = Recv.PointIndices(iPoint);
      LinearPartition.NRanks(iLinearPoint) = Rank;
    }
  }

  for (auto &Pair : MRecvs) {
    send_recv &Recv = Pair.second;
    for (long long iPoint = 0; iPoint < Recv.Count; ++iPoint) {
      long long iLinearPoint = Recv.PointIndices(iPoint);
      Recv.Ranks(iPoint) = LinearPartition.NRanks(iLinearPoint);
    }
  }

  for (auto &Pair : NRecvs) {
    send_recv &Recv = Pair.second;
    for (long long iPoint = 0; iPoint < Recv.Count; ++iPoint) {
      long long iLinearPoint = Recv.PointIndices(iPoint);
      Recv.Ranks(iPoint) = LinearPartition.MRanks(iLinearPoint);
    }
  }

  for (auto &Pair : MSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Irecv(Send.Ranks.Data(), Send.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : NSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Irecv(Send.Ranks.Data(), Send.Count, MPI_INT, Rank, 1, Comm);
  }

  for (auto &Pair : MRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    Isend(Recv.Ranks.Data(), Recv.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : NRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    Isend(Recv.Ranks.Data(), Recv.Count, MPI_INT, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  MRecvs.clear();
  NRecvs.clear();

  LinearPartition = linear_partition();

  for (auto &Pair : MSends) {
    send_recv &Send = Pair.second;
    // Reuse count for unpacking
    Send.Count = 0;
  }

  for (auto &Pair : NSends) {
    send_recv &Recv = Pair.second;
    // Reuse count for unpacking
    Recv.Count = 0;
  }

  for (auto &IDPair : ConnectivityMGridIDs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid &MGrid = Domain.Grid(MGridID);
    const grid_info &NGridInfo = Domain.GridInfo(NGridID);
    range_indexer NGridGlobalIndexer(NGridInfo.Cart().Range());
    local_m &LocalM = LocalMs_(MGridID,NGridID);
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
        send_recv &Send = MSends.at(iLinearPartition);
        DestinationRanks(iDonor) = Send.Ranks(Send.Count);
        ++Send.Count;
      }
    }
  }

  for (auto &IDPair : ConnectivityNGridIDs) {
    int MGridID = IDPair(0);
    int NGridID = IDPair(1);
    const grid &NGrid = Domain.Grid(NGridID);
    range_indexer NGridGlobalIndexer(NGrid.GlobalRange());
    local_n &LocalN = LocalNs_(MGridID,NGridID);
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
        send_recv &Send = NSends.at(iLinearPartition);
        SourceRanks(iReceiver) = Send.Ranks(Send.Count);
        ++Send.Count;
      }
    }
  }

  if (OVK_DEBUG) {
    for (auto &IDPair : ConnectivityMGridIDs) {
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      const grid &MGrid = Domain.Grid(MGridID);
      const local_m &LocalM = LocalMs_(MGridID,NGridID);
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
      int MGridID = IDPair(0);
      int NGridID = IDPair(1);
      const grid &NGrid = Domain.Grid(NGridID);
      const local_n &LocalN = LocalNs_(MGridID,NGridID);
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

void exchanger::PurgeExchanges_() {

  const domain &Domain = *Domain_;

  for (auto &EventEntry : ConnectivityEventFlags_) {
    int MGridID = EventEntry.Key(0);
    int NGridID = EventEntry.Key(1);
    connectivity_event_flags Flags = EventEntry.Value();
    if ((Flags & connectivity_event_flags::ALL_EDITS) != connectivity_event_flags::NONE) {
      const grid_info &MGridInfo = Domain.GridInfo(MGridID);
      const grid_info &NGridInfo = Domain.GridInfo(NGridID);
      if (MGridInfo.IsLocal()) {
        const grid &MGrid = Domain.Grid(MGridID);
        local_m &LocalM = LocalMs_(MGridID,NGridID);
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
        local_n &LocalN = LocalNs_(MGridID,NGridID);
        const connectivity_n &ConnectivityN = *LocalN.Connectivity;
        LocalN.Recvs.Clear();
        LocalN.Disperses.Clear();
        array<long long> Order = GetSendRecvOrder(ConnectivityN.Points(), NGrid.GlobalRange());
        LocalN.RecvMap = core::recv_map(ConnectivityN.Count(), std::move(Order),
          LocalN.SourceRanks);
      }
    }
  }

}

const id_set<1> &exchanger::CollectIDs(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const local_m &LocalM = LocalMs_(MGridID,NGridID);
  const id_map<1,core::collect> &Collects = LocalM.Collects;

  return Collects.Keys();

}

bool exchanger::CollectExists(int MGridID, int NGridID, int CollectID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const local_m &LocalM = LocalMs_(MGridID,NGridID);
  const id_map<1,core::collect> &Collects = LocalM.Collects;

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

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const grid &MGrid = Domain.Grid(MGridID);
  const core::comm &GridComm = MGrid.core_Comm();

  OVK_DEBUG_ASSERT(GridValuesRange.Includes(MGrid.LocalRange()), "Invalid grid values range.");

  MPI_Barrier(GridComm);

  local_m &LocalM = LocalMs_(MGridID,NGridID);
  const connectivity_m &ConnectivityM = *LocalM.Connectivity;
  const core::collect_map &CollectMap = LocalM.CollectMap;
  id_map<1,core::collect> &Collects = LocalM.Collects;

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
    Collect = core::CreateCollectInterp(Domain.SharedContext(), GridComm, Cart, LocalRange,
      CollectMap, ValueType, Count, GridValuesRange, GridValuesLayout, ConnectivityM.InterpCoefs());
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

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);

  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const grid &MGrid = Domain.Grid(MGridID);
  const core::comm &GridComm = MGrid.core_Comm();

  MPI_Barrier(GridComm);

  local_m &LocalM = LocalMs_(MGridID,NGridID);
  id_map<1,core::collect> &Collects = LocalM.Collects;

  OVK_DEBUG_ASSERT(Collects.Contains(CollectID), "Collect %i does not exist.", CollectID);

  Collects.Erase(CollectID);

  MPI_Barrier(GridComm);

}

void exchanger::Collect(int MGridID, int NGridID, int CollectID, const void * const *GridValues,
  void **DonorValues) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const grid &MGrid = Domain.Grid(MGridID);
  const core::comm &GridComm = MGrid.core_Comm();

  MPI_Barrier(GridComm);

  local_m &LocalM = LocalMs_(MGridID,NGridID);
  id_map<1,core::collect> &Collects = LocalM.Collects;

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.Start(COLLECT_TIME);

  OVK_DEBUG_ASSERT(Collects.Contains(CollectID), "Collect %i does not exist.", CollectID);

  core::collect &Collect = Collects(CollectID);

  Collect.Collect(GridValues, DonorValues);

  Profiler.Stop(COLLECT_TIME);

  MPI_Barrier(GridComm);

}

const id_set<1> &exchanger::SendIDs(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const local_m &LocalM = LocalMs_(MGridID,NGridID);
  const id_map<1,core::send> &Sends = LocalM.Sends;

  return Sends.Keys();

}

bool exchanger::SendExists(int MGridID, int NGridID, int SendID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  const local_m &LocalM = LocalMs_(MGridID,NGridID);
  const id_map<1,core::send> &Sends = LocalM.Sends;

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

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  local_m &LocalM = LocalMs_(MGridID,NGridID);
  const core::send_map &SendMap = LocalM.SendMap;
  id_map<1,core::send> &Sends = LocalM.Sends;

  OVK_DEBUG_ASSERT(!Sends.Contains(SendID), "Send %i already exists.", SendID);

  const id_set<2> &ConnectivityIDs = ConnectivityComponent.ConnectivityIDs();
  int GlobalTagMultiplier = ConnectivityIDs.Count();
  int GlobalTagOffset = ConnectivityIDs.Find(MGridID,NGridID) - ConnectivityIDs.Begin();
  int GlobalTag = GlobalTagMultiplier*Tag + GlobalTagOffset;

  core::send Send = core::CreateSend(Context_, Domain.core_Comm(), SendMap, ValueType, Count,
    GlobalTag);

  Sends.Insert(SendID, std::move(Send));

}

void exchanger::DestroySend(int MGridID, int NGridID, int SendID) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  local_m &LocalM = LocalMs_(MGridID,NGridID);
  id_map<1,core::send> &Sends = LocalM.Sends;

  OVK_DEBUG_ASSERT(Sends.Contains(SendID), "Send %i does not exist.", SendID);

  Sends.Erase(SendID);

}

request exchanger::Send(int MGridID, int NGridID, int SendID, const void * const *DonorValues) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(MGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(MGridID).Name());

  local_m &LocalM = LocalMs_(MGridID,NGridID);
  id_map<1,core::send> &Sends = LocalM.Sends;

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.Start(SEND_RECV_TIME);

  OVK_DEBUG_ASSERT(Sends.Contains(SendID), "Send %i does not exist.", SendID);

  core::send &Send = Sends(SendID);

  request Request = Send.Send(DonorValues);

  Profiler.Stop(SEND_RECV_TIME);

  return Request;

}

const id_set<1> &exchanger::ReceiveIDs(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const local_n &LocalN = LocalNs_(MGridID,NGridID);
  const id_map<1,core::recv> &Recvs = LocalN.Recvs;

  return Recvs.Keys();

}

bool exchanger::ReceiveExists(int MGridID, int NGridID, int RecvID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const local_n &LocalN = LocalNs_(MGridID,NGridID);
  const id_map<1,core::recv> &Recvs = LocalN.Recvs;

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

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_(MGridID,NGridID);
  const core::recv_map &RecvMap = LocalN.RecvMap;
  id_map<1,core::recv> &Recvs = LocalN.Recvs;

  OVK_DEBUG_ASSERT(!Recvs.Contains(RecvID), "Receive %i already exists.", RecvID);

  const id_set<2> &ConnectivityIDs = ConnectivityComponent.ConnectivityIDs();
  int GlobalTagMultiplier = ConnectivityIDs.Count();
  int GlobalTagOffset = ConnectivityIDs.Find(MGridID,NGridID) - ConnectivityIDs.Begin();
  int GlobalTag = GlobalTagMultiplier*Tag + GlobalTagOffset;

  core::recv Recv = core::CreateRecv(Context_, Domain.core_Comm(), RecvMap, ValueType, Count,
    GlobalTag);

  Recvs.Insert(RecvID, std::move(Recv));

}

void exchanger::DestroyReceive(int MGridID, int NGridID, int RecvID) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_(MGridID,NGridID);
  id_map<1,core::recv> &Recvs = LocalN.Recvs;

  OVK_DEBUG_ASSERT(Recvs.Contains(RecvID), "Receive %i does not exist.", RecvID);

  Recvs.Erase(RecvID);

}

request exchanger::Receive(int MGridID, int NGridID, int RecvID, void **ReceiverValues) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_(MGridID,NGridID);
  id_map<1,core::recv> &Recvs = LocalN.Recvs;

  core::profiler &Profiler = Context_->core_Profiler();

  Profiler.Start(SEND_RECV_TIME);

  OVK_DEBUG_ASSERT(Recvs.Contains(RecvID), "Receive %i does not exist.", RecvID);

  core::recv &Recv = Recvs(RecvID);

  request Request = Recv.Recv(ReceiverValues);

  Profiler.Stop(SEND_RECV_TIME);

  return Request;

}

const id_set<1> &exchanger::DisperseIDs(int MGridID, int NGridID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const local_n &LocalN = LocalNs_(MGridID,NGridID);
  const id_map<1,core::disperse> &Disperses = LocalN.Disperses;

  return Disperses.Keys();

}

bool exchanger::DisperseExists(int MGridID, int NGridID, int DisperseID) const {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const local_n &LocalN = LocalNs_(MGridID,NGridID);
  const id_map<1,core::disperse> &Disperses = LocalN.Disperses;

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

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  const grid &NGrid = Domain.Grid(NGridID);

  OVK_DEBUG_ASSERT(GridValuesRange.Includes(NGrid.LocalRange()), "Invalid grid values range.");

  local_n &LocalN = LocalNs_(MGridID,NGridID);
  const connectivity_n &ConnectivityN = *LocalN.Connectivity;
  id_map<1,core::disperse> &Disperses = LocalN.Disperses;

  OVK_DEBUG_ASSERT(!Disperses.Contains(DisperseID), "Disperse %i already exists.", DisperseID);

  core::disperse Disperse;

  switch (DisperseOp) {
  case disperse_op::OVERWRITE:
    Disperse = core::CreateDisperseOverwrite(Context_, ConnectivityN.Points(), ValueType, Count,
      GridValuesRange, GridValuesLayout);
    break;
  case disperse_op::APPEND:
    Disperse = core::CreateDisperseAppend(Context_, ConnectivityN.Points(), ValueType, Count,
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

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_(MGridID,NGridID);
  id_map<1,core::disperse> &Disperses = LocalN.Disperses;

  OVK_DEBUG_ASSERT(Disperses.Contains(DisperseID), "Disperse %i does not exist.", DisperseID);

  Disperses.Erase(DisperseID);

}

void exchanger::Disperse(int MGridID, int NGridID, int DisperseID, const void * const
  *ReceiverValues, void **GridValues) {

  OVK_DEBUG_ASSERT(Domain_, "Exchanger is not bound to a domain.");

  const domain &Domain = *Domain_;

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");

  const connectivity_component &ConnectivityComponent = *ConnectivityComponent_;

  OVK_DEBUG_ASSERT(ConnectivityComponent.ConnectivityExists(MGridID, NGridID), "Connectivity "
    "(%i,%i) does not exist.", MGridID, NGridID);
  OVK_DEBUG_ASSERT(Domain.GridIsLocal(NGridID), "Grid %s is not local to rank @rank@.",
    Domain.GridInfo(NGridID).Name());

  local_n &LocalN = LocalNs_(MGridID,NGridID);
  id_map<1,core::disperse> &Disperses = LocalN.Disperses;

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

namespace {

long long BinDivide(long long N, int NumBins) {

  long long NumBinsLongLong = (long long)(NumBins);

  return (N + NumBinsLongLong - 1)/NumBinsLongLong;

}

array<long long> GetSendRecvOrder(const array<int,2> &ReceiverPoints, const range
  &ReceiverGridGlobalRange) {

  long long NumReceivers = ReceiverPoints.Size(1);

  array<long long> Order({NumReceivers});

  using range_indexer = indexer<long long, int, MAX_DIMS, array_layout::GRID>;
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
