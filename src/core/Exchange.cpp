// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Exchange.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Collect.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Disperse.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Recv.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/Send.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

namespace ovk {

namespace {

void UpdateCollectSendInfo(exchange &Exchange);
void UpdateCollectReceiveInfo(exchange &Exchange);

void UpdateDestRanks(exchange &Exchange);
void UpdateSourceRanks(exchange &Exchange);

void UpdateDonorsSorted(exchange &Exchange);
void UpdateReceiversSorted(exchange &Exchange);

void UpdateSendInfo(exchange &Exchange);
void UpdateReceiveInfo(exchange &Exchange);

}

namespace core {

void CreateExchange(exchange &Exchange, const connectivity &Connectivity, logger &Logger,
  error_handler &ErrorHandler, profiler &Profiler) {

  Exchange.Connectivity_ = &Connectivity;

  Exchange.Logger_ = &Logger;
  Exchange.ErrorHandler_ = &ErrorHandler;

  Exchange.Profiler_ = &Profiler;
  AddProfilerTimer(Profiler, "Collect");
  AddProfilerTimer(Profiler, "Collect::MemAlloc");
  AddProfilerTimer(Profiler, "Collect::MPI");
  AddProfilerTimer(Profiler, "Collect::Pack");
  AddProfilerTimer(Profiler, "Collect::Reduce");
  AddProfilerTimer(Profiler, "SendRecv");
  AddProfilerTimer(Profiler, "SendRecv::MemAlloc");
  AddProfilerTimer(Profiler, "SendRecv::Pack");
  AddProfilerTimer(Profiler, "SendRecv::MPI");
  AddProfilerTimer(Profiler, "SendRecv::Unpack");
  AddProfilerTimer(Profiler, "Disperse");

  GetConnectivityDimension(Connectivity, Exchange.NumDims_);

  Exchange.Comm_ = core::GetConnectivityComm(Connectivity);

  MPI_Barrier(Exchange.Comm_);

  const grid_info *DonorGridInfo;
  GetConnectivityDonorGridInfo(Connectivity, DonorGridInfo);

  const grid_info *ReceiverGridInfo;
  GetConnectivityReceiverGridInfo(Connectivity, ReceiverGridInfo);

  const connectivity_d *Donors = nullptr;
  const grid *DonorGrid = nullptr;
  if (RankHasConnectivityDonorSide(Connectivity)) {
    GetConnectivityDonorSide(Connectivity, Donors);
    GetConnectivityDonorSideGrid(*Donors, DonorGrid);
  }

  const connectivity_r *Receivers = nullptr;
  const grid *ReceiverGrid = nullptr;
  if (RankHasConnectivityReceiverSide(Connectivity)) {
    GetConnectivityReceiverSide(Connectivity, Receivers);
    GetConnectivityReceiverSideGrid(*Receivers, ReceiverGrid);
  }

  range DonorGridGlobalRange;
  GetGridInfoGlobalRange(*DonorGridInfo, DonorGridGlobalRange);

  range DonorGridLocalRange = MakeEmptyRange(Exchange.NumDims_);
  if (DonorGrid) {
    GetGridLocalRange(*DonorGrid, DonorGridLocalRange);
  }

  CreatePartitionHash(Exchange.SourceHash_, Exchange.NumDims_, Exchange.Comm_, DonorGridGlobalRange,
    DonorGridLocalRange);

  range ReceiverGridGlobalRange;
  GetGridInfoGlobalRange(*ReceiverGridInfo, ReceiverGridGlobalRange);

  range ReceiverGridLocalRange = MakeEmptyRange(Exchange.NumDims_);
  if (ReceiverGrid) {
    GetGridLocalRange(*ReceiverGrid, ReceiverGridLocalRange);
  }

  CreatePartitionHash(Exchange.DestinationHash_, Exchange.NumDims_, Exchange.Comm_,
    ReceiverGridGlobalRange, ReceiverGridLocalRange);

  MPI_Barrier(Exchange.Comm_);

  core::LogStatus(*Exchange.Logger_, Exchange.Comm_.Rank() == 0, 0, "Created exchange %s.",
    Connectivity.Name_);

}

void DestroyExchange(exchange &Exchange) {

  MPI_Barrier(Exchange.Comm_);

  const connectivity &Connectivity = *Exchange.Connectivity_;

  Exchange.Sends_.Clear();
  Exchange.Recvs_.Clear();

  Exchange.NumRemoteDonorPoints_.Clear();

  Exchange.RemoteDonorPointsData_.Clear();
  Exchange.RemoteDonorPointCollectRecvsData_.Clear();
  Exchange.RemoteDonorPointCollectRecvBufferIndicesData_.Clear();

  Exchange.CollectSends_.Clear();
  Exchange.CollectRecvs_.Clear();

  DestroyPartitionHash(Exchange.SourceHash_);
  DestroyPartitionHash(Exchange.DestinationHash_);

  Exchange.DonorsSorted_.Clear();
  Exchange.ReceiversSorted_.Clear();

  Exchange.DonorDestRanks_.Clear();
  Exchange.ReceiverSourceRanks_.Clear();

  MPI_Barrier(Exchange.Comm_);

  LogStatus(*Exchange.Logger_, Exchange.Comm_.Rank() == 0, 0, "Destroyed exchange %s.",
    Connectivity.Name_);

  Exchange.Comm_.Reset();

}

void CreateExchangeInfo(exchange_info &Info, const exchange *Exchange, comm_view Comm) {

  bool IsLocal = Exchange != nullptr;
  bool IsRoot = false;
  if (IsLocal) {
    IsRoot = Exchange->Comm_.Rank() == 0;
  }

  const connectivity *Connectivity = nullptr;
  if (IsLocal) {
    Connectivity = Exchange->Connectivity_;
  }

  int RootRank;
  if (IsRoot) RootRank = Comm.Rank();
  core::BroadcastAnySource(&RootRank, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) {
    Info.DonorGridID_ = Connectivity->DonorGridID_;
    Info.ReceiverGridID_ = Connectivity->ReceiverGridID_;
    Info.NumDims_ = Connectivity->NumDims_;
  }
  MPI_Bcast(&Info.DonorGridID_, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info.ReceiverGridID_, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info.NumDims_, 1, MPI_INT, RootRank, Comm);

  int NameLength;
  if (IsRoot) NameLength = Connectivity->Name_.length();
  MPI_Bcast(&NameLength, 1, MPI_INT, RootRank, Comm);
  array<char> NameChars({NameLength});
  if (IsRoot) NameChars.Fill(Connectivity->Name_.begin());
  MPI_Bcast(NameChars.Data(), NameLength, MPI_CHAR, RootRank, Comm);
  Info.Name_.assign(NameChars.LinearBegin(), NameChars.LinearEnd());

  Info.RootRank_ = RootRank;

  Info.IsLocal_ = IsLocal;

}

void DestroyExchangeInfo(exchange_info &Info) {

  Info.Name_.clear();

}

}

bool RankHasExchangeDonorSide(const exchange &Exchange) {

  return RankHasConnectivityDonorSide(*Exchange.Connectivity_);

}

bool RankHasExchangeReceiverSide(const exchange &Exchange) {

  return RankHasConnectivityReceiverSide(*Exchange.Connectivity_);

}

namespace core {

void UpdateExchange(exchange &Exchange) {

  MPI_Barrier(Exchange.Comm_);

  const connectivity &Connectivity = *Exchange.Connectivity_;

  core::LogStatus(*Exchange.Logger_, Exchange.Comm_.Rank() == 0, 0, "Updating exchange %s...",
    Connectivity.Name_);

  const connectivity::edits *Edits;
  core::GetConnectivityEdits(Connectivity, Edits);

  bool NeedToUpdateDonorsSorted = false;
  bool NeedToUpdateReceiversSorted = false;
  bool NeedToUpdateDestSourceRanks = false;
  bool NeedToUpdateCollectInfo = false;
  bool NeedToUpdateSendRecvInfo = false;

  if (Edits->NumDonors_) {
    NeedToUpdateDonorsSorted = true;
    NeedToUpdateDestSourceRanks = true;
    NeedToUpdateCollectInfo = true;
    NeedToUpdateSendRecvInfo = true;
  }

  if (Edits->DonorExtents_) {
    NeedToUpdateCollectInfo = true;
  }

  if (Edits->DonorDestinations_) {
    NeedToUpdateDonorsSorted = true;
    NeedToUpdateDestSourceRanks = true;
    NeedToUpdateSendRecvInfo = true;
  }

  if (Edits->NumReceivers_) {
    NeedToUpdateReceiversSorted = true;
    NeedToUpdateDestSourceRanks = true;
    NeedToUpdateSendRecvInfo = true;
  }

  if (Edits->ReceiverSources_) {
    NeedToUpdateReceiversSorted = true;
    NeedToUpdateDestSourceRanks = true;
    NeedToUpdateSendRecvInfo = true;
  }

  if (NeedToUpdateDonorsSorted) UpdateDonorsSorted(Exchange);
  if (NeedToUpdateReceiversSorted) UpdateReceiversSorted(Exchange);

  if (NeedToUpdateDestSourceRanks) {
    UpdateDestRanks(Exchange);
    UpdateSourceRanks(Exchange);
  }

  if (NeedToUpdateCollectInfo) {
    UpdateCollectSendInfo(Exchange);
    UpdateCollectReceiveInfo(Exchange);
  }

  if (NeedToUpdateSendRecvInfo) {
    UpdateSendInfo(Exchange);
    UpdateReceiveInfo(Exchange);
  }

  MPI_Barrier(Exchange.Comm_);

  core::LogStatus(*Exchange.Logger_, Exchange.Comm_.Rank() == 0, 0, "Done updating exchange %s.",
    Connectivity.Name_);

}

}

namespace {

void UpdateCollectSendInfo(exchange &Exchange) {

  Exchange.CollectSends_.Clear();

  int NumDims = Exchange.NumDims_;
  const connectivity &Connectivity = *Exchange.Connectivity_;

  const connectivity_d *Donors;
  long long NumDonors = 0;
  if (RankHasConnectivityDonorSide(Connectivity)) {
    GetConnectivityDonorSide(Connectivity, Donors);
    GetConnectivityDonorSideCount(*Donors, NumDonors);
  }

  if (NumDonors > 0) {

    const grid *Grid;
    GetConnectivityDonorSideGrid(*Donors, Grid);

    range GlobalRange, LocalRange;
    GetGridGlobalRange(*Grid, GlobalRange);
    GetGridLocalRange(*Grid, LocalRange);

    const array<core::grid_neighbor> &GridNeighbors = core::GetGridNeighbors(*Grid);
    int NumNeighbors = GridNeighbors.Count();

    cart Cart;
    GetGridCart(*Grid, Cart);

    array<range> SendToNeighborRanges({NumNeighbors});
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      SendToNeighborRanges(iNeighbor) = MakeEmptyRange(NumDims);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      range DonorRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorRange.Begin(iDim) = Donors->Extents_(0,iDim,iDonor);
        DonorRange.End(iDim) = Donors->Extents_(1,iDim,iDonor);
      }
      bool AwayFromEdge = GlobalRange.Includes(DonorRange);
      for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
        if (AwayFromEdge) {
          bool Overlaps = RangesOverlap(GridNeighbors(iNeighbor).LocalRange, DonorRange);
          if (Overlaps) {
            SendToNeighborRanges(iNeighbor) = UnionRanges(SendToNeighborRanges(iNeighbor),
              IntersectRanges(LocalRange, DonorRange));
          }
        } else {
          bool Overlaps = false;
          for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
            for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
              for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
                elem<int,MAX_DIMS> Point = Cart.PeriodicAdjust({i,j,k});
                if (GridNeighbors(iNeighbor).LocalRange.Contains(Point)) {
                  Overlaps = true;
                  goto done_checking_for_overlap1;
                }
              }
            }
          }
          done_checking_for_overlap1:;
          if (Overlaps) {
            for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
              for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
                for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
                  elem<int,MAX_DIMS> Point = Cart.PeriodicAdjust({i,j,k});
                  if (LocalRange.Contains(Point)) {
                    ExtendRange(SendToNeighborRanges(iNeighbor), Point);
                  }
                }
              }
            }
          }
        }
      }
    }

    using range_indexer = indexer<long long, int, MAX_DIMS, array_layout::GRID>;
    array<range_indexer> SendToNeighborIndexers({NumNeighbors});
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      const range &SendToNeighborRange = SendToNeighborRanges(iNeighbor);
      SendToNeighborIndexers(iNeighbor) = range_indexer(SendToNeighborRange);
    }

    array<int> CollectSendIndexToNeighbor;
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!SendToNeighborRanges(iNeighbor).Empty()) {
        CollectSendIndexToNeighbor.Append(iNeighbor);
      }
    }
    int NumCollectSends = CollectSendIndexToNeighbor.Count();

    Exchange.CollectSends_.Resize({NumCollectSends});

    for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      int iNeighbor = CollectSendIndexToNeighbor(iCollectSend);
      Exchange.CollectSends_(iCollectSend).Rank = GridNeighbors(iNeighbor).Rank;
    }

    array<array<bool>> CollectSendMasks({NumCollectSends});
    for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      int iNeighbor = CollectSendIndexToNeighbor(iCollectSend);
      CollectSendMasks(iCollectSend).Resize({SendToNeighborRanges(iNeighbor).Count()}, false);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      range DonorRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorRange.Begin(iDim) = Donors->Extents_(0,iDim,iDonor);
        DonorRange.End(iDim) = Donors->Extents_(1,iDim,iDonor);
      }
      bool AwayFromEdge = GlobalRange.Includes(DonorRange);
      for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
        int iNeighbor = CollectSendIndexToNeighbor(iCollectSend);
        const range_indexer &Indexer = SendToNeighborIndexers(iNeighbor);
        if (AwayFromEdge) {
          bool Overlaps = RangesOverlap(GridNeighbors(iNeighbor).LocalRange, DonorRange);
          if (Overlaps) {
            range LocalDonorRange = IntersectRanges(LocalRange, DonorRange);
            for (int k = LocalDonorRange.Begin(2); k < LocalDonorRange.End(2); ++k) {
              for (int j = LocalDonorRange.Begin(1); j < LocalDonorRange.End(1); ++j) {
                for (int i = LocalDonorRange.Begin(0); i < LocalDonorRange.End(0); ++i) {
                  long long iPoint = Indexer.ToIndex(i,j,k);
                  CollectSendMasks(iCollectSend)(iPoint) = true;
                }
              }
            }
          }
        } else {
          bool Overlaps = false;
          for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
            for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
              for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
                elem<int,MAX_DIMS> Point = Cart.PeriodicAdjust({i,j,k});
                if (GridNeighbors(iNeighbor).LocalRange.Contains(Point)) {
                  Overlaps = true;
                  goto done_checking_for_overlap2;
                }
              }
            }
          }
          done_checking_for_overlap2:;
          if (Overlaps) {
            for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
              for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
                for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
                  elem<int,MAX_DIMS> Point = Cart.PeriodicAdjust({i,j,k});
                  if (LocalRange.Contains(Point)) {
                    long long iPoint = Indexer.ToIndex(Point);
                    CollectSendMasks(iCollectSend)(iPoint) = true;
                  }
                }
              }
            }
          }
        }
      }
    }

    for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      exchange::collect_send &CollectSend = Exchange.CollectSends_(iCollectSend);
      int iNeighbor = CollectSendIndexToNeighbor(iCollectSend);
      CollectSend.NumPoints = 0;
      for (long long iPoint = 0; iPoint < SendToNeighborRanges(iNeighbor).Count(); ++iPoint) {
        if (CollectSendMasks(iCollectSend)(iPoint)) {
          ++CollectSend.NumPoints;
        }
      }
    }

    for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      exchange::collect_send &CollectSend = Exchange.CollectSends_(iCollectSend);
      int iNeighbor = CollectSendIndexToNeighbor(iCollectSend);
      CollectSend.Points.Resize({{MAX_DIMS,CollectSend.NumPoints}});
      const range_indexer &Indexer = SendToNeighborIndexers(iNeighbor);
      long long iCollectSendPoint = 0;
      for (long long iPoint = 0; iPoint < SendToNeighborRanges(iNeighbor).Count(); ++iPoint) {
        if (CollectSendMasks(iCollectSend)(iPoint)) {
          elem<int,MAX_DIMS> Point = Indexer.ToTuple(iPoint);
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            CollectSend.Points(iDim,iCollectSendPoint) = Point[iDim];
          }
          ++iCollectSendPoint;
        }
      }
    }

  }

}

void UpdateCollectReceiveInfo(exchange &Exchange) {

  Exchange.CollectRecvs_.Clear();

  Exchange.RemoteDonorPoints_.Clear();
  Exchange.RemoteDonorPointsData_.Clear();
  Exchange.RemoteDonorPointCollectRecvs_.Clear();
  Exchange.RemoteDonorPointCollectRecvsData_.Clear();
  Exchange.RemoteDonorPointCollectRecvBufferIndices_.Clear();
  Exchange.RemoteDonorPointCollectRecvBufferIndicesData_.Clear();

  int NumDims = Exchange.NumDims_;
  const connectivity &Connectivity = *Exchange.Connectivity_;

  const connectivity_d *Donors;
  long long NumDonors = 0;
  int MaxSize = 0;
  if (RankHasConnectivityDonorSide(Connectivity)) {
    GetConnectivityDonorSide(Connectivity, Donors);
    GetConnectivityDonorSideCount(*Donors, NumDonors);
    GetConnectivityDonorSideMaxSize(*Donors, MaxSize);
  }

  if (NumDonors > 0) {

    const grid *Grid;
    GetConnectivityDonorSideGrid(*Donors, Grid);

    range GlobalRange, LocalRange;
    GetGridGlobalRange(*Grid, GlobalRange);
    GetGridLocalRange(*Grid, LocalRange);

    const array<core::grid_neighbor> &GridNeighbors = core::GetGridNeighbors(*Grid);
    int NumNeighbors = GridNeighbors.Count();

    cart Cart;
    GetGridCart(*Grid, Cart);

    array<range> RecvFromNeighborRanges({NumNeighbors});
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      RecvFromNeighborRanges(iNeighbor) = MakeEmptyRange(NumDims);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      range DonorRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorRange.Begin(iDim) = Donors->Extents_(0,iDim,iDonor);
        DonorRange.End(iDim) = Donors->Extents_(1,iDim,iDonor);
      }
      bool AwayFromEdge = GlobalRange.Includes(DonorRange);
      for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
        if (AwayFromEdge) {
          RecvFromNeighborRanges(iNeighbor) = UnionRanges(RecvFromNeighborRanges(iNeighbor),
            IntersectRanges(GridNeighbors(iNeighbor).LocalRange, DonorRange));
        } else {
          for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
            for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
              for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
                elem<int,MAX_DIMS> Point = Cart.PeriodicAdjust({i,j,k});
                if (GridNeighbors(iNeighbor).LocalRange.Contains(Point)) {
                  ExtendRange(RecvFromNeighborRanges(iNeighbor), Point);
                }
              }
            }
          }
        }
      }
    }

    using range_indexer = indexer<long long, int, MAX_DIMS, array_layout::GRID>;
    array<range_indexer> RecvFromNeighborIndexers({NumNeighbors});
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      const range &RecvFromNeighborRange = RecvFromNeighborRanges(iNeighbor);
      RecvFromNeighborIndexers(iNeighbor) = range_indexer(RecvFromNeighborRange);
    }

    array<int> CollectRecvIndexToNeighbor;
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!RecvFromNeighborRanges(iNeighbor).Empty()) {
        CollectRecvIndexToNeighbor.Append(iNeighbor);
      }
    }
    int NumCollectRecvs = CollectRecvIndexToNeighbor.Count();

    Exchange.CollectRecvs_.Resize({NumCollectRecvs});

    for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      int iNeighbor = CollectRecvIndexToNeighbor(iCollectRecv);
      Exchange.CollectRecvs_(iCollectRecv).Rank = GridNeighbors(iNeighbor).Rank;
    }

    array<array<bool>> CollectRecvMasks({NumCollectRecvs});
    for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      int iNeighbor = CollectRecvIndexToNeighbor(iCollectRecv);
      CollectRecvMasks(iCollectRecv).Resize({RecvFromNeighborRanges(iNeighbor).Count()}, false);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      range DonorRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorRange.Begin(iDim) = Donors->Extents_(0,iDim,iDonor);
        DonorRange.End(iDim) = Donors->Extents_(1,iDim,iDonor);
      }
      bool AwayFromEdge = GlobalRange.Includes(DonorRange);
      for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
        int iNeighbor = CollectRecvIndexToNeighbor(iCollectRecv);
        const range_indexer &Indexer = RecvFromNeighborIndexers(iNeighbor);
        if (AwayFromEdge) {
          range RemoteDonorRange = IntersectRanges(GridNeighbors(iNeighbor).LocalRange, DonorRange);
          for (int k = RemoteDonorRange.Begin(2); k < RemoteDonorRange.End(2); ++k) {
            for (int j = RemoteDonorRange.Begin(1); j < RemoteDonorRange.End(1); ++j) {
              for (int i = RemoteDonorRange.Begin(0); i < RemoteDonorRange.End(0); ++i) {
                long long iPoint = Indexer.ToIndex(i,j,k);
                CollectRecvMasks(iCollectRecv)(iPoint) = true;
              }
            }
          }
        } else {
          for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
            for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
              for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
                elem<int,MAX_DIMS> Point = Cart.PeriodicAdjust({i,j,k});
                if (GridNeighbors(iNeighbor).LocalRange.Contains(Point)) {
                  long long iPoint = Indexer.ToIndex(Point);
                  CollectRecvMasks(iCollectRecv)(iPoint) = true;
                }
              }
            }
          }
        }
      }
    }

    for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      exchange::collect_recv &CollectRecv = Exchange.CollectRecvs_(iCollectRecv);
      int iNeighbor = CollectRecvIndexToNeighbor(iCollectRecv);
      CollectRecv.NumPoints = 0;
      for (long long iPoint = 0; iPoint < RecvFromNeighborRanges(iNeighbor).Count(); ++iPoint) {
        if (CollectRecvMasks(iCollectRecv)(iPoint)) {
          ++CollectRecv.NumPoints;
        }
      }
    }

    Exchange.NumRemoteDonorPoints_.Resize({NumDonors}, 0);
    Exchange.RemoteDonorPoints_.Resize({NumDonors}, nullptr);
    Exchange.RemoteDonorPointCollectRecvs_.Resize({NumDonors}, nullptr);
    Exchange.RemoteDonorPointCollectRecvBufferIndices_.Resize({NumDonors}, nullptr);

    long long TotalRemoteDonorPoints = 0;
    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      range DonorRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorRange.Begin(iDim) = Donors->Extents_(0,iDim,iDonor);
        DonorRange.End(iDim) = Donors->Extents_(1,iDim,iDonor);
      }
      bool AwayFromEdge = GlobalRange.Includes(DonorRange);
      int NumRemoteDonorPoints;
      if (AwayFromEdge) {
        range LocalDonorRange = IntersectRanges(LocalRange, DonorRange);
        NumRemoteDonorPoints = DonorRange.Count() - LocalDonorRange.Count();
      } else {
        NumRemoteDonorPoints = 0;
        for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
          for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
            for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
              elem<int,MAX_DIMS> Point = Cart.PeriodicAdjust({i,j,k});
              if (!LocalRange.Contains(Point)) {
                ++NumRemoteDonorPoints;
              }
            }
          }
        }
      }
      Exchange.NumRemoteDonorPoints_(iDonor) = NumRemoteDonorPoints;
      TotalRemoteDonorPoints += NumRemoteDonorPoints;
    }

    Exchange.RemoteDonorPointsData_.Resize({TotalRemoteDonorPoints});
    Exchange.RemoteDonorPointCollectRecvsData_.Resize({TotalRemoteDonorPoints});
    Exchange.RemoteDonorPointCollectRecvBufferIndicesData_.Resize({TotalRemoteDonorPoints});

    long long Offset = 0;
    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      long long NumRemoteDonorPoints = Exchange.NumRemoteDonorPoints_(iDonor);
      Exchange.RemoteDonorPoints_(iDonor) = Exchange.RemoteDonorPointsData_.Data(Offset);
      Exchange.RemoteDonorPointCollectRecvs_(iDonor) =
        Exchange.RemoteDonorPointCollectRecvsData_.Data(Offset);
      Exchange.RemoteDonorPointCollectRecvBufferIndices_(iDonor) =
        Exchange.RemoteDonorPointCollectRecvBufferIndicesData_.Data(Offset);
      Offset += NumRemoteDonorPoints;
    }

    array<array<long long>> CollectRecvBufferIndices({NumCollectRecvs});
    for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      int iNeighbor = CollectRecvIndexToNeighbor(iCollectRecv);
      long long NumPoints = RecvFromNeighborRanges(iNeighbor).Count();
      CollectRecvBufferIndices(iCollectRecv).Resize({NumPoints}, -1);
      long long iRemotePoint = 0;
      for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
        if (CollectRecvMasks(iCollectRecv)(iPoint)) {
          CollectRecvBufferIndices(iCollectRecv)(iPoint) = iRemotePoint;
          ++iRemotePoint;
        }
      }
    }

    int MaxPointsInCell = 1;
    for (int iDim = 0; iDim <= NumDims; ++iDim) {
      MaxPointsInCell *= MaxSize;
    }
    array<int> CellCollectRecvs({MaxPointsInCell});
    array<long long> CellCollectRecvBufferIndices({MaxPointsInCell});

    using donor_indexer = indexer<int, int, MAX_DIMS, array_layout::GRID>;

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      for (int iPointInCell = 0; iPointInCell < MaxPointsInCell; ++iPointInCell) {
        CellCollectRecvs(iPointInCell) = -1;
        CellCollectRecvBufferIndices(iPointInCell) = -1;
      }
      range DonorRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorRange.Begin(iDim) = Donors->Extents_(0,iDim,iDonor);
        DonorRange.End(iDim) = Donors->Extents_(1,iDim,iDonor);
      }
      donor_indexer DonorIndexer(DonorRange);
      bool AwayFromEdge = GlobalRange.Includes(DonorRange);
      for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
        int iNeighbor = CollectRecvIndexToNeighbor(iCollectRecv);
        const range_indexer &RecvFromNeighborIndexer = RecvFromNeighborIndexers(iNeighbor);
        if (AwayFromEdge) {
          range RemoteDonorRange = IntersectRanges(GridNeighbors(iNeighbor).LocalRange, DonorRange);
          for (int k = RemoteDonorRange.Begin(2); k < RemoteDonorRange.End(2); ++k) {
            for (int j = RemoteDonorRange.Begin(1); j < RemoteDonorRange.End(1); ++j) {
              for (int i = RemoteDonorRange.Begin(0); i < RemoteDonorRange.End(0); ++i) {
                int iPointInCell = DonorIndexer.ToIndex(i,j,k);
                long long iPoint = RecvFromNeighborIndexer.ToIndex(i,j,k);
                CellCollectRecvs(iPointInCell) = iCollectRecv;
                CellCollectRecvBufferIndices(iPointInCell) = CollectRecvBufferIndices(iCollectRecv)
                  (iPoint);
              }
            }
          }
        } else {
          int iPointInCell = 0;
          for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
            for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
              for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
                elem<int,MAX_DIMS> Point = Cart.PeriodicAdjust({i,j,k});
                if (GridNeighbors(iNeighbor).LocalRange.Contains(Point)) {
                  long long iPoint = RecvFromNeighborIndexer.ToIndex(Point);
                  CellCollectRecvs(iPointInCell) = iCollectRecv;
                  CellCollectRecvBufferIndices(iPointInCell) = CollectRecvBufferIndices
                    (iCollectRecv)(iPoint);
                }
                ++iPointInCell;
              }
            }
          }
        }
      }
      int iRemoteDonorPoint = 0;
      int iPointInCell = 0;
      for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
        for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
          for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
            if (CellCollectRecvs(iPointInCell) >= 0) {
              Exchange.RemoteDonorPoints_(iDonor)[iRemoteDonorPoint] = iPointInCell;
              Exchange.RemoteDonorPointCollectRecvs_(iDonor)[iRemoteDonorPoint] =
                CellCollectRecvs(iPointInCell);
              Exchange.RemoteDonorPointCollectRecvBufferIndices_(iDonor)[iRemoteDonorPoint]
                = CellCollectRecvBufferIndices(iPointInCell);
              ++iRemoteDonorPoint;
            }
            ++iPointInCell;
          }
        }
      }
    }

  }

}

void UpdateDonorsSorted(exchange &Exchange) {

  Exchange.DonorsSorted_.Clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  if (RankHasConnectivityDonorSide(Connectivity)) {

    const connectivity_d *Donors;
    GetConnectivityDonorSide(Connectivity, Donors);

    long long NumDonors;
    GetConnectivityDonorSideCount(*Donors, NumDonors);

    if (NumDonors > 0) {

      Exchange.DonorsSorted_.Resize({NumDonors});

      const grid_info *ReceiverGridInfo;
      GetConnectivityReceiverGridInfo(Connectivity, ReceiverGridInfo);

      range ReceiverGridGlobalRange;
      GetGridInfoGlobalRange(*ReceiverGridInfo, ReceiverGridGlobalRange);

      using range_indexer = indexer<long long, int, MAX_DIMS, array_layout::GRID>;
      range_indexer ReceiverGridGlobalIndexer(ReceiverGridGlobalRange);

      array<long long> DestinationIndices({NumDonors});

      for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
        elem<int,MAX_DIMS> DestinationPoint = {
          Donors->Destinations_(0,iDonor),
          Donors->Destinations_(1,iDonor),
          Donors->Destinations_(2,iDonor)
        };
        DestinationIndices(iDonor) = ReceiverGridGlobalIndexer.ToIndex(DestinationPoint);
      }

      bool Sorted = true;

      // Check if they're already sorted
      long long PrevIndex = 0;
      for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
        if (DestinationIndices(iDonor) < PrevIndex) {
          Sorted = false;
          break;
        }
        PrevIndex = DestinationIndices(iDonor);
      }

      if (Sorted) {
        for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
          Exchange.DonorsSorted_(iDonor) = iDonor;
        }
      } else {
        core::SortPermutation(DestinationIndices, Exchange.DonorsSorted_);
      }

    }

  }

}

void UpdateReceiversSorted(exchange &Exchange) {

  Exchange.ReceiversSorted_.Clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  if (RankHasConnectivityReceiverSide(Connectivity)) {

    const connectivity_r *Receivers;
    GetConnectivityReceiverSide(Connectivity, Receivers);

    long long NumReceivers = 0;
    GetConnectivityReceiverSideCount(*Receivers, NumReceivers);

    if (NumReceivers > 0) {

      Exchange.ReceiversSorted_.Resize({NumReceivers});

      const grid *ReceiverGrid;
      GetConnectivityReceiverSideGrid(*Receivers, ReceiverGrid);

      range GlobalRange;
      GetGridGlobalRange(*ReceiverGrid, GlobalRange);

      using range_indexer = indexer<long long, int, MAX_DIMS, array_layout::GRID>;
      range_indexer GlobalIndexer(GlobalRange);

      array<long long> PointIndices({NumReceivers});

      for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
        elem<int,MAX_DIMS> Point = {
          Receivers->Points_(0,iReceiver),
          Receivers->Points_(1,iReceiver),
          Receivers->Points_(2,iReceiver)
        };
        PointIndices(iReceiver) = GlobalIndexer.ToIndex(Point);
      }

      bool Sorted = true;

      // Check if they're already sorted
      long long PrevIndex = 0;
      for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
        if (PointIndices(iReceiver) < PrevIndex) {
          Sorted = false;
          break;
        }
        PrevIndex = PointIndices(iReceiver);
      }

      if (Sorted) {
        for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
          Exchange.ReceiversSorted_(iReceiver) = iReceiver;
        }
      } else {
        core::SortPermutation(PointIndices, Exchange.ReceiversSorted_);
      }

    }

  }

}

void UpdateDestRanks(exchange &Exchange) {

  Exchange.DonorDestRanks_.Clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  bool DonorGridIsLocal = RankHasConnectivityDonorSide(Connectivity);

  long long NumDonors = 0;

  const connectivity_d *Donors;
  if (DonorGridIsLocal) {
    GetConnectivityDonorSide(Connectivity, Donors);
    GetConnectivityDonorSideCount(*Donors, NumDonors);
  }

  Exchange.DonorDestRanks_.Resize({NumDonors}, -1);

  std::map<int, core::partition_bin> Bins;

  array<int> DestinationBinIndices;
  if (DonorGridIsLocal) {
    DestinationBinIndices.Resize({NumDonors});
    core::MapToPartitionBins(Exchange.DestinationHash_, Donors->Destinations_,
      DestinationBinIndices);
    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int BinIndex = DestinationBinIndices(iDonor);
      auto Iter = Bins.lower_bound(BinIndex);
      if (Iter == Bins.end() || Iter->first > BinIndex) {
        Bins.emplace_hint(Iter, BinIndex, core::partition_bin());
      }
    }
  }

  core::RetrievePartitionBins(Exchange.DestinationHash_, Bins);

  if (DonorGridIsLocal) {
    core::FindPartitions(Exchange.DestinationHash_, Bins, Donors->Destinations_,
      DestinationBinIndices, Exchange.DonorDestRanks_);
  }

}

void UpdateSourceRanks(exchange &Exchange) {

  Exchange.ReceiverSourceRanks_.Clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  bool ReceiverGridIsLocal = RankHasConnectivityReceiverSide(Connectivity);

  long long NumReceivers = 0;

  const connectivity_r *Receivers;
  if (ReceiverGridIsLocal) {
    GetConnectivityReceiverSide(Connectivity, Receivers);
    GetConnectivityReceiverSideCount(*Receivers, NumReceivers);
  }

  Exchange.ReceiverSourceRanks_.Resize({NumReceivers}, -1);

  std::map<int, core::partition_bin> Bins;

  array<int> SourceBinIndices;
  if (ReceiverGridIsLocal) {
    SourceBinIndices.Resize({NumReceivers});
    core::MapToPartitionBins(Exchange.SourceHash_, Receivers->Sources_, SourceBinIndices);
    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      int BinIndex = SourceBinIndices(iReceiver);
      auto Iter = Bins.lower_bound(BinIndex);
      if (Iter == Bins.end() || Iter->first > BinIndex) {
        Bins.emplace_hint(Iter, BinIndex, core::partition_bin());
      }
    }
  }

  core::RetrievePartitionBins(Exchange.SourceHash_, Bins);

  if (ReceiverGridIsLocal) {
    core::FindPartitions(Exchange.SourceHash_, Bins, Receivers->Sources_, SourceBinIndices,
      Exchange.ReceiverSourceRanks_);
  }

}

void UpdateSendInfo(exchange &Exchange) {

  Exchange.Sends_.Clear();

  Exchange.DonorSendIndices_.Clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  const connectivity_d *Donors;
  long long NumDonors = 0;
  cart Cart;
  range LocalRange;
  if (RankHasConnectivityDonorSide(Connectivity)) {
    GetConnectivityDonorSide(Connectivity, Donors);
    GetConnectivityDonorSideCount(*Donors, NumDonors);
    const grid *DonorGrid;
    GetConnectivityDonorSideGrid(*Donors, DonorGrid);
    GetGridCart(*DonorGrid, Cart);
    GetGridLocalRange(*DonorGrid, LocalRange);
  }

  if (NumDonors > 0) {

    Exchange.DonorSendIndices_.Resize({NumDonors}, -1);

    array<bool> DonorCommunicates({NumDonors});

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      bool Communicates = Exchange.DonorDestRanks_(iDonor) >= 0;
      if (Communicates) {
        elem<int,MAX_DIMS> DonorCell = Cart.PeriodicAdjust({
          Donors->Extents_(0,0,iDonor),
          Donors->Extents_(0,1,iDonor),
          Donors->Extents_(0,2,iDonor)
        });
        Communicates = LocalRange.Contains(DonorCell);
      }
      DonorCommunicates(iDonor) = Communicates;
    }

    std::map<int, long long> SendCounts;

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      if (DonorCommunicates(iDonor)) {
        int DestRank = Exchange.DonorDestRanks_(iDonor);
        auto Iter = SendCounts.lower_bound(DestRank);
        if (Iter != SendCounts.end() && Iter->first == DestRank) {
          ++Iter->second;
        } else {
          SendCounts.emplace_hint(Iter, DestRank, 1);
        }
      }
    }

    int NumSends = SendCounts.size();

    for (auto &Pair : SendCounts) {
      exchange::send &Send = Exchange.Sends_.Append();
      Send.Rank = Pair.first;
      Send.Count = Pair.second;
    }

    SendCounts.clear();

    std::map<int, int> RankToSendIndex;

    for (int iSend = 0; iSend < NumSends; ++iSend) {
      RankToSendIndex.emplace(Exchange.Sends_(iSend).Rank, iSend);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      if (DonorCommunicates(iDonor)) {
        int DestRank = Exchange.DonorDestRanks_(iDonor);
        Exchange.DonorSendIndices_(iDonor) = RankToSendIndex[DestRank];
      }
    }

  }

}

void UpdateReceiveInfo(exchange &Exchange) {

  Exchange.Recvs_.Clear();

  Exchange.ReceiverRecvIndices_.Clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  const connectivity_r *Receivers;
  long long NumReceivers = 0;
  if (RankHasConnectivityReceiverSide(Connectivity)) {
    GetConnectivityReceiverSide(Connectivity, Receivers);
    GetConnectivityReceiverSideCount(*Receivers, NumReceivers);
  }

  if (NumReceivers > 0) {

    Exchange.ReceiverRecvIndices_.Resize({NumReceivers}, -1);

    std::map<int, long long> RecvCounts;

    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      if (Exchange.ReceiverSourceRanks_(iReceiver) >= 0) {
        int SourceRank = Exchange.ReceiverSourceRanks_(iReceiver);
        auto Iter = RecvCounts.lower_bound(SourceRank);
        if (Iter != RecvCounts.end() && Iter->first == SourceRank) {
          ++Iter->second;
        } else {
          RecvCounts.emplace_hint(Iter, SourceRank, 1);
        }
      }
    }

    int NumRecvs = RecvCounts.size();

    for (auto &Pair : RecvCounts) {
      exchange::recv &Recv = Exchange.Recvs_.Append();
      Recv.Rank = Pair.first;
      Recv.Count = Pair.second;
    }

    RecvCounts.clear();

    std::map<int, int> RankToRecvIndex;

    for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
      RankToRecvIndex.emplace(Exchange.Recvs_(iRecv).Rank, iRecv);
    }

    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      int SourceRank = Exchange.ReceiverSourceRanks_(iReceiver);
      if (SourceRank >= 0) {
        Exchange.ReceiverRecvIndices_(iReceiver) = RankToRecvIndex[SourceRank];
      }
    }

  }

}

}

namespace core {

void Collect(const exchange &Exchange, data_type ValueType, int Count, collect_op CollectOp,
  const range &GridValuesRange, array_layout GridValuesLayout, const void * const *GridValues,
  void **DonorValues) {

  collect Collect_ = MakeCollect(CollectOp, ValueType, GridValuesLayout);

  Collect_.Initialize(Exchange, Count, GridValuesRange);
  Collect_.Collect(GridValues, DonorValues);

}

request Send(const exchange &Exchange, data_type ValueType, int Count, const void * const
  *DonorValues, int Tag) {

  send Send_ = MakeSend(ValueType);

  Send_.Initialize(Exchange, Count, Tag);

  return Send_.Send(DonorValues);

}

request Receive(const exchange &Exchange, data_type ValueType, int Count, void **ReceiverValues,
  int Tag) {

  recv Recv = MakeRecv(ValueType);

  Recv.Initialize(Exchange, Count, Tag);

  return Recv.Recv(ReceiverValues);

}

void Disperse(const exchange &Exchange, data_type ValueType, int Count, disperse_op DisperseOp,
  const void * const *ReceiverValues, const range &GridValuesRange, array_layout GridValuesLayout,
  void **GridValues) {

  disperse Disperse_ = MakeDisperse(DisperseOp, ValueType, GridValuesLayout);

  Disperse_.Initialize(Exchange, Count, GridValuesRange);
  Disperse_.Disperse(ReceiverValues, GridValues);

}

}

void GetExchangeInfoDonorGridID(const exchange_info &Info, int &DonorGridID) {

  DonorGridID = Info.DonorGridID_;

}

void GetExchangeInfoReceiverGridID(const exchange_info &Info, int &ReceiverGridID) {

  ReceiverGridID = Info.ReceiverGridID_;

}

void GetExchangeInfoName(const exchange_info &Info, std::string &Name) {

  Name = Info.Name_;

}

void GetExchangeInfoDimension(const exchange_info &Info, int &NumDims) {

  NumDims = Info.NumDims_;

}

void GetExchangeInfoRootRank(const exchange_info &Info, int &RootRank) {

  RootRank = Info.RootRank_;

}

void GetExchangeInfoIsLocal(const exchange_info &Info, bool &IsLocal) {

  IsLocal = Info.IsLocal_;

}

}
