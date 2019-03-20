// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Exchange.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Request.hpp"
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

class collect {

public:

  collect() = default;
  collect(const collect &) = delete;
  collect(collect &&) = default;

  template <typename T> explicit collect(T &&Collect):
    Collect_(new model<T>(std::forward<T>(Collect)))
  {}

  template <typename T> collect &operator=(T &&Collect) {
    Collect_.reset(new model<T>(std::forward<T>(Collect)));
    return *this;
  }

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {
    Collect_->Initialize(Exchange, Count, GridValuesRange);
  }

  void Collect(const void * const *GridValues, void **DonorValues) {
    Collect_->Collect(GridValues, DonorValues);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) = 0;
    virtual void Collect(const void * const *GridValues, void **DonorValues) = 0;
  };

  template <typename T> class model : public concept {
  public:
    using value_type = typename T::value_type;
    explicit model(T Collect):
      Collect_(std::move(Collect))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange)
      override {
      Collect_.Initialize(Exchange, Count, GridValuesRange);
      GridValuesRange_ = GridValuesRange;
      const connectivity_d *Donors;
      GetConnectivityDonorSide(*Exchange.Connectivity_, Donors);
      GetConnectivityDonorSideCount(*Donors, NumDonors_);
      GridValues_.Resize({Count});
      DonorValues_.Resize({Count});
    }
    virtual void Collect(const void * const *GridValuesVoid, void **DonorValuesVoid) override {
      OVK_DEBUG_ASSERT(GridValuesVoid || GridValues_.Count() == 0, "Invalid grid values pointer.");
      OVK_DEBUG_ASSERT(DonorValuesVoid || DonorValues_.Count() == 0, "Invalid donor values pointer.");
      for (int iCount = 0; iCount < GridValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(GridValuesVoid[iCount] || NumDonors_ == 0, "Invalid grid values pointer.");
        GridValues_(iCount) = {static_cast<const value_type *>(GridValuesVoid[iCount]),
          {GridValuesRange_.Count()}};
      }
      for (int iCount = 0; iCount < DonorValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(DonorValuesVoid[iCount] || NumDonors_ == 0, "Invalid donor values "
          "pointer.");
        DonorValues_(iCount) = {static_cast<value_type *>(DonorValuesVoid[iCount]),
          {NumDonors_}};
      }
      Collect_.Collect(GridValues_, DonorValues_);
    }
  private:
    T Collect_;
    range GridValuesRange_;
    long long NumDonors_;
    array<array_view<const value_type>> GridValues_;
    array<array_view<value_type>> DonorValues_;
  };

  std::unique_ptr<concept> Collect_;

};

template <typename T, array_layout Layout> class collect_base {

  using mpi_value_type = core::mpi_compatible_type<T>;

public:

  using value_type = T;

  collect_base() = default;
  collect_base(const collect_base &) = delete;
  collect_base(collect_base &&) = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    int NumDims = Exchange.NumDims_;
    Count_ = Count;

    const connectivity &Connectivity = *Exchange.Connectivity_;

    GetConnectivityDonorSide(Connectivity, Donors_);
    const connectivity_d &Donors = *Donors_;

    GetConnectivityDonorSideGrid(Donors, Grid_);
    const grid &Grid = *Grid_;

    Profiler_ = Exchange.Profiler_;
    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");

    GetGridCart(Grid, Cart_);
    GetGridGlobalRange(Grid, GlobalRange_);
    GetGridLocalRange(Grid, LocalRange_);

    OVK_DEBUG_ASSERT(GridValuesRange.Includes(LocalRange_), "Invalid grid values range.");

    int MaxSize;
    GetConnectivityDonorSideCount(Donors, NumDonors_);
    GetConnectivityDonorSideMaxSize(Donors, MaxSize);

    MaxPointsInCell_ = 1;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      MaxPointsInCell_ *= MaxSize;
    }

    GridValuesIndexer_ = range_indexer(GridValuesRange);

    Sends_ = Exchange.CollectSends_;
    Recvs_ = Exchange.CollectRecvs_;

    int NumSends = Sends_.Count();
    int NumRecvs = Recvs_.Count();

    core::StartProfile(*Profiler_, MemAllocTime);

    SendBuffers_.Resize({NumSends});
    for (int iSend = 0; iSend < NumSends; ++iSend) {
      SendBuffers_(iSend).Resize({{Count_,Exchange.CollectSends_(iSend).NumPoints}});
    }

    if (!std::is_same<value_type, mpi_value_type>::value) {
      RecvBuffers_.Resize({NumRecvs});
      for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
        RecvBuffers_(iRecv).Resize({{Count_,Exchange.CollectRecvs_(iRecv).NumPoints}});
      }
    }

    Requests_.Reserve(NumSends+NumRecvs);

    core::EndProfile(*Profiler_, MemAllocTime);

    NumRemoteDonorPoints_ = Exchange.NumRemoteDonorPoints_;
    RemoteDonorPoints_ = Exchange.RemoteDonorPoints_;
    RemoteDonorPointCollectRecvs_ = Exchange.RemoteDonorPointCollectRecvs_;
    RemoteDonorPointCollectRecvBufferIndices_ = Exchange.RemoteDonorPointCollectRecvBufferIndices_;

    LocalDonorPointIndices_.Resize({MaxPointsInCell_});
    LocalDonorPointGridValuesIndices_.Resize({MaxPointsInCell_});

  }

protected:

  using range_indexer = indexer<long long, int, MAX_DIMS, Layout>;
  using donor_indexer = indexer<int, int, MAX_DIMS, Layout>;

  const connectivity_d *Donors_;
  const grid *Grid_;
  core::profiler *Profiler_;
  int Count_;
  long long NumDonors_;
  int MaxPointsInCell_;
  range_indexer GridValuesIndexer_;

  void AllocateRemoteDonorValues(array<array<value_type,2>> &RemoteDonorValues) {

    int NumRecvs = Recvs_.Count();

    RemoteDonorValues.Resize({NumRecvs});
    for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
      RemoteDonorValues(iRecv).Resize({{Count_,Recvs_(iRecv).NumPoints}});
    }

  }

  void RetrieveRemoteDonorValues(array_view<array_view<const value_type>> GridValues, array<array<
    value_type,2>> &RemoteDonorValues) {

    MPI_Comm Comm;
    GetGridComm(*Grid_, Comm);

    MPI_Datatype MPIDataType = core::GetMPIDataType<mpi_value_type>();

    int MPITime = core::GetProfilerTimerID(*Profiler_, "Collect::MPI");
    int PackTime = core::GetProfilerTimerID(*Profiler_, "Collect::Pack");

    core::StartProfile(*Profiler_, MPITime);

    if (std::is_same<value_type, mpi_value_type>::value) {
      for (int iRecv = 0; iRecv < Recvs_.Count(); ++iRecv) {
        const exchange::collect_recv &Recv = Recvs_(iRecv);
        MPI_Request &Request = Requests_.Append();
        MPI_Irecv(RemoteDonorValues(iRecv).Data(), Count_*Recv.NumPoints, MPIDataType, Recv.Rank, 0,
          Comm, &Request);
      }
    } else {
      for (int iRecv = 0; iRecv < Recvs_.Count(); ++iRecv) {
        const exchange::collect_recv &Recv = Recvs_(iRecv);
        MPI_Request &Request = Requests_.Append();
        MPI_Irecv(RecvBuffers_(iRecv).Data(), Count_*Recv.NumPoints, MPIDataType, Recv.Rank, 0,
          Comm, &Request);
      }
    }

    core::EndProfile(*Profiler_, MPITime);
    core::StartProfile(*Profiler_, PackTime);

    for (int iSend = 0; iSend < Sends_.Count(); ++iSend) {
      const exchange::collect_send &Send = Sends_(iSend);
      for (long long iSendPoint = 0; iSendPoint < Send.NumPoints; ++iSendPoint) {
        elem<int,MAX_DIMS> Point = {
          Send.Points(0,iSendPoint),
          Send.Points(1,iSendPoint),
          Send.Points(2,iSendPoint)
        };
        long long iGridPoint = GridValuesIndexer_.ToIndex(Point);
        for (int iCount = 0; iCount < Count_; ++iCount) {
          SendBuffers_(iSend)(iCount,iSendPoint) = mpi_value_type(GridValues(iCount)(iGridPoint));
        }
      }
    }

    core::EndProfile(*Profiler_, PackTime);
    core::StartProfile(*Profiler_, MPITime);

    for (int iSend = 0; iSend < Sends_.Count(); ++iSend) {
      const exchange::collect_send &Send = Sends_(iSend);
      MPI_Request &Request = Requests_.Append();
      MPI_Isend(SendBuffers_(iSend).Data(), Count_*Send.NumPoints, MPIDataType, Send.Rank, 0, Comm,
        &Request);
    }

    MPI_Waitall(Requests_.Count(), Requests_.Data(), MPI_STATUSES_IGNORE);

    if (!std::is_same<value_type, mpi_value_type>::value) {
      for (int iRecv = 0; iRecv < RecvBuffers_.Count(); ++iRecv) {
        for (long long iBuffer = 0; iBuffer < RecvBuffers_(iRecv).Count(); ++iBuffer) {
          RemoteDonorValues(iRecv)[iBuffer] = value_type(RecvBuffers_(iRecv)[iBuffer]);
        }
      }
    }

    core::EndProfile(*Profiler_, MPITime);

    Requests_.Clear();

  }

  void AssembleDonorPointValues(array_view<array_view<const value_type>> GridValues, const
    array<array<value_type,2>> &RemoteDonorValues, long long iDonor, elem<int,MAX_DIMS>
    &DonorSize, array_view<value_type,2> DonorPointValues) {

    const connectivity_d &Donors = *Donors_;

    range DonorRange;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      DonorRange.Begin(iDim) = Donors.Extents_(0,iDim,iDonor);
      DonorRange.End(iDim) = Donors.Extents_(1,iDim,iDonor);
    }

    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      DonorSize[iDim] = DonorRange.Size(iDim);
    }

    using donor_indexer = indexer<int, int, MAX_DIMS, array_layout::GRID>;
    donor_indexer DonorIndexer(DonorRange);

    bool AwayFromEdge = GlobalRange_.Includes(DonorRange);

    // Fill in the local data
    int NumLocalDonorPoints;
    if (AwayFromEdge) {
      range LocalDonorRange = IntersectRanges(LocalRange_, DonorRange);
      int iLocalDonorPoint = 0;
      for (int k = LocalDonorRange.Begin(2); k < LocalDonorRange.End(2); ++k) {
        for (int j = LocalDonorRange.Begin(1); j < LocalDonorRange.End(1); ++j) {
          for (int i = LocalDonorRange.Begin(0); i < LocalDonorRange.End(0); ++i) {
            LocalDonorPointIndices_(iLocalDonorPoint) = DonorIndexer.ToIndex(i,j,k);
            LocalDonorPointGridValuesIndices_(iLocalDonorPoint) = GridValuesIndexer_.ToIndex(i,j,k);
            ++iLocalDonorPoint;
          }
        }
      }
      NumLocalDonorPoints = iLocalDonorPoint;
    } else {
      int iLocalDonorPoint = 0;
      for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
        for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
          for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
            elem<int,MAX_DIMS> Point = {i,j,k};
            elem<int,MAX_DIMS> AdjustedPoint = Cart_.PeriodicAdjust(Point);
            if (LocalRange_.Contains(AdjustedPoint)) {
              LocalDonorPointIndices_(iLocalDonorPoint) = DonorIndexer.ToIndex(Point);
              LocalDonorPointGridValuesIndices_(iLocalDonorPoint) = GridValuesIndexer_.ToIndex(
                AdjustedPoint);
              ++iLocalDonorPoint;
            }
          }
        }
      }
      NumLocalDonorPoints = iLocalDonorPoint;
    }
    for (int iCount = 0; iCount < Count_; ++iCount) {
      for (int iLocalDonorPoint = 0; iLocalDonorPoint < NumLocalDonorPoints; ++iLocalDonorPoint) {
        DonorPointValues(iCount,LocalDonorPointIndices_(iLocalDonorPoint)) = GridValues(iCount)(
          LocalDonorPointGridValuesIndices_(iLocalDonorPoint));
      }
    }

    // Fill in the remote data
    for (int iRemotePoint = 0; iRemotePoint < NumRemoteDonorPoints_(iDonor); ++iRemotePoint) {
      int iPointInCell = RemoteDonorPoints_(iDonor)[iRemotePoint];
      int iRecv = RemoteDonorPointCollectRecvs_(iDonor)[iRemotePoint];
      long long iBuffer = RemoteDonorPointCollectRecvBufferIndices_(iDonor)[iRemotePoint];
      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorPointValues(iCount,iPointInCell) = RemoteDonorValues(iRecv)(iCount,iBuffer);
      }
    }

  }

private:

  cart Cart_;
  range GlobalRange_;
  range LocalRange_;
  array_view<const exchange::collect_send> Sends_;
  array_view<const exchange::collect_recv> Recvs_;
  array<array<mpi_value_type,2>> SendBuffers_;
  array<array<mpi_value_type,2>> RecvBuffers_;
  array<MPI_Request> Requests_;
  array_view<const int> NumRemoteDonorPoints_;
  array_view<const long long * const> RemoteDonorPoints_;
  array_view<const int * const> RemoteDonorPointCollectRecvs_;
  array_view<const long long * const> RemoteDonorPointCollectRecvBufferIndices_;
  array<int> LocalDonorPointIndices_;
  array<long long> LocalDonorPointGridValuesIndices_;

};

template <typename T, array_layout Layout> class collect_none : public collect_base<T, Layout> {

protected:

  using parent_type = collect_base<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_none() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");

    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");

    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(true);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) = DonorValues(iCount)(iDonor) && !DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;

};

template <typename T, array_layout Layout> class collect_any : public collect_base<T, Layout> {

protected:

  using parent_type = collect_base<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_any() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");
    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(false);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) = DonorValues(iCount)(iDonor) || DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;

};

template <typename T, array_layout Layout> class collect_not_all : public collect_base<T, Layout> {

protected:

  using parent_type = collect_base<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_not_all() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");
    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(false);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) = DonorValues(iCount)(iDonor) || !DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;

};

template <typename T, array_layout Layout> class collect_all : public collect_base<T, Layout> {

protected:

  using parent_type = collect_base<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_all() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");
    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(true);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) = DonorValues(iCount)(iDonor) && DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;

};

template <typename T, array_layout Layout> class collect_interp : public collect_base<T, Layout> {

protected:

  using parent_type = collect_base<T, Layout>;

  using typename parent_type::donor_indexer;
  using parent_type::Donors_;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::value_type;

  collect_interp() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

    InterpCoefs_ = Donors_->InterpCoefs_;

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.Resize({{Count_,MaxPointsInCell_}});
    DonorPointCoefs_.Resize({MaxPointsInCell_});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Collect(array_view<array_view<const value_type>> GridValues, array_view<array_view<
    value_type>> DonorValues) {

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    int ReduceTime = core::GetProfilerTimerID(*Profiler_, "Collect::Reduce");
    core::StartProfile(*Profiler_, ReduceTime);

    for (long long iDonor = 0; iDonor < NumDonors_; ++iDonor) {

      elem<int,MAX_DIMS> DonorSize;
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_);
      int NumDonorPoints = DonorSize[0]*DonorSize[1]*DonorSize[2];

      donor_indexer Indexer(DonorSize);

      for (int k = 0; k < DonorSize[2]; ++k) {
        for (int j = 0; j < DonorSize[1]; ++j) {
          for (int i = 0; i < DonorSize[0]; ++i) {
            int iPointInCell = Indexer.ToIndex(i,j,k);
            DonorPointCoefs_(iPointInCell) =
              InterpCoefs_(0,i,iDonor) *
              InterpCoefs_(1,j,iDonor) *
              InterpCoefs_(2,k,iDonor);
          }
        }
      }

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues(iCount)(iDonor) = value_type(0);
        for (int iPointInCell = 0; iPointInCell < NumDonorPoints; ++iPointInCell) {
          DonorValues(iCount)(iDonor) += DonorPointCoefs_(iPointInCell)*DonorPointValues_(iCount,
            iPointInCell);
        }
      }

    }

    core::EndProfile(*Profiler_, ReduceTime);

  }

private:

  array_view<const double,3> InterpCoefs_;
  array<array<value_type,2>> RemoteDonorValues_;
  array<value_type,2> DonorPointValues_;
  array<double> DonorPointCoefs_;

};

template <array_layout Layout> void CollectLayout(const exchange &Exchange, data_type ValueType,
  int Count, collect_op CollectOp, const range &GridValuesRange, const void * const *GridValues,
  void **DonorValues) {

  collect CollectWrapper;

  switch (CollectOp) {
  case collect_op::NONE:
    switch (ValueType) {
    case data_type::BOOL: CollectWrapper = collect_none<bool, Layout>(); break;
    case data_type::BYTE: CollectWrapper = collect_none<unsigned char, Layout>(); break;
    case data_type::INT: CollectWrapper = collect_none<int, Layout>(); break;
    case data_type::LONG: CollectWrapper = collect_none<long, Layout>(); break;
    case data_type::LONG_LONG: CollectWrapper = collect_none<long long, Layout>(); break;
    case data_type::UNSIGNED_INT: CollectWrapper = collect_none<unsigned int, Layout>(); break;
    case data_type::UNSIGNED_LONG: CollectWrapper = collect_none<unsigned long, Layout>(); break;
    case data_type::UNSIGNED_LONG_LONG: CollectWrapper = collect_none<unsigned long long, Layout>(); break;
    case data_type::FLOAT: CollectWrapper = collect_none<float, Layout>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_none<double, Layout>(); break;
    }
    break;
  case collect_op::ANY:
    switch (ValueType) {
    case data_type::BOOL: CollectWrapper = collect_any<bool, Layout>(); break;
    case data_type::BYTE: CollectWrapper = collect_any<unsigned char, Layout>(); break;
    case data_type::INT: CollectWrapper = collect_any<int, Layout>(); break;
    case data_type::LONG: CollectWrapper = collect_any<long, Layout>(); break;
    case data_type::LONG_LONG: CollectWrapper = collect_any<long long, Layout>(); break;
    case data_type::UNSIGNED_INT: CollectWrapper = collect_any<unsigned int, Layout>(); break;
    case data_type::UNSIGNED_LONG: CollectWrapper = collect_any<unsigned long, Layout>(); break;
    case data_type::UNSIGNED_LONG_LONG: CollectWrapper = collect_any<unsigned long long, Layout>(); break;
    case data_type::FLOAT: CollectWrapper = collect_any<float, Layout>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_any<double, Layout>(); break;
    }
    break;
  case collect_op::NOT_ALL:
    switch (ValueType) {
    case data_type::BOOL: CollectWrapper = collect_not_all<bool, Layout>(); break;
    case data_type::BYTE: CollectWrapper = collect_not_all<unsigned char, Layout>(); break;
    case data_type::INT: CollectWrapper = collect_not_all<int, Layout>(); break;
    case data_type::LONG: CollectWrapper = collect_not_all<long, Layout>(); break;
    case data_type::LONG_LONG: CollectWrapper = collect_not_all<long long, Layout>(); break;
    case data_type::UNSIGNED_INT: CollectWrapper = collect_not_all<unsigned int, Layout>(); break;
    case data_type::UNSIGNED_LONG: CollectWrapper = collect_not_all<unsigned long, Layout>(); break;
    case data_type::UNSIGNED_LONG_LONG: CollectWrapper = collect_not_all<unsigned long long, Layout>(); break;
    case data_type::FLOAT: CollectWrapper = collect_not_all<float, Layout>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_not_all<double, Layout>(); break;
    }
    break;
  case collect_op::ALL:
    switch (ValueType) {
    case data_type::BOOL: CollectWrapper = collect_all<bool, Layout>(); break;
    case data_type::BYTE: CollectWrapper = collect_all<unsigned char, Layout>(); break;
    case data_type::INT: CollectWrapper = collect_all<int, Layout>(); break;
    case data_type::LONG: CollectWrapper = collect_all<long, Layout>(); break;
    case data_type::LONG_LONG: CollectWrapper = collect_all<long long, Layout>(); break;
    case data_type::UNSIGNED_INT: CollectWrapper = collect_all<unsigned int, Layout>(); break;
    case data_type::UNSIGNED_LONG: CollectWrapper = collect_all<unsigned long, Layout>(); break;
    case data_type::UNSIGNED_LONG_LONG: CollectWrapper = collect_all<unsigned long long, Layout>(); break;
    case data_type::FLOAT: CollectWrapper = collect_all<float, Layout>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_all<double, Layout>(); break;
    }
    break;
  case collect_op::INTERPOLATE:
    switch (ValueType) {
    case data_type::FLOAT: CollectWrapper = collect_interp<float, Layout>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_interp<double, Layout>(); break;
    default:
      OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
      break;
    }
    break;
  }

  CollectWrapper.Initialize(Exchange, Count, GridValuesRange);
  CollectWrapper.Collect(GridValues, DonorValues);

}

}

namespace core {

void Collect(const exchange &Exchange, data_type ValueType, int Count, collect_op CollectOp,
  const range &GridValuesRange, array_layout GridValuesLayout, const void * const *GridValues,
  void **DonorValues) {

  switch (GridValuesLayout) {
  case array_layout::ROW_MAJOR:
    CollectLayout<array_layout::ROW_MAJOR>(Exchange, ValueType, Count, CollectOp, GridValuesRange,
      GridValues, DonorValues); break;
  case array_layout::COLUMN_MAJOR:
    CollectLayout<array_layout::COLUMN_MAJOR>(Exchange, ValueType, Count, CollectOp,
      GridValuesRange, GridValues, DonorValues); break;
  }

}

}

namespace {

class send {

public:

  send() = default;
  send(const send &) = delete;
  send(send &&) = default;

  template <typename T> explicit send(T &&Send):
    Send_(new model<T>(std::forward<T>(Send)))
  {}

  template <typename T> send &operator=(T &&Send) {
    Send_.reset(new model<T>(std::forward<T>(Send)));
    return *this;
  }

  void Initialize(const exchange &Exchange, int Count, int Tag) {
    Send_->Initialize(Exchange, Count, Tag);
  }

  void Send(const void * const *DonorValues, request &Request) {
    Send_->Send(DonorValues, Request);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Initialize(const exchange &Exchange, int Count, int Tag) = 0;
    virtual void Send(const void * const *DonorValues, request &Request) = 0;
  };

  template <typename T> class model : public concept {
  public:
    using value_type = typename T::value_type;
    explicit model(T Send):
      Send_(std::move(Send))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, int Tag) override {
      Send_.Initialize(Exchange, Count, Tag);
      const connectivity_d *Donors;
      GetConnectivityDonorSide(*Exchange.Connectivity_, Donors);
      GetConnectivityDonorSideCount(*Donors, NumDonors_);
      DonorValues_.Resize({Count});
    }
    virtual void Send(const void * const *DonorValuesVoid, request &Request) override {
      OVK_DEBUG_ASSERT(DonorValuesVoid || DonorValues_.Count() == 0, "Invalid donor values "
        "pointer.");
      for (int iCount = 0; iCount < DonorValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(DonorValuesVoid[iCount] || NumDonors_ == 0, "Invalid donor values "
          "pointer.");
        DonorValues_(iCount) = {static_cast<const value_type *>(DonorValuesVoid[iCount]),
          {NumDonors_}};
      }
      Send_.Send(DonorValues_, Request);
    }
  private:
    T Send_;
    long long NumDonors_;
    array<array_view<const value_type>> DonorValues_;
  };

  std::unique_ptr<concept> Send_;

};

template <typename T> class send_request {
public:
  using value_type = T;
  using mpi_value_type = core::mpi_compatible_type<value_type>;
  send_request(const exchange &Exchange, int Count, int NumSends, array<array<mpi_value_type,2>>
    Buffers, array<MPI_Request> MPIRequests):
    Exchange_(&Exchange),
    Count_(Count),
    NumSends_(NumSends),
    Buffers_(std::move(Buffers)),
    MPIRequests_(std::move(MPIRequests)),
    Profiler_(Exchange.Profiler_),
    MemAllocTime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc")),
    MPITime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::MPI"))
  {
    GetConnectivityDonorSide(*Exchange.Connectivity_, Donors_);
  }
  array_view<MPI_Request> MPIRequests() { return MPIRequests_; }
  void Finish(int) { /* Nothing to finish */ }
  void Wait();
  void StartProfileMemAlloc() const { core::StartProfile(*Profiler_, MemAllocTime_); }
  void EndProfileMemAlloc() const { core::EndProfile(*Profiler_, MemAllocTime_); }
  void StartProfileMPI() const { core::StartProfile(*Profiler_, MPITime_); }
  void EndProfileMPI() const { core::EndProfile(*Profiler_, MPITime_); }
private:
  const exchange *Exchange_;
  const connectivity_d *Donors_;
  int Count_;
  int NumSends_;
  array<array<mpi_value_type,2>> Buffers_;
  array<MPI_Request> MPIRequests_;
  mutable core::profiler *Profiler_;
  int MemAllocTime_;
  int MPITime_;
};

template <typename T> class send_impl {

public:

  using value_type = T;
  using mpi_value_type = core::mpi_compatible_type<value_type>;
  using request_type = send_request<T>;

  send_impl() = default;
  send_impl(const send_impl &) = delete;
  send_impl(send_impl &&) = default;

  void Initialize(const exchange &Exchange, int Count, int Tag) {

    Exchange_ = &Exchange;
    Count_ = Count;
    Tag_ = Tag;

    const connectivity &Connectivity = *Exchange.Connectivity_;

    GetConnectivityDonorSide(Connectivity, Donors_);

    Profiler_ = Exchange.Profiler_;

    NumSends_ = Exchange.Sends_.Count();

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    NextBufferEntry_.Resize({NumSends_});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  void Send(array_view<array_view<const value_type>> DonorValues, request &Request) {

    const exchange &Exchange = *Exchange_;
    const connectivity_d &Donors = *Donors_;

    long long NumDonors;
    GetConnectivityDonorSideCount(Donors, NumDonors);

    MPI_Datatype MPIDataType = core::GetMPIDataType<mpi_value_type>();

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc");
    int PackTime = core::GetProfilerTimerID(*Profiler_, "SendRecv::Pack");
    int MPITime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MPI");

    core::StartProfile(*Profiler_, MemAllocTime);

    // Will be moved into request object
    array<array<mpi_value_type,2>> Buffers({NumSends_});
    for (int iSend = 0; iSend < NumSends_; ++iSend) {
      const exchange::send &Send = Exchange.Sends_(iSend);
      Buffers(iSend).Resize({{Count_,Send.Count}});
    }

    core::EndProfile(*Profiler_, MemAllocTime);
    core::StartProfile(*Profiler_, PackTime);

    for (auto &iBuffer : NextBufferEntry_) {
      iBuffer = 0;
    }

    for (long long iDonorOrder = 0; iDonorOrder < NumDonors; ++iDonorOrder) {
      long long iDonor = Exchange.DonorsSorted_(iDonorOrder);
      int iSend = Exchange.DonorSendIndices_(iDonor);
      if (iSend >= 0) {
        long long iBuffer = NextBufferEntry_(iSend);
        for (int iCount = 0; iCount < Count_; ++iCount) {
          Buffers(iSend)(iCount,iBuffer) = mpi_value_type(DonorValues(iCount)(iDonor));
        }
        ++NextBufferEntry_(iSend);
      }
    }

    core::EndProfile(*Profiler_, PackTime);
    core::StartProfile(*Profiler_, MemAllocTime);

    // Will be moved into request object
    array<MPI_Request> MPIRequests({NumSends_});

    core::EndProfile(*Profiler_, MemAllocTime);
    core::StartProfile(*Profiler_, MPITime);

    for (int iSend = 0; iSend < NumSends_; ++iSend) {
      const exchange::send &Send = Exchange.Sends_(iSend);
      MPI_Isend(Buffers(iSend).Data(), Count_*Send.Count, MPIDataType, Send.Rank, Tag_,
        Exchange.Comm_, &MPIRequests(iSend));
    }

    core::EndProfile(*Profiler_, MPITime);

    Request = request_type(Exchange, Count_, NumSends_, std::move(Buffers), std::move(MPIRequests));

  }

private:

  const exchange *Exchange_;
  const connectivity_d *Donors_;
  core::profiler *Profiler_;
  int Count_;
  int Tag_;
  int NumSends_;
  array<long long> NextBufferEntry_;

};

template <typename T> void send_request<T>::Wait() {

  core::StartProfile(*Profiler_, MPITime_);

  MPI_Waitall(NumSends_, MPIRequests_.Data(), MPI_STATUSES_IGNORE);

  core::EndProfile(*Profiler_, MPITime_);

  Buffers_.Clear();
  MPIRequests_.Clear();

}

}

namespace core {

void Send(const exchange &Exchange, data_type ValueType, int Count, const void * const *DonorValues,
  int Tag, request &Request) {

  send SendWrapper;

  switch (ValueType) {
  case data_type::BOOL: SendWrapper = send_impl<bool>(); break;
  case data_type::BYTE: SendWrapper = send_impl<unsigned char>(); break;
  case data_type::INT: SendWrapper = send_impl<int>(); break;
  case data_type::LONG: SendWrapper = send_impl<long>(); break;
  case data_type::LONG_LONG: SendWrapper = send_impl<long long>(); break;
  case data_type::UNSIGNED_INT: SendWrapper = send_impl<unsigned int>(); break;
  case data_type::UNSIGNED_LONG: SendWrapper = send_impl<unsigned long>(); break;
  case data_type::UNSIGNED_LONG_LONG: SendWrapper = send_impl<unsigned long long>(); break;
  case data_type::FLOAT: SendWrapper = send_impl<float>(); break;
  case data_type::DOUBLE: SendWrapper = send_impl<double>(); break;
  }

  SendWrapper.Initialize(Exchange, Count, Tag);
  SendWrapper.Send(DonorValues, Request);

}

}

namespace {

class recv {

public:

  recv() = default;
  recv(const recv &) = delete;
  recv(recv &&) = default;

  template <typename T> explicit recv(T &&Recv):
    Recv_(new model<T>(std::forward<T>(Recv)))
  {}

  template <typename T> recv &operator=(T &&Recv) {
    Recv_.reset(new model<T>(std::forward<T>(Recv)));
    return *this;
  }

  void Initialize(const exchange &Exchange, int Count, int Tag) {
    Recv_->Initialize(Exchange, Count, Tag);
  }

  void Recv(void **ReceiverValues, request &Request) {
    Recv_->Recv(ReceiverValues, Request);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Initialize(const exchange &Exchange, int Count, int Tag) = 0;
    virtual void Recv(void **ReceiverValues, request &Request) = 0;
  };

  template <typename T> class model : public concept {
  public:
    using value_type = typename T::value_type;
    explicit model(T Recv):
      Recv_(std::move(Recv))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, int Tag) override {
      Recv_.Initialize(Exchange, Count, Tag);
      const connectivity_r *Receivers;
      GetConnectivityReceiverSide(*Exchange.Connectivity_, Receivers);
      GetConnectivityReceiverSideCount(*Receivers, NumReceivers_);
      ReceiverValues_.Resize({Count});
    }
    virtual void Recv(void **ReceiverValuesVoid, request &Request) override {
      OVK_DEBUG_ASSERT(ReceiverValuesVoid || ReceiverValues_.Count() == 0, "Invalid receiver values "
        "pointer.");
      for (int iCount = 0; iCount < ReceiverValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(ReceiverValuesVoid[iCount] || NumReceivers_ == 0, "Invalid receiver data "
          "pointer.");
        ReceiverValues_(iCount) = {static_cast<value_type *>(ReceiverValuesVoid[iCount]),
          {NumReceivers_}};
      }
      Recv_.Recv(ReceiverValues_, Request);
    }
  private:
    T Recv_;
    long long NumReceivers_;
    array<array_view<value_type>> ReceiverValues_;
  };

  std::unique_ptr<concept> Recv_;

};

template <typename T> class recv_request {
public:
  using value_type = T;
  using mpi_value_type = core::mpi_compatible_type<value_type>;
  recv_request(const exchange &Exchange, int Count, int NumRecvs, array<array<mpi_value_type,2>>
    Buffers, array<MPI_Request> MPIRequests, array<array_view<value_type>> ReceiverValues):
    Exchange_(&Exchange),
    Count_(Count),
    NumRecvs_(NumRecvs),
    Buffers_(std::move(Buffers)),
    MPIRequests_(std::move(MPIRequests)),
    ReceiverValues_(std::move(ReceiverValues)),
    Profiler_(Exchange.Profiler_),
    MemAllocTime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc")),
    MPITime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::MPI")),
    UnpackTime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::Unpack"))
  {
    GetConnectivityReceiverSide(*Exchange.Connectivity_, Receivers_);
  }
  array_view<MPI_Request> MPIRequests() { return MPIRequests_; }
  void Finish(int) { /* Can't finish until all requests are done */ }
  void Wait();
  void StartProfileMemAlloc() const { core::StartProfile(*Profiler_, MemAllocTime_); }
  void EndProfileMemAlloc() const { core::EndProfile(*Profiler_, MemAllocTime_); }
  void StartProfileMPI() const { core::StartProfile(*Profiler_, MPITime_); }
  void EndProfileMPI() const { core::EndProfile(*Profiler_, MPITime_); }
private:
  const exchange *Exchange_;
  const connectivity_r *Receivers_;
  int Count_;
  int NumRecvs_;
  array<array<mpi_value_type,2>> Buffers_;
  array<MPI_Request> MPIRequests_;
  array<array_view<value_type>> ReceiverValues_;
  mutable core::profiler *Profiler_;
  int MemAllocTime_;
  int MPITime_;
  int UnpackTime_;
};

template <typename T> class recv_impl {

public:

  using value_type = T;
  using mpi_value_type = core::mpi_compatible_type<value_type>;
  using request_type = recv_request<T>;

  recv_impl() = default;
  recv_impl(const recv_impl &) = delete;
  recv_impl(recv_impl &&) = default;

  void Initialize(const exchange &Exchange, int Count, int Tag) {

    Exchange_ = &Exchange;
    Count_ = Count;
    Tag_ = Tag;

    const connectivity &Connectivity = *Exchange.Connectivity_;

    GetConnectivityReceiverSide(Connectivity, Receivers_);

    Profiler_ = Exchange.Profiler_;

    NumRecvs_ = Exchange.Recvs_.Count();

  }

  void Recv(array_view<array_view<value_type>> ReceiverValues, request &Request) {

    const exchange &Exchange = *Exchange_;

    MPI_Datatype MPIDataType = core::GetMPIDataType<mpi_value_type>();

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc");
    int MPITime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MPI");

    core::StartProfile(*Profiler_, MemAllocTime);

    // Will be moved into request object
    array<array<mpi_value_type,2>> Buffers({NumRecvs_});
    for (int iRecv = 0; iRecv < NumRecvs_; ++iRecv) {
      const exchange::recv &Recv = Exchange.Recvs_(iRecv);
      Buffers(iRecv).Resize({{Count_,Recv.Count}});
    }

    // Will be moved into request object
    array<MPI_Request> MPIRequests({NumRecvs_});

    core::EndProfile(*Profiler_, MemAllocTime);
    core::StartProfile(*Profiler_, MPITime);

    for (int iRecv = 0; iRecv < NumRecvs_; ++iRecv) {
      const exchange::recv &Recv = Exchange.Recvs_(iRecv);
      MPI_Irecv(Buffers(iRecv).Data(), Count_*Recv.Count, MPIDataType, Recv.Rank, Tag_,
        Exchange.Comm_, &MPIRequests(iRecv));
    }

    core::EndProfile(*Profiler_, MPITime);
    core::StartProfile(*Profiler_, MemAllocTime);

    // Will be moved into request object
    array<array_view<value_type>> ReceiverValuesSaved({ReceiverValues});

    core::EndProfile(*Profiler_, MemAllocTime);

    Request = request_type(Exchange, Count_, NumRecvs_, std::move(Buffers), std::move(MPIRequests),
      std::move(ReceiverValuesSaved));

  }

private:

  const exchange *Exchange_;
  const connectivity_r *Receivers_;
  core::profiler *Profiler_;
  int Count_;
  int Tag_;
  int NumRecvs_;

};

template <typename T> void recv_request<T>::Wait() {

  const exchange &Exchange = *Exchange_;

  const connectivity_r &Receivers = *Receivers_;

  long long NumReceivers;
  GetConnectivityReceiverSideCount(Receivers, NumReceivers);

  core::StartProfile(*Profiler_, MemAllocTime_);

  array<long long> NextBufferEntry({NumRecvs_}, 0);

  core::EndProfile(*Profiler_, MemAllocTime_);
  core::StartProfile(*Profiler_, MPITime_);

  MPI_Waitall(NumRecvs_, MPIRequests_.Data(), MPI_STATUSES_IGNORE);

  core::EndProfile(*Profiler_, MPITime_);
  core::StartProfile(*Profiler_, UnpackTime_);

  for (long long iReceiverOrder = 0; iReceiverOrder < NumReceivers; ++iReceiverOrder) {
    long long iReceiver = Exchange.ReceiversSorted_(iReceiverOrder);
    int iRecv = Exchange.ReceiverRecvIndices_(iReceiver);
    if (iRecv >= 0) {
      long long iBuffer = NextBufferEntry(iRecv);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        ReceiverValues_(iCount)(iReceiver) = value_type(Buffers_(iRecv)(iCount,iBuffer));
      }
      ++NextBufferEntry(iRecv);
    }
  }

  core::EndProfile(*Profiler_, UnpackTime_);

  Buffers_.Clear();
  MPIRequests_.Clear();
  ReceiverValues_.Clear();

}

}

namespace core {

void Receive(const exchange &Exchange, data_type ValueType, int Count, void **ReceiverValues,
  int Tag, request &Request) {

  recv RecvWrapper;

  switch (ValueType) {
  case data_type::BOOL: RecvWrapper = recv_impl<bool>(); break;
  case data_type::BYTE: RecvWrapper = recv_impl<unsigned char>(); break;
  case data_type::INT: RecvWrapper = recv_impl<int>(); break;
  case data_type::LONG: RecvWrapper = recv_impl<long>(); break;
  case data_type::LONG_LONG: RecvWrapper = recv_impl<long long>(); break;
  case data_type::UNSIGNED_INT: RecvWrapper = recv_impl<unsigned int>(); break;
  case data_type::UNSIGNED_LONG: RecvWrapper = recv_impl<unsigned long>(); break;
  case data_type::UNSIGNED_LONG_LONG: RecvWrapper = recv_impl<unsigned long long>(); break;
  case data_type::FLOAT: RecvWrapper = recv_impl<float>(); break;
  case data_type::DOUBLE: RecvWrapper = recv_impl<double>(); break;
  }

  RecvWrapper.Initialize(Exchange, Count, Tag);
  RecvWrapper.Recv(ReceiverValues, Request);

}

}

namespace {

class disperse {

public:

  disperse() = default;
  disperse(const disperse &) = delete;
  disperse(disperse &&) = default;

  template <typename T> explicit disperse(T &&Disperse):
    Disperse_(new model<T>(std::forward<T>(Disperse)))
  {}

  template <typename T> disperse &operator=(T &&Disperse) {
    Disperse_.reset(new model<T>(std::forward<T>(Disperse)));
    return *this;
  }

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {
    Disperse_->Initialize(Exchange, Count, GridValuesRange);
  }

  void Disperse(const void * const *ReceiverValues, void **GridValues) {
    Disperse_->Disperse(ReceiverValues, GridValues);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) = 0;
    virtual void Disperse(const void * const *ReceiverValues, void **GridValues) = 0;
  };

  template <typename T> class model : public concept {
  public:
    using value_type = typename T::value_type;
    explicit model(T Disperse):
      Disperse_(std::move(Disperse))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange)
      override {
      Disperse_.Initialize(Exchange, Count, GridValuesRange);
      const connectivity_r *Receivers;
      GetConnectivityReceiverSide(*Exchange.Connectivity_, Receivers);
      GetConnectivityReceiverSideCount(*Receivers, NumReceivers_);
      GridValuesRange_ = GridValuesRange;
      ReceiverValues_.Resize({Count});
      GridValues_.Resize({Count});
    }
    virtual void Disperse(const void * const *ReceiverValuesVoid, void **GridValuesVoid) override {
      OVK_DEBUG_ASSERT(ReceiverValuesVoid || ReceiverValues_.Count() == 0, "Invalid receiver values "
        "pointer.");
      OVK_DEBUG_ASSERT(GridValuesVoid || GridValues_.Count() == 0, "Invalid grid values pointer.");
      for (int iCount = 0; iCount < ReceiverValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(ReceiverValuesVoid[iCount] || NumReceivers_ == 0, "Invalid receiver "
          "values pointer.");
        ReceiverValues_(iCount) = {static_cast<const value_type *>(ReceiverValuesVoid[iCount]),
          {NumReceivers_}};
      }
      for (int iCount = 0; iCount < GridValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(GridValuesVoid[iCount] || NumReceivers_ == 0, "Invalid grid values "
          "pointer.");
        GridValues_(iCount) = {static_cast<value_type *>(GridValuesVoid[iCount]),
          {GridValuesRange_.Count()}};
      }
      Disperse_.Disperse(ReceiverValues_, GridValues_);
    }
  private:
    T Disperse_;
    long long NumReceivers_;
    range GridValuesRange_;
    array<array_view<const value_type>> ReceiverValues_;
    array<array_view<value_type>> GridValues_;
  };

  std::unique_ptr<concept> Disperse_;

};

template <typename T, array_layout Layout> class disperse_base {

public:

  using value_type = T;

  disperse_base() = default;
  disperse_base(const disperse_base &) = delete;
  disperse_base(disperse_base &&) = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    Count_ = Count;
    GridValuesRange_ = GridValuesRange;
    GridValuesIndexer_ = range_indexer(GridValuesRange);

    const connectivity &Connectivity = *Exchange.Connectivity_;

    GetConnectivityReceiverSide(Connectivity, Receivers_);
    const connectivity_r &Receivers = *Receivers_;

    GetConnectivityReceiverSideCount(Receivers, NumReceivers_);

    GetConnectivityReceiverSideGrid(Receivers, Grid_);
    const grid &Grid = *Grid_;

    if (OVK_DEBUG) {
      range LocalRange;
      GetGridLocalRange(Grid, LocalRange);
      OVK_DEBUG_ASSERT(GridValuesRange.Includes(LocalRange), "Invalid grid values range.");
    }

    Points_ = Receivers.Points_;

  }

protected:

  using range_indexer = indexer<long long, int, MAX_DIMS, Layout>;

  const connectivity_r *Receivers_;
  const grid *Grid_;
  int Count_;
  long long NumReceivers_;
  range GridValuesRange_;
  range_indexer GridValuesIndexer_;
  array_view<const int,2> Points_;

};

template <typename T, array_layout Layout> class disperse_overwrite : public disperse_base<T,
  Layout> {

protected:

  using parent_type = disperse_base<T, Layout>;

  using parent_type::Receivers_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::NumReceivers_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesIndexer_;
  using parent_type::Points_;

public:

  using typename parent_type::value_type;

  disperse_overwrite() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

  }

  void Disperse(array_view<array_view<const value_type>> ReceiverValues, array_view<
    array_view<value_type>> GridValues) {

    for (long long iReceiver = 0; iReceiver < NumReceivers_; ++iReceiver) {
      elem<int,MAX_DIMS> Point = {
        Points_(0,iReceiver),
        Points_(1,iReceiver),
        Points_(2,iReceiver)
      };
      long long iPoint = GridValuesIndexer_.ToIndex(Point);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        GridValues(iCount)(iPoint) = ReceiverValues(iCount)(iReceiver);
      }
    }

  }

};

template <typename T, array_layout Layout> class disperse_append : public disperse_base<T, Layout> {

protected:

  using parent_type = disperse_base<T, Layout>;

  using parent_type::Receivers_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::NumReceivers_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesIndexer_;
  using parent_type::Points_;

public:

  using typename parent_type::value_type;

  disperse_append() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

  }

  void Disperse(array_view<array_view<const value_type>> ReceiverValues, array_view<
    array_view<value_type>> GridValues) {

    for (long long iReceiver = 0; iReceiver < NumReceivers_; ++iReceiver) {
      elem<int,MAX_DIMS> Point = {
        Points_(0,iReceiver),
        Points_(1,iReceiver),
        Points_(2,iReceiver)
      };
      long long iPoint = GridValuesIndexer_.ToIndex(Point);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        GridValues(iCount)(iPoint) += ReceiverValues(iCount)(iReceiver);
      }
    }

  }

};

template <array_layout Layout> void DisperseLayout(const exchange &Exchange, data_type ValueType,
  int Count, disperse_op DisperseOp, const void * const *ReceiverValues, const range
  &GridValuesRange, void **GridValues) {

  disperse DisperseWrapper;

  switch (DisperseOp) {
  case disperse_op::OVERWRITE:
    switch (ValueType) {
    case data_type::BOOL: DisperseWrapper = disperse_overwrite<bool, Layout>(); break;
    case data_type::BYTE: DisperseWrapper = disperse_overwrite<unsigned char, Layout>(); break;
    case data_type::INT: DisperseWrapper = disperse_overwrite<int, Layout>(); break;
    case data_type::LONG: DisperseWrapper = disperse_overwrite<long, Layout>(); break;
    case data_type::LONG_LONG: DisperseWrapper = disperse_overwrite<long long, Layout>(); break;
    case data_type::UNSIGNED_INT: DisperseWrapper = disperse_overwrite<unsigned int, Layout>(); break;
    case data_type::UNSIGNED_LONG: DisperseWrapper = disperse_overwrite<unsigned long, Layout>(); break;
    case data_type::UNSIGNED_LONG_LONG: DisperseWrapper = disperse_overwrite<unsigned long long, Layout>(); break;
    case data_type::FLOAT: DisperseWrapper = disperse_overwrite<float, Layout>(); break;
    case data_type::DOUBLE: DisperseWrapper = disperse_overwrite<double, Layout>(); break;
    }
    break;
  case disperse_op::APPEND:
    switch (ValueType) {
    case data_type::BYTE: DisperseWrapper = disperse_append<unsigned char, Layout>(); break;
    case data_type::INT: DisperseWrapper = disperse_append<int, Layout>(); break;
    case data_type::LONG: DisperseWrapper = disperse_append<long, Layout>(); break;
    case data_type::LONG_LONG: DisperseWrapper = disperse_append<long long, Layout>(); break;
    case data_type::UNSIGNED_INT: DisperseWrapper = disperse_append<unsigned int, Layout>(); break;
    case data_type::UNSIGNED_LONG: DisperseWrapper = disperse_append<unsigned long, Layout>(); break;
    case data_type::UNSIGNED_LONG_LONG: DisperseWrapper = disperse_append<unsigned long long, Layout>(); break;
    case data_type::FLOAT: DisperseWrapper = disperse_append<float, Layout>(); break;
    case data_type::DOUBLE: DisperseWrapper = disperse_append<double, Layout>(); break;
    default:
      OVK_DEBUG_ASSERT(false, "Invalid data type for append disperse operation.");
      break;
    }
    break;
  }

  DisperseWrapper.Initialize(Exchange, Count, GridValuesRange);
  DisperseWrapper.Disperse(ReceiverValues, GridValues);

}

}

namespace core {

void Disperse(const exchange &Exchange, data_type ValueType, int Count, disperse_op DisperseOp,
  const void * const *ReceiverValues, const range &GridValuesRange, array_layout GridValuesLayout,
  void **GridValues) {

  switch (GridValuesLayout) {
  case array_layout::ROW_MAJOR:
    DisperseLayout<array_layout::ROW_MAJOR>(Exchange, ValueType, Count, DisperseOp,
      ReceiverValues, GridValuesRange, GridValues); break;
  case array_layout::COLUMN_MAJOR:
    DisperseLayout<array_layout::COLUMN_MAJOR>(Exchange, ValueType, Count, DisperseOp,
      ReceiverValues, GridValuesRange, GridValues); break;
  }

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
