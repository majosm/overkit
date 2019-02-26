// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Exchange.hpp"

#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

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
  error_handler &ErrorHandler) {

  Exchange.Connectivity_ = &Connectivity;

  Exchange.Logger_ = &Logger;
  Exchange.ErrorHandler_ = &ErrorHandler;

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

  range DonorGridLocalRange(Exchange.NumDims_);
  if (DonorGrid) {
    GetGridLocalRange(*DonorGrid, DonorGridLocalRange);
  }

  CreatePartitionHash(Exchange.SourceHash_, Exchange.NumDims_, Exchange.Comm_, DonorGridGlobalRange,
    DonorGridLocalRange);

  range ReceiverGridGlobalRange;
  GetGridInfoGlobalRange(*ReceiverGridInfo, ReceiverGridGlobalRange);

  range ReceiverGridLocalRange(Exchange.NumDims_);
  if (ReceiverGrid) {
    GetGridLocalRange(*ReceiverGrid, ReceiverGridLocalRange);
  }

  CreatePartitionHash(Exchange.DestinationHash_, Exchange.NumDims_, Exchange.Comm_,
    ReceiverGridGlobalRange, ReceiverGridLocalRange);

  Exchange.RemoteDonorPoints_ = nullptr;
  Exchange.RemoteDonorPointCollectRecvs_ = nullptr;
  Exchange.RemoteDonorPointCollectRecvBufferIndices_ = nullptr;

  MPI_Barrier(Exchange.Comm_);

  core::LogStatus(*Exchange.Logger_, Exchange.Comm_.Rank() == 0, 0, "Created exchange %s.",
    Connectivity.Name_);

}

void DestroyExchange(exchange &Exchange) {

  MPI_Barrier(Exchange.Comm_);

  const connectivity &Connectivity = *Exchange.Connectivity_;

  Exchange.Sends_.clear();
  Exchange.Recvs_.clear();

  Exchange.NumRemoteDonorPoints_.clear();

  Exchange.RemoteDonorPointsPtrs_.clear();
  Exchange.RemoteDonorPointsData_.clear();

  Exchange.RemoteDonorPointCollectRecvsPtrs_.clear();
  Exchange.RemoteDonorPointCollectRecvsData_.clear();

  Exchange.RemoteDonorPointCollectRecvBufferIndicesPtrs_.clear();
  Exchange.RemoteDonorPointCollectRecvBufferIndicesData_.clear();

  Exchange.CollectSends_.clear();
  Exchange.CollectRecvs_.clear();

  DestroyPartitionHash(Exchange.SourceHash_);
  DestroyPartitionHash(Exchange.DestinationHash_);

  Exchange.DonorsSorted_.clear();
  Exchange.ReceiversSorted_.clear();

  Exchange.DonorDestRanks_.clear();
  Exchange.ReceiverSourceRanks_.clear();

  MPI_Barrier(Exchange.Comm_);

  LogStatus(*Exchange.Logger_, Exchange.Comm_.Rank() == 0, 0, "Destroyed exchange %s.",
    Connectivity.Name_);

  Exchange.Comm_.Reset();

}

void CreateExchangeInfo(exchange_info &Info, const exchange *Exchange, const comm &Comm) {

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
  std::vector<char> NameChars(NameLength);
  if (IsRoot) NameChars.assign(Connectivity->Name_.begin(), Connectivity->Name_.end());
  MPI_Bcast(NameChars.data(), NameLength, MPI_CHAR, RootRank, Comm);
  Info.Name_.assign(NameChars.begin(), NameChars.end());

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

  Exchange.CollectSends_.clear();

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

    const std::vector<core::grid_neighbor> &GridNeighbors = core::GetGridNeighbors(*Grid);
    int NumNeighbors = GridNeighbors.size();

    cart Cart;
    GetGridCart(*Grid, Cart);

    std::vector<range> SendToNeighborRanges(NumNeighbors);
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      SendToNeighborRanges[iNeighbor] = range(NumDims);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS], DonorEnd[MAX_DIMS];
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorBegin[iDim] = Donors->Extents_[0][iDim][iDonor];
        DonorEnd[iDim] = Donors->Extents_[1][iDim][iDonor];
      }
      range DonorRange(NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = RangeIncludes(GlobalRange, DonorRange);
      for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
        if (AwayFromEdge) {
          bool Overlaps = RangesOverlap(GridNeighbors[iNeighbor].LocalRange, DonorRange);
          if (Overlaps) {
            SendToNeighborRanges[iNeighbor] = UnionRanges(SendToNeighborRanges[iNeighbor],
              IntersectRanges(LocalRange, DonorRange));
          }
        } else {
          bool Overlaps = false;
          for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
            for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
              for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                CartPeriodicAdjust(Cart, Point, Point);
                if (RangeContains(GridNeighbors[iNeighbor].LocalRange, Point)) {
                  Overlaps = true;
                  goto done_checking_for_overlap1;
                }
              }
            }
          }
          done_checking_for_overlap1:;
          if (Overlaps) {
            for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
              for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
                for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
                  int Point[MAX_DIMS] = {i, j, k};
                  CartPeriodicAdjust(Cart, Point, Point);
                  if (RangeContains(LocalRange, Point)) {
                    ExtendRange(SendToNeighborRanges[iNeighbor], Point);
                  }
                }
              }
            }
          }
        }
      }
    }

    std::vector<int> CollectSendIndexToNeighbor;
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!SendToNeighborRanges[iNeighbor].Empty()) {
        CollectSendIndexToNeighbor.push_back(iNeighbor);
      }
    }
    int NumCollectSends = CollectSendIndexToNeighbor.size();

    Exchange.CollectSends_.resize(NumCollectSends);

    for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      int iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
      Exchange.CollectSends_[iCollectSend].Rank = GridNeighbors[iNeighbor].Rank;
    }

    std::vector<std::vector<char>> CollectSendMasks(NumCollectSends);
    for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      int iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
      CollectSendMasks[iCollectSend].resize(SendToNeighborRanges[iNeighbor].Count(), false);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS], DonorEnd[MAX_DIMS];
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorBegin[iDim] = Donors->Extents_[0][iDim][iDonor];
        DonorEnd[iDim] = Donors->Extents_[1][iDim][iDonor];
      }
      range DonorRange(NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = RangeIncludes(GlobalRange, DonorRange);
      for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
        int iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
        if (AwayFromEdge) {
          bool Overlaps = RangesOverlap(GridNeighbors[iNeighbor].LocalRange, DonorRange);
          if (Overlaps) {
            range LocalDonorRange = IntersectRanges(LocalRange, DonorRange);
            for (int k = LocalDonorRange.Begin(2); k < LocalDonorRange.End(2); ++k) {
              for (int j = LocalDonorRange.Begin(1); j < LocalDonorRange.End(1); ++j) {
                for (int i = LocalDonorRange.Begin(0); i < LocalDonorRange.End(0); ++i) {
                  int Point[MAX_DIMS] = {i, j, k};
                  long long iPoint = RangeTupleToIndex(SendToNeighborRanges[iNeighbor],
                    array_layout::COLUMN_MAJOR, Point);
                  CollectSendMasks[iCollectSend][iPoint] = true;
                }
              }
            }
          }
        } else {
          bool Overlaps = false;
          for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
            for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
              for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                CartPeriodicAdjust(Cart, Point, Point);
                if (RangeContains(GridNeighbors[iNeighbor].LocalRange, Point)) {
                  Overlaps = true;
                  goto done_checking_for_overlap2;
                }
              }
            }
          }
          done_checking_for_overlap2:;
          if (Overlaps) {
            for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
              for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
                for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
                  int Point[MAX_DIMS] = {i, j, k};
                  CartPeriodicAdjust(Cart, Point, Point);
                  if (RangeContains(LocalRange, Point)) {
                    long long iPoint = RangeTupleToIndex(SendToNeighborRanges[iNeighbor],
                      array_layout::COLUMN_MAJOR, Point);
                    CollectSendMasks[iCollectSend][iPoint] = true;
                  }
                }
              }
            }
          }
        }
      }
    }

    for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      exchange::collect_send &CollectSend = Exchange.CollectSends_[iCollectSend];
      int iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
      CollectSend.NumPoints = 0;
      for (long long iPoint = 0; iPoint < SendToNeighborRanges[iNeighbor].Count(); ++iPoint) {
        if (CollectSendMasks[iCollectSend][iPoint]) {
          ++CollectSend.NumPoints;
        }
      }
    }

    for (int iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      exchange::collect_send &CollectSend = Exchange.CollectSends_[iCollectSend];
      int iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
      CollectSend.PointsData.resize(MAX_DIMS*CollectSend.NumPoints);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        CollectSend.Points[iDim] = CollectSend.PointsData.data() + iDim*CollectSend.NumPoints;
      }
      long long iCollectSendPoint = 0;
      for (long long iPoint = 0; iPoint < SendToNeighborRanges[iNeighbor].Count(); ++iPoint) {
        if (CollectSendMasks[iCollectSend][iPoint]) {
          std::array<int,MAX_DIMS> Point = RangeIndexToTuple(SendToNeighborRanges[iNeighbor],
            array_layout::COLUMN_MAJOR, iPoint);
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            CollectSend.Points[iDim][iCollectSendPoint] = Point[iDim];
          }
          ++iCollectSendPoint;
        }
      }
    }

  }

}

void UpdateCollectReceiveInfo(exchange &Exchange) {

  Exchange.CollectRecvs_.clear();

  Exchange.RemoteDonorPoints_ = nullptr;
  Exchange.RemoteDonorPointsPtrs_.clear();
  Exchange.RemoteDonorPointsData_.clear();
  Exchange.RemoteDonorPointCollectRecvs_ = nullptr;
  Exchange.RemoteDonorPointCollectRecvsPtrs_.clear();
  Exchange.RemoteDonorPointCollectRecvsData_.clear();
  Exchange.RemoteDonorPointCollectRecvBufferIndices_ = nullptr;
  Exchange.RemoteDonorPointCollectRecvBufferIndicesPtrs_.clear();
  Exchange.RemoteDonorPointCollectRecvBufferIndicesData_.clear();
  
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

    const std::vector<core::grid_neighbor> &GridNeighbors = core::GetGridNeighbors(*Grid);
    int NumNeighbors = GridNeighbors.size();

    cart Cart;
    GetGridCart(*Grid, Cart);

    std::vector<range> RecvFromNeighborRanges(NumNeighbors);
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      RecvFromNeighborRanges[iNeighbor] = range(NumDims);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS], DonorEnd[MAX_DIMS];
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorBegin[iDim] = Donors->Extents_[0][iDim][iDonor];
        DonorEnd[iDim] = Donors->Extents_[1][iDim][iDonor];
      }
      range DonorRange(NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = RangeIncludes(GlobalRange, DonorRange);
      for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
        if (AwayFromEdge) {
          RecvFromNeighborRanges[iNeighbor] = UnionRanges(RecvFromNeighborRanges[iNeighbor],
            IntersectRanges(GridNeighbors[iNeighbor].LocalRange, DonorRange));
        } else {
          for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
            for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
              for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                CartPeriodicAdjust(Cart, Point, Point);
                if (RangeContains(GridNeighbors[iNeighbor].LocalRange, Point)) {
                  ExtendRange(RecvFromNeighborRanges[iNeighbor], Point);
                }
              }
            }
          }
        }
      }
    }

    std::vector<int> CollectRecvIndexToNeighbor;
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!RecvFromNeighborRanges[iNeighbor].Empty()) {
        CollectRecvIndexToNeighbor.push_back(iNeighbor);
      }
    }
    int NumCollectRecvs = CollectRecvIndexToNeighbor.size();

    Exchange.CollectRecvs_.resize(NumCollectRecvs);

    for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      int iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      Exchange.CollectRecvs_[iCollectRecv].Rank = GridNeighbors[iNeighbor].Rank;
    }

    std::vector<std::vector<char>> CollectRecvMasks(NumCollectRecvs);
    for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      int iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      CollectRecvMasks[iCollectRecv].resize(RecvFromNeighborRanges[iNeighbor].Count(), false);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS], DonorEnd[MAX_DIMS];
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorBegin[iDim] = Donors->Extents_[0][iDim][iDonor];
        DonorEnd[iDim] = Donors->Extents_[1][iDim][iDonor];
      }
      range DonorRange(NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = RangeIncludes(GlobalRange, DonorRange);
      for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
        int iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
        if (AwayFromEdge) {
          range RemoteDonorRange = IntersectRanges(GridNeighbors[iNeighbor].LocalRange, DonorRange);
          for (int k = RemoteDonorRange.Begin(2); k < RemoteDonorRange.End(2); ++k) {
            for (int j = RemoteDonorRange.Begin(1); j < RemoteDonorRange.End(1); ++j) {
              for (int i = RemoteDonorRange.Begin(0); i < RemoteDonorRange.End(0); ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                long long iPoint = RangeTupleToIndex(RecvFromNeighborRanges[iNeighbor],
                  array_layout::COLUMN_MAJOR, Point);
                CollectRecvMasks[iCollectRecv][iPoint] = true;
              }
            }
          }
        } else {
          for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
            for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
              for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                CartPeriodicAdjust(Cart, Point, Point);
                if (RangeContains(GridNeighbors[iNeighbor].LocalRange, Point)) {
                  long long iPoint = RangeTupleToIndex(RecvFromNeighborRanges[iNeighbor],
                    array_layout::COLUMN_MAJOR, Point);
                  CollectRecvMasks[iCollectRecv][iPoint] = true;
                }
              }
            }
          }
        }
      }
    }

    for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      exchange::collect_recv &CollectRecv = Exchange.CollectRecvs_[iCollectRecv];
      int iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      CollectRecv.NumPoints = 0;
      for (long long iPoint = 0; iPoint < RecvFromNeighborRanges[iNeighbor].Count(); ++iPoint) {
        if (CollectRecvMasks[iCollectRecv][iPoint]) {
          ++CollectRecv.NumPoints;
        }
      }
    }

    for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      exchange::collect_recv &CollectRecv = Exchange.CollectRecvs_[iCollectRecv];
      int iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      CollectRecv.PointsData.resize(MAX_DIMS*CollectRecv.NumPoints);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        CollectRecv.Points[iDim] = CollectRecv.PointsData.data() + iDim*CollectRecv.NumPoints;
      }
      long long iCollectRecvPoint = 0;
      for (long long iPoint = 0; iPoint < RecvFromNeighborRanges[iNeighbor].Count(); ++iPoint) {
        if (CollectRecvMasks[iCollectRecv][iPoint]) {
          std::array<int,MAX_DIMS> Point = RangeIndexToTuple(RecvFromNeighborRanges[iNeighbor],
            array_layout::COLUMN_MAJOR, iPoint);
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            CollectRecv.Points[iDim][iCollectRecvPoint] = Point[iDim];
          }
          ++iCollectRecvPoint;
        }
      }

    }

    Exchange.NumRemoteDonorPoints_.resize(NumDonors, 0);
    Exchange.RemoteDonorPointsPtrs_.resize(NumDonors, nullptr);
    Exchange.RemoteDonorPoints_ = Exchange.RemoteDonorPointsPtrs_.data();
    Exchange.RemoteDonorPointCollectRecvsPtrs_.resize(NumDonors, nullptr);
    Exchange.RemoteDonorPointCollectRecvs_ = Exchange.RemoteDonorPointCollectRecvsPtrs_.data();
    Exchange.RemoteDonorPointCollectRecvBufferIndicesPtrs_.resize(NumDonors, nullptr);
    Exchange.RemoteDonorPointCollectRecvBufferIndices_ =
      Exchange.RemoteDonorPointCollectRecvBufferIndicesPtrs_.data();

    long long TotalRemoteDonorPoints = 0;
    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS], DonorEnd[MAX_DIMS];
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorBegin[iDim] = Donors->Extents_[0][iDim][iDonor];
        DonorEnd[iDim] = Donors->Extents_[1][iDim][iDonor];
      }
      range DonorRange(NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = RangeIncludes(GlobalRange, DonorRange);
      int NumRemoteDonorPoints;
      if (AwayFromEdge) {
        range LocalDonorRange = IntersectRanges(LocalRange, DonorRange);
        NumRemoteDonorPoints = DonorRange.Count() - LocalDonorRange.Count();
      } else {
        NumRemoteDonorPoints = 0;
        for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
          for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
            for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
              int Point[MAX_DIMS] = {i, j, k};
              CartPeriodicAdjust(Cart, Point, Point);
              if (!RangeContains(LocalRange, Point)) {
                ++NumRemoteDonorPoints;
              }
            }
          }
        }
      }
      Exchange.NumRemoteDonorPoints_[iDonor] = NumRemoteDonorPoints;
      TotalRemoteDonorPoints += NumRemoteDonorPoints;
    }

    Exchange.RemoteDonorPointsData_.resize(TotalRemoteDonorPoints);
    Exchange.RemoteDonorPointCollectRecvsData_.resize(TotalRemoteDonorPoints);
    Exchange.RemoteDonorPointCollectRecvBufferIndicesData_.resize(TotalRemoteDonorPoints);

    long long RemoteDonorPointOffset = 0;
    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      long long NumRemoteDonorPoints = Exchange.NumRemoteDonorPoints_[iDonor];
      Exchange.RemoteDonorPoints_[iDonor] = Exchange.RemoteDonorPointsData_.data() +
        RemoteDonorPointOffset;
      Exchange.RemoteDonorPointCollectRecvs_[iDonor] =
        Exchange.RemoteDonorPointCollectRecvsData_.data() + RemoteDonorPointOffset;
      Exchange.RemoteDonorPointCollectRecvBufferIndices_[iDonor] =
        Exchange.RemoteDonorPointCollectRecvBufferIndicesData_.data() + RemoteDonorPointOffset;
      RemoteDonorPointOffset += NumRemoteDonorPoints;
    }

    std::vector<std::vector<long long>> CollectRecvBufferIndices(NumCollectRecvs);
    for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      int iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      long long NumPoints = RecvFromNeighborRanges[iNeighbor].Count();
      CollectRecvBufferIndices[iCollectRecv].resize(NumPoints, -1);
      long long iRemotePoint = 0;
      for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
        if (CollectRecvMasks[iCollectRecv][iPoint]) {
          CollectRecvBufferIndices[iCollectRecv][iPoint] = iRemotePoint;
          ++iRemotePoint;
        }
      }
    }

    int MaxPointsInCell = 1;
    for (int iDim = 0; iDim <= NumDims; ++iDim) {
      MaxPointsInCell *= MaxSize;
    }
    std::vector<int> CellCollectRecvs(MaxPointsInCell);
    std::vector<long long> CellCollectRecvBufferIndices(MaxPointsInCell);

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      for (int iPointInCell = 0; iPointInCell < MaxPointsInCell; ++iPointInCell) {
        CellCollectRecvs[iPointInCell] = -1;
        CellCollectRecvBufferIndices[iPointInCell] = -1;
      }
      int DonorBegin[MAX_DIMS], DonorEnd[MAX_DIMS];
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorBegin[iDim] = Donors->Extents_[0][iDim][iDonor];
        DonorEnd[iDim] = Donors->Extents_[1][iDim][iDonor];
      }
      range DonorRange(NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = RangeIncludes(GlobalRange, DonorRange);
      for (int iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
        int iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
        if (AwayFromEdge) {
          range RemoteDonorRange = IntersectRanges(GridNeighbors[iNeighbor].LocalRange, DonorRange);
          for (int k = RemoteDonorRange.Begin(2); k < RemoteDonorRange.End(2); ++k) {
            for (int j = RemoteDonorRange.Begin(1); j < RemoteDonorRange.End(1); ++j) {
              for (int i = RemoteDonorRange.Begin(0); i < RemoteDonorRange.End(0); ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                int iPointInCell = RangeTupleToIndex<int>(DonorRange, array_layout::COLUMN_MAJOR,
                  Point);
                long long iPoint = RangeTupleToIndex(RecvFromNeighborRanges[iNeighbor],
                  array_layout::COLUMN_MAJOR, Point);
                CellCollectRecvs[iPointInCell] = iCollectRecv;
                CellCollectRecvBufferIndices[iPointInCell] = CollectRecvBufferIndices[iCollectRecv]
                  [iPoint];
              }
            }
          }
        } else {
          int iPointInCell = 0;
          for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
            for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
              for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                CartPeriodicAdjust(Cart, Point, Point);
                if (RangeContains(GridNeighbors[iNeighbor].LocalRange, Point)) {
                  long long iPoint = RangeTupleToIndex(RecvFromNeighborRanges[iNeighbor],
                    array_layout::COLUMN_MAJOR, Point);
                  CellCollectRecvs[iPointInCell] = iCollectRecv;
                  CellCollectRecvBufferIndices[iPointInCell] = CollectRecvBufferIndices[iCollectRecv]
                    [iPoint];
                }
                ++iPointInCell;
              }
            }
          }
        }
      }
      int iRemoteDonorPoint = 0;
      int iPointInCell = 0;
      for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
        for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
          for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
            if (CellCollectRecvs[iPointInCell] >= 0) {
              Exchange.RemoteDonorPoints_[iDonor][iRemoteDonorPoint] = iPointInCell;
              Exchange.RemoteDonorPointCollectRecvs_[iDonor][iRemoteDonorPoint] =
                CellCollectRecvs[iPointInCell];
              Exchange.RemoteDonorPointCollectRecvBufferIndices_[iDonor][iRemoteDonorPoint]
                = CellCollectRecvBufferIndices[iPointInCell];
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

  Exchange.DonorsSorted_.clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  if (RankHasConnectivityDonorSide(Connectivity)) {

    const connectivity_d *Donors;
    GetConnectivityDonorSide(Connectivity, Donors);

    long long NumDonors;
    GetConnectivityDonorSideCount(*Donors, NumDonors);

    if (NumDonors > 0) {

      Exchange.DonorsSorted_.resize(NumDonors);

      const grid_info *ReceiverGridInfo;
      GetConnectivityReceiverGridInfo(Connectivity, ReceiverGridInfo);

      range ReceiverGridGlobalRange;
      GetGridInfoGlobalRange(*ReceiverGridInfo, ReceiverGridGlobalRange);

      std::vector<long long> DestinationIndices(NumDonors);

      for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
        int DestinationPoint[MAX_DIMS] = {
          Donors->Destinations_[0][iDonor],
          Donors->Destinations_[1][iDonor],
          Donors->Destinations_[2][iDonor]
        };
        DestinationIndices[iDonor] = RangeTupleToIndex(ReceiverGridGlobalRange,
          array_layout::COLUMN_MAJOR, DestinationPoint);
      }

      bool Sorted = true;

      // Check if they're already sorted
      long long PrevIndex = 0;
      for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
        if (DestinationIndices[iDonor] < PrevIndex) {
          Sorted = false;
          break;
        }
        PrevIndex = DestinationIndices[iDonor];
      }

      if (Sorted) {
        for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
          Exchange.DonorsSorted_[iDonor] = iDonor;
        }
      } else {
        core::SortPermutation(NumDonors, DestinationIndices.data(), Exchange.DonorsSorted_.data());
      }

    }

  }

}

void UpdateReceiversSorted(exchange &Exchange) {

  Exchange.ReceiversSorted_.clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  if (RankHasConnectivityReceiverSide(Connectivity)) {

    const connectivity_r *Receivers;
    GetConnectivityReceiverSide(Connectivity, Receivers);

    long long NumReceivers = 0;
    GetConnectivityReceiverSideCount(*Receivers, NumReceivers);

    if (NumReceivers > 0) {

      Exchange.ReceiversSorted_.resize(NumReceivers);

      const grid *ReceiverGrid;
      GetConnectivityReceiverSideGrid(*Receivers, ReceiverGrid);

      range GlobalRange;
      GetGridGlobalRange(*ReceiverGrid, GlobalRange);

      std::vector<long long> PointIndices(NumReceivers);

      for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
        int Point[MAX_DIMS] = {
          Receivers->Points_[0][iReceiver],
          Receivers->Points_[1][iReceiver],
          Receivers->Points_[2][iReceiver]
        };
        PointIndices[iReceiver] = RangeTupleToIndex(GlobalRange, array_layout::COLUMN_MAJOR, Point);
      }

      bool Sorted = true;

      // Check if they're already sorted
      long long PrevIndex = 0;
      for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
        if (PointIndices[iReceiver] < PrevIndex) {
          Sorted = false;
          break;
        }
        PrevIndex = PointIndices[iReceiver];
      }

      if (Sorted) {
        for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
          Exchange.ReceiversSorted_[iReceiver] = iReceiver;
        }
      } else {
        core::SortPermutation(NumReceivers, PointIndices.data(), Exchange.ReceiversSorted_.data());
      }

    }

  }

}

void UpdateDestRanks(exchange &Exchange) {

  Exchange.DonorDestRanks_.clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  bool DonorGridIsLocal = RankHasConnectivityDonorSide(Connectivity);

  long long NumDonors = 0;

  const connectivity_d *Donors;
  if (DonorGridIsLocal) {
    GetConnectivityDonorSide(Connectivity, Donors);
    GetConnectivityDonorSideCount(*Donors, NumDonors);
  }

  Exchange.DonorDestRanks_.resize(NumDonors, -1);

  std::map<int, core::partition_bin> Bins;

  std::vector<int> DestinationBinIndices;
  if (DonorGridIsLocal) {
    DestinationBinIndices.resize(NumDonors);
    core::MapToPartitionBins(Exchange.DestinationHash_, NumDonors, Donors->Destinations_,
      DestinationBinIndices.data());
    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int BinIndex = DestinationBinIndices[iDonor];
      auto Iter = Bins.lower_bound(BinIndex);
      if (Iter == Bins.end() || Iter->first > BinIndex) {
        Bins.emplace_hint(Iter, BinIndex, core::partition_bin());
      }
    }
  }

  core::RetrievePartitionBins(Exchange.DestinationHash_, Bins);

  if (DonorGridIsLocal) {
    core::FindPartitions(Exchange.DestinationHash_, Bins, NumDonors, Donors->Destinations_,
      DestinationBinIndices.data(), Exchange.DonorDestRanks_.data());
  }

}

void UpdateSourceRanks(exchange &Exchange) {

  Exchange.ReceiverSourceRanks_.clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  bool ReceiverGridIsLocal = RankHasConnectivityReceiverSide(Connectivity);

  long long NumReceivers = 0;

  const connectivity_r *Receivers;
  if (ReceiverGridIsLocal) {
    GetConnectivityReceiverSide(Connectivity, Receivers);
    GetConnectivityReceiverSideCount(*Receivers, NumReceivers);
  }

  Exchange.ReceiverSourceRanks_.resize(NumReceivers, -1);

  std::map<int, core::partition_bin> Bins;

  std::vector<int> SourceBinIndices;
  if (ReceiverGridIsLocal) {
    SourceBinIndices.resize(NumReceivers);
    core::MapToPartitionBins(Exchange.SourceHash_, NumReceivers, Receivers->Sources_,
      SourceBinIndices.data());
    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      int BinIndex = SourceBinIndices[iReceiver];
      auto Iter = Bins.lower_bound(BinIndex);
      if (Iter == Bins.end() || Iter->first > BinIndex) {
        Bins.emplace_hint(Iter, BinIndex, core::partition_bin());
      }
    }
  }

  core::RetrievePartitionBins(Exchange.SourceHash_, Bins);

  if (ReceiverGridIsLocal) {
    core::FindPartitions(Exchange.SourceHash_, Bins, NumReceivers, Receivers->Sources_,
      SourceBinIndices.data(), Exchange.ReceiverSourceRanks_.data());
  }

}

void UpdateSendInfo(exchange &Exchange) {

  Exchange.Sends_.clear();

  Exchange.DonorSendIndices_.clear();

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

    Exchange.DonorSendIndices_.resize(NumDonors, -1);

    std::vector<char> DonorCommunicates(NumDonors);

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      bool Communicates = Exchange.DonorDestRanks_[iDonor] >= 0;
      if (Communicates) {
        int DonorCell[MAX_DIMS] = {
          Donors->Extents_[0][0][iDonor],
          Donors->Extents_[0][1][iDonor],
          Donors->Extents_[0][2][iDonor]
        };
        CartPeriodicAdjust(Cart, DonorCell, DonorCell);
        Communicates = RangeContains(LocalRange, DonorCell);
      }
      DonorCommunicates[iDonor] = Communicates;
    }

    std::map<int, long long> SendCounts;

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      if (DonorCommunicates[iDonor]) {
        int DestRank = Exchange.DonorDestRanks_[iDonor];
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
      Exchange.Sends_.emplace_back();
      exchange::send &Send = Exchange.Sends_.back();
      Send.Rank = Pair.first;
      Send.Count = Pair.second;
    }

    SendCounts.clear();

    std::map<int, int> RankToSendIndex;

    for (int iSend = 0; iSend < NumSends; ++iSend) {
      RankToSendIndex.emplace(Exchange.Sends_[iSend].Rank, iSend);
    }

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      if (DonorCommunicates[iDonor]) {
        int DestRank = Exchange.DonorDestRanks_[iDonor];
        Exchange.DonorSendIndices_[iDonor] = RankToSendIndex[DestRank];
      }
    }

  }

}

void UpdateReceiveInfo(exchange &Exchange) {

  Exchange.Recvs_.clear();

  Exchange.ReceiverRecvIndices_.clear();

  const connectivity &Connectivity = *Exchange.Connectivity_;

  const connectivity_r *Receivers;
  long long NumReceivers = 0;
  if (RankHasConnectivityReceiverSide(Connectivity)) {
    GetConnectivityReceiverSide(Connectivity, Receivers);
    GetConnectivityReceiverSideCount(*Receivers, NumReceivers);
  }

  if (NumReceivers > 0) {

    Exchange.ReceiverRecvIndices_.resize(NumReceivers, -1);

    std::map<int, long long> RecvCounts;

    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      if (Exchange.ReceiverSourceRanks_[iReceiver] >= 0) {
        int SourceRank = Exchange.ReceiverSourceRanks_[iReceiver];
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
      Exchange.Recvs_.emplace_back();
      exchange::recv &Recv = Exchange.Recvs_.back();
      Recv.Rank = Pair.first;
      Recv.Count = Pair.second;
    }

    RecvCounts.clear();

    std::map<int, int> RankToRecvIndex;

    for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
      RankToRecvIndex.emplace(Exchange.Recvs_[iRecv].Rank, iRecv);
    }

    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      int SourceRank = Exchange.ReceiverSourceRanks_[iReceiver];
      if (SourceRank >= 0) {
        Exchange.ReceiverRecvIndices_[iReceiver] = RankToRecvIndex[SourceRank];
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

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {
    Collect_->Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);
  }

  void Collect(const void * const *GridValues, void **DonorValues) {
    Collect_->Collect(GridValues, DonorValues);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange,
      array_layout GridValuesLayout) = 0;
    virtual void Collect(const void * const *GridValues, void **DonorValues) = 0;
  };

  template <typename T> class model : public concept {
  public:
    using user_value_type = typename T::user_value_type;
    explicit model(T Collect):
      Collect_(std::move(Collect))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange,
      array_layout GridValuesLayout) override {
      Collect_.Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);
      GridValues_.resize(Count);
      DonorValues_.resize(Count);
    }
    virtual void Collect(const void * const *GridValuesVoid, void **DonorValuesVoid) override {
      OVK_DEBUG_ASSERT(GridValuesVoid || GridValues_.size() == 0, "Invalid grid values pointer.");
      OVK_DEBUG_ASSERT(DonorValuesVoid || DonorValues_.size() == 0, "Invalid donor values pointer.");
      for (int iCount = 0; iCount < int(GridValues_.size()); ++iCount) {
        GridValues_[iCount] = static_cast<const user_value_type *>(GridValuesVoid[iCount]);
      }
      for (int iCount = 0; iCount < int(DonorValues_.size()); ++iCount) {
        DonorValues_[iCount] = static_cast<user_value_type *>(DonorValuesVoid[iCount]);
      }
      Collect_.Collect(GridValues_.data(), DonorValues_.data());
    }
    T Collect_;
    std::vector<const user_value_type *> GridValues_;
    std::vector<user_value_type *> DonorValues_;
  };

  std::unique_ptr<concept> Collect_;

};

// Use unsigned char in place of bool for MPI sends/recvs
namespace no_bool_internal {
template <typename T> struct helper { using type = T; };
template <> struct helper<bool> { using type = unsigned char; };
}
template <typename T> using no_bool = typename no_bool_internal::helper<T>::type;

template <typename T> class collect_base {

public:

  using user_value_type = T;
  using value_type = no_bool<T>;

  collect_base() = default;
  collect_base(const collect_base &) = delete;
  collect_base(collect_base &&) = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {

    NumDims_ = Exchange.NumDims_;
    Count_ = Count;
    GridValuesRange_ = GridValuesRange;
    GridValuesLayout_ = GridValuesLayout;

    const connectivity &Connectivity = *Exchange.Connectivity_;

    GetConnectivityDonorSide(Connectivity, Donors_);
    const connectivity_d &Donors = *Donors_;

    GetConnectivityDonorSideGrid(Donors, Grid_);
    const grid &Grid = *Grid_;

    GetGridCart(Grid, Cart_);
    GetGridGlobalRange(Grid, GlobalRange_);
    GetGridLocalRange(Grid, LocalRange_);

    OVK_DEBUG_ASSERT(RangeIncludes(GridValuesRange, LocalRange_), "Invalid grid values range.");

    long long NumDonors;
    int MaxSize;
    GetConnectivityDonorSideCount(Donors, NumDonors);
    GetConnectivityDonorSideMaxSize(Donors, MaxSize);

    MaxPointsInCell_ = 1;
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      MaxPointsInCell_ *= MaxSize;
    }

    NumSends_ = Exchange.CollectSends_.size();
    Sends_ = Exchange.CollectSends_.data();
    NumRecvs_ = Exchange.CollectRecvs_.size();
    Recvs_ = Exchange.CollectRecvs_.data();

    SendBuffers_.resize(NumSends_);
    for (int iSend = 0; iSend < NumSends_; ++iSend) {
      SendBuffers_[iSend].resize(Count_);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        SendBuffers_[iSend][iCount].resize(Exchange.CollectSends_[iSend].NumPoints);
      }
    }

    Requests_.resize(Count_*(NumSends_+NumRecvs_));

    NumRemoteDonorPoints_ = Exchange.NumRemoteDonorPoints_.data();
    RemoteDonorPoints_ = Exchange.RemoteDonorPoints_;
    RemoteDonorPointCollectRecvs_ = Exchange.RemoteDonorPointCollectRecvs_;
    RemoteDonorPointCollectRecvBufferIndices_ = Exchange.RemoteDonorPointCollectRecvBufferIndices_;

  }

protected:

  using remote_donor_values = std::vector<std::vector<std::vector<value_type>>>;

  const connectivity_d *Donors_;
  const grid *Grid_;
  int Count_;
  range GridValuesRange_;
  array_layout GridValuesLayout_;
  int MaxPointsInCell_;

  void CheckValuesPointers(const user_value_type * const *GridValues, const user_value_type * const
    *DonorValues) {

    const connectivity_d &Donors = *Donors_;
    long long NumDonors;
    GetConnectivityDonorSideCount(Donors, NumDonors);

    if (NumDonors > 0) {
      for (int iCount = 0; iCount < Count_; ++iCount) {
        OVK_DEBUG_ASSERT(GridValues[iCount], "Invalid grid values pointer.");
      }
      for (int iCount = 0; iCount < Count_; ++iCount) {
        OVK_DEBUG_ASSERT(DonorValues[iCount], "Invalid donor values pointer.");
      }
    }

  }

  void AllocateRemoteDonorValues(remote_donor_values &RemoteDonorValues) {

    RemoteDonorValues.resize(NumRecvs_);
    for (int iRecv = 0; iRecv < NumRecvs_; ++iRecv) {
      RemoteDonorValues[iRecv].resize(Count_);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        RemoteDonorValues[iRecv][iCount].resize(Recvs_[iRecv].NumPoints);
      }
    }
    
  }

  void RetrieveRemoteDonorValues(const user_value_type * const *GridValues, remote_donor_values
    &RemoteDonorValues) {

    MPI_Comm Comm;
    GetGridComm(*Grid_, Comm);

    MPI_Datatype MPIDataType = core::GetMPIDataType<value_type>();

    int iRequest = 0;

    for (int iRecv = 0; iRecv < NumRecvs_; ++iRecv) {
      const exchange::collect_recv &Recv = Recvs_[iRecv];
      for (int iCount = 0; iCount < Count_; ++iCount) {
        MPI_Irecv(RemoteDonorValues[iRecv][iCount].data(), Recv.NumPoints, MPIDataType, Recv.Rank,
          iCount, Comm, &Requests_[iRequest]);
        ++iRequest;
      }
    }

    for (int iSend = 0; iSend < NumSends_; ++iSend) {
      const exchange::collect_send &Send = Sends_[iSend];
      for (long long iSendPoint = 0; iSendPoint < Send.NumPoints; ++iSendPoint) {
        int Point[MAX_DIMS] = {
          Send.Points[0][iSendPoint],
          Send.Points[1][iSendPoint],
          Send.Points[2][iSendPoint]
        };
        long long iGridPoint = RangeTupleToIndex(GridValuesRange_, GridValuesLayout_, Point);
        for (int iCount = 0; iCount < Count_; ++iCount) {
          SendBuffers_[iSend][iCount][iSendPoint] = value_type(GridValues[iCount][iGridPoint]);
        }
      }
    }

    for (int iSend = 0; iSend < NumSends_; ++iSend) {
      const exchange::collect_send &Send = Sends_[iSend];
      for (int iCount = 0; iCount < Count_; ++iCount) {
        MPI_Isend(SendBuffers_[iSend][iCount].data(), Send.NumPoints, MPIDataType, Send.Rank, iCount,
          Comm, &Requests_[iRequest]);
        ++iRequest;
      }
    }

    MPI_Waitall(Requests_.size(), Requests_.data(), MPI_STATUSES_IGNORE);

  }

  void AssembleDonorPointValues(const user_value_type * const *GridValues, const remote_donor_values
    &RemoteDonorValues, long long iDonor, int *DonorSize, value_type *DonorPointValues) {

    const connectivity_d &Donors = *Donors_;

    int DonorBegin[MAX_DIMS], DonorEnd[MAX_DIMS];
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      DonorBegin[iDim] = Donors.Extents_[0][iDim][iDonor];
      DonorEnd[iDim] = Donors.Extents_[1][iDim][iDonor];
    }

    range DonorRange(NumDims_, DonorBegin, DonorEnd);

    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      DonorSize[iDim] = DonorRange.Size(iDim);
    }
    int NumPointsInCell = DonorRange.Count();

    bool AwayFromEdge = RangeIncludes(GlobalRange_, DonorRange);

    // Fill in the local data
    if (AwayFromEdge) {
      range LocalDonorRange = IntersectRanges(LocalRange_, DonorRange);
      for (int k = LocalDonorRange.Begin(2); k < LocalDonorRange.End(2); ++k) {
        for (int j = LocalDonorRange.Begin(1); j < LocalDonorRange.End(1); ++j) {
          for (int i = LocalDonorRange.Begin(0); i < LocalDonorRange.End(0); ++i) {
            int Point[MAX_DIMS] = {i, j, k};
            int iPointInCell = RangeTupleToIndex<int>(DonorRange, array_layout::COLUMN_MAJOR,
              Point);
            long long iGridPoint = RangeTupleToIndex(GridValuesRange_, GridValuesLayout_, Point);
            for (int iCount = 0; iCount < Count_; ++iCount) {
              DonorPointValues[iPointInCell+iCount*NumPointsInCell] = value_type(GridValues[iCount]
                [iGridPoint]);
            }
          }
        }
      }
    } else {
      int iPointInCell = 0;
      for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
        for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
          for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
            int Point[MAX_DIMS] = {i, j, k};
            CartPeriodicAdjust(Cart_, Point, Point);
            if (RangeContains(LocalRange_, Point)) {
              long long iGridPoint = RangeTupleToIndex(GridValuesRange_, GridValuesLayout_, Point);
              for (int iCount = 0; iCount < Count_; ++iCount) {
                DonorPointValues[iPointInCell+iCount*NumPointsInCell] = value_type(GridValues[iCount]
                  [iGridPoint]);
              }
            }
            ++iPointInCell;
          }
        }
      }
    }

    // Fill in the remote data
    for (int iRemotePoint = 0; iRemotePoint < NumRemoteDonorPoints_[iDonor]; ++iRemotePoint) {
      int iPointInCell = RemoteDonorPoints_[iDonor][iRemotePoint];
      int iRecv = RemoteDonorPointCollectRecvs_[iDonor][iRemotePoint];
      long long iBuffer = RemoteDonorPointCollectRecvBufferIndices_[iDonor][iRemotePoint];
      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorPointValues[iPointInCell+iCount*NumPointsInCell] = RemoteDonorValues[iRecv][iCount]
          [iBuffer];
      }
    }

  }

private:

  int NumDims_;
  cart Cart_;
  range GlobalRange_;
  range LocalRange_;
  int NumSends_;
  const exchange::collect_send *Sends_;
  int NumRecvs_;
  const exchange::collect_recv *Recvs_;
  std::vector<std::vector<std::vector<value_type>>> SendBuffers_;
  std::vector<MPI_Request> Requests_;
  const int *NumRemoteDonorPoints_;
  const long long * const *RemoteDonorPoints_;
  const int * const *RemoteDonorPointCollectRecvs_;
  const long long * const *RemoteDonorPointCollectRecvBufferIndices_;
  
};

template <typename T> class collect_none : public collect_base<T> {

protected:

  using parent_type = collect_base<T>;

  using parent_type::Donors_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesLayout_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::user_value_type;
  using typename parent_type::value_type;

  collect_none() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {

    parent_type::Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.resize(Count_*MaxPointsInCell_);

  }

  void Collect(const user_value_type * const *GridValues, user_value_type **DonorValues) {

    const connectivity_d &Donors = *Donors_;

    int NumDims;
    long long NumDonors;
    GetConnectivityDonorSideDimension(Donors, NumDims);
    GetConnectivityDonorSideCount(Donors, NumDonors);

    if (OVK_DEBUG) parent_type::CheckValuesPointers(GridValues, DonorValues);

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {

      int DonorSize[MAX_DIMS];
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_.data());
      int NumPointsInCell = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues[iCount][iDonor] = user_value_type(true);
      }
      for (int k = 0; k < DonorSize[2]; ++k) {
        for (int j = 0; j < DonorSize[1]; ++j) {
          for (int i = 0; i < DonorSize[0]; ++i) {
            int iPoint = i + DonorSize[0]*(j + DonorSize[1]*k);
            for (int iCount = 0; iCount < Count_; ++iCount) {
              DonorValues[iCount][iDonor] = DonorValues[iCount][iDonor] &&
                !user_value_type(DonorPointValues_[iCount*NumPointsInCell+iPoint]);
            }
          }
        }
      }

    }

  }

private:

  typename parent_type::remote_donor_values RemoteDonorValues_;
  std::vector<value_type> DonorPointValues_;

};

template <typename T> class collect_any : public collect_base<T> {

protected:

  using parent_type = collect_base<T>;

  using parent_type::Donors_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesLayout_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::user_value_type;
  using typename parent_type::value_type;

  collect_any() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {

    parent_type::Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.resize(Count_*MaxPointsInCell_);

  }

  void Collect(const user_value_type * const *GridValues, user_value_type **DonorValues) {

    const connectivity_d &Donors = *Donors_;

    int NumDims;
    long long NumDonors;
    GetConnectivityDonorSideDimension(Donors, NumDims);
    GetConnectivityDonorSideCount(Donors, NumDonors);

    if (OVK_DEBUG) parent_type::CheckValuesPointers(GridValues, DonorValues);

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {

      int DonorSize[MAX_DIMS];
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_.data());
      int NumPointsInCell = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues[iCount][iDonor] = user_value_type(false);
      }
      for (int k = 0; k < DonorSize[2]; ++k) {
        for (int j = 0; j < DonorSize[1]; ++j) {
          for (int i = 0; i < DonorSize[0]; ++i) {
            int iPoint = i + DonorSize[0]*(j + DonorSize[1]*k);
            for (int iCount = 0; iCount < Count_; ++iCount) {
              DonorValues[iCount][iDonor] = DonorValues[iCount][iDonor] ||
                user_value_type(DonorPointValues_[iCount*NumPointsInCell+iPoint]);
            }
          }
        }
      }

    }

  }

private:

  typename parent_type::remote_donor_values RemoteDonorValues_;
  std::vector<value_type> DonorPointValues_;

};

template <typename T> class collect_not_all : public collect_base<T> {

protected:

  using parent_type = collect_base<T>;

  using parent_type::Donors_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesLayout_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::user_value_type;
  using typename parent_type::value_type;

  collect_not_all() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {

    parent_type::Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.resize(Count_*MaxPointsInCell_);

  }

  void Collect(const user_value_type * const *GridValues, user_value_type **DonorValues) {

    const connectivity_d &Donors = *Donors_;

    int NumDims;
    long long NumDonors;
    GetConnectivityDonorSideDimension(Donors, NumDims);
    GetConnectivityDonorSideCount(Donors, NumDonors);

    if (OVK_DEBUG) parent_type::CheckValuesPointers(GridValues, DonorValues);

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {

      int DonorSize[MAX_DIMS];
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_.data());
      int NumPointsInCell = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues[iCount][iDonor] = user_value_type(false);
      }
      for (int k = 0; k < DonorSize[2]; ++k) {
        for (int j = 0; j < DonorSize[1]; ++j) {
          for (int i = 0; i < DonorSize[0]; ++i) {
            int iPoint = i + DonorSize[0]*(j + DonorSize[1]*k);
            for (int iCount = 0; iCount < Count_; ++iCount) {
              DonorValues[iCount][iDonor] = DonorValues[iCount][iDonor] ||
                !user_value_type(DonorPointValues_[iCount*NumPointsInCell+iPoint]);
            }
          }
        }
      }

    }

  }

private:

  typename parent_type::remote_donor_values RemoteDonorValues_;
  std::vector<value_type> DonorPointValues_;

};

template <typename T> class collect_all : public collect_base<T> {

protected:

  using parent_type = collect_base<T>;

  using parent_type::Donors_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesLayout_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::user_value_type;
  using typename parent_type::value_type;

  collect_all() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {

    parent_type::Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.resize(Count_*MaxPointsInCell_);

  }

  void Collect(const user_value_type * const *GridValues, user_value_type **DonorValues) {

    const connectivity_d &Donors = *Donors_;

    int NumDims;
    long long NumDonors;
    GetConnectivityDonorSideDimension(Donors, NumDims);
    GetConnectivityDonorSideCount(Donors, NumDonors);

    if (OVK_DEBUG) parent_type::CheckValuesPointers(GridValues, DonorValues);

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {

      int DonorSize[MAX_DIMS];
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_.data());
      int NumPointsInCell = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues[iCount][iDonor] = user_value_type(true);
      }
      for (int k = 0; k < DonorSize[2]; ++k) {
        for (int j = 0; j < DonorSize[1]; ++j) {
          for (int i = 0; i < DonorSize[0]; ++i) {
            int iPoint = i + DonorSize[0]*(j + DonorSize[1]*k);
            for (int iCount = 0; iCount < Count_; ++iCount) {
              DonorValues[iCount][iDonor] = DonorValues[iCount][iDonor] &&
                user_value_type(DonorPointValues_[iCount*NumPointsInCell+iPoint]);
            }
          }
        }
      }

    }

  }

private:

  typename parent_type::remote_donor_values RemoteDonorValues_;
  std::vector<value_type> DonorPointValues_;

};

template <typename T> class collect_interp : public collect_base<T> {

protected:

  using parent_type = collect_base<T>;

  using parent_type::Donors_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesLayout_;
  using parent_type::MaxPointsInCell_;

public:

  using typename parent_type::user_value_type;
  using typename parent_type::value_type;

  collect_interp() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {

    parent_type::Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);

    InterpCoefs_ = Donors_->InterpCoefs_;

    parent_type::AllocateRemoteDonorValues(RemoteDonorValues_);

    DonorPointValues_.resize(Count_*MaxPointsInCell_);

  }

  void Collect(const user_value_type * const *GridValues, user_value_type **DonorValues) {

    const connectivity_d &Donors = *Donors_;

    int NumDims;
    long long NumDonors;
    GetConnectivityDonorSideDimension(Donors, NumDims);
    GetConnectivityDonorSideCount(Donors, NumDonors);

    if (OVK_DEBUG) parent_type::CheckValuesPointers(GridValues, DonorValues);

    parent_type::RetrieveRemoteDonorValues(GridValues, RemoteDonorValues_);

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {

      int DonorSize[MAX_DIMS];
      parent_type::AssembleDonorPointValues(GridValues, RemoteDonorValues_, iDonor, DonorSize,
        DonorPointValues_.data());
      int NumPointsInCell = DonorSize[0]*DonorSize[1]*DonorSize[2];

      for (int iCount = 0; iCount < Count_; ++iCount) {
        DonorValues[iCount][iDonor] = user_value_type(0);
      }
      for (int k = 0; k < DonorSize[2]; ++k) {
        for (int j = 0; j < DonorSize[1]; ++j) {
          for (int i = 0; i < DonorSize[0]; ++i) {
            int Point[MAX_DIMS] = {i, j, k};
            int iPoint = i + DonorSize[0]*(j + DonorSize[1]*k);
            double Coef = 1.;
            for (int iDim = 0; iDim < NumDims; ++iDim) {
              Coef *= InterpCoefs_[iDim][Point[iDim]][iDonor];
            }
            for (int iCount = 0; iCount < Count_; ++iCount) {
              DonorValues[iCount][iDonor] += user_value_type(Coef*DonorPointValues_[iCount*
                NumPointsInCell+iPoint]);
            }
          }
        }
      }

    }

  }

private:

  const double * const * const *InterpCoefs_;
  typename parent_type::remote_donor_values RemoteDonorValues_;
  std::vector<value_type> DonorPointValues_;

};

}

namespace core {

void Collect(const exchange &Exchange, data_type ValueType, int Count, collect_op CollectOp,
  const range &GridValuesRange, array_layout GridValuesLayout, const void * const *GridValues,
  void **DonorValues) {

  collect CollectWrapper;

  switch (CollectOp) {
  case collect_op::NONE:
    switch (ValueType) {
    case data_type::BOOL: CollectWrapper = collect_none<bool>(); break;
    case data_type::BYTE: CollectWrapper = collect_none<unsigned char>(); break;
    case data_type::INT: CollectWrapper = collect_none<int>(); break;
    case data_type::LONG: CollectWrapper = collect_none<long>(); break;
    case data_type::LONG_LONG: CollectWrapper = collect_none<long long>(); break;
    case data_type::UNSIGNED_INT: CollectWrapper = collect_none<unsigned int>(); break;
    case data_type::UNSIGNED_LONG: CollectWrapper = collect_none<unsigned long>(); break;
    case data_type::UNSIGNED_LONG_LONG: CollectWrapper = collect_none<unsigned long long>(); break;
    case data_type::FLOAT: CollectWrapper = collect_none<float>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_none<double>(); break;
    }
    break;
  case collect_op::ANY:
    switch (ValueType) {
    case data_type::BOOL: CollectWrapper = collect_any<bool>(); break;
    case data_type::BYTE: CollectWrapper = collect_any<unsigned char>(); break;
    case data_type::INT: CollectWrapper = collect_any<int>(); break;
    case data_type::LONG: CollectWrapper = collect_any<long>(); break;
    case data_type::LONG_LONG: CollectWrapper = collect_any<long long>(); break;
    case data_type::UNSIGNED_INT: CollectWrapper = collect_any<unsigned int>(); break;
    case data_type::UNSIGNED_LONG: CollectWrapper = collect_any<unsigned long>(); break;
    case data_type::UNSIGNED_LONG_LONG: CollectWrapper = collect_any<unsigned long long>(); break;
    case data_type::FLOAT: CollectWrapper = collect_any<float>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_any<double>(); break;
    }
    break;
  case collect_op::NOT_ALL:
    switch (ValueType) {
    case data_type::BOOL: CollectWrapper = collect_not_all<bool>(); break;
    case data_type::BYTE: CollectWrapper = collect_not_all<unsigned char>(); break;
    case data_type::INT: CollectWrapper = collect_not_all<int>(); break;
    case data_type::LONG: CollectWrapper = collect_not_all<long>(); break;
    case data_type::LONG_LONG: CollectWrapper = collect_not_all<long long>(); break;
    case data_type::UNSIGNED_INT: CollectWrapper = collect_not_all<unsigned int>(); break;
    case data_type::UNSIGNED_LONG: CollectWrapper = collect_not_all<unsigned long>(); break;
    case data_type::UNSIGNED_LONG_LONG: CollectWrapper = collect_not_all<unsigned long long>(); break;
    case data_type::FLOAT: CollectWrapper = collect_not_all<float>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_not_all<double>(); break;
    }
    break;
  case collect_op::ALL:
    switch (ValueType) {
    case data_type::BOOL: CollectWrapper = collect_all<bool>(); break;
    case data_type::BYTE: CollectWrapper = collect_all<unsigned char>(); break;
    case data_type::INT: CollectWrapper = collect_all<int>(); break;
    case data_type::LONG: CollectWrapper = collect_all<long>(); break;
    case data_type::LONG_LONG: CollectWrapper = collect_all<long long>(); break;
    case data_type::UNSIGNED_INT: CollectWrapper = collect_all<unsigned int>(); break;
    case data_type::UNSIGNED_LONG: CollectWrapper = collect_all<unsigned long>(); break;
    case data_type::UNSIGNED_LONG_LONG: CollectWrapper = collect_all<unsigned long long>(); break;
    case data_type::FLOAT: CollectWrapper = collect_all<float>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_all<double>(); break;
    }
    break;
  case collect_op::INTERPOLATE:
    switch (ValueType) {
    case data_type::FLOAT: CollectWrapper = collect_interp<float>(); break;
    case data_type::DOUBLE: CollectWrapper = collect_interp<double>(); break;
    default:
      OVK_DEBUG_ASSERT(false, "Invalid data type for interpolation collect operation.");
      break;
    }
    break;
  }

  CollectWrapper.Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);
  CollectWrapper.Collect(GridValues, DonorValues);

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
    using user_value_type = typename T::user_value_type;
    explicit model(T Send):
      Send_(std::move(Send))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, int Tag) override {
      Send_.Initialize(Exchange, Count, Tag);
      DonorValues_.resize(Count);
    }
    virtual void Send(const void * const *DonorValuesVoid, request &Request) override {
      OVK_DEBUG_ASSERT(DonorValuesVoid || DonorValues_.size() == 0, "Invalid donor values pointer.");
      for (int iCount = 0; iCount < int(DonorValues_.size()); ++iCount) {
        DonorValues_[iCount] = static_cast<const user_value_type *>(DonorValuesVoid[iCount]);
      }
      Send_.Send(DonorValues_.data(), Request);
    }
    T Send_;
    std::vector<const user_value_type *> DonorValues_;
  };

  std::unique_ptr<concept> Send_;

};

template <typename T> class send_request {
public:
  using user_value_type = T;
  using value_type = no_bool<T>;
  send_request(const exchange &Exchange, int Count, int NumSends,
    std::vector<std::vector<value_type>> Buffers, std::vector<MPI_Request> MPIRequests):
    Exchange_(&Exchange),
    Count_(Count),
    NumSends_(NumSends),
    Buffers_(std::move(Buffers)),
    MPIRequests_(std::move(MPIRequests))
  {
    const connectivity &Connectivity = *Exchange.Connectivity_;
    GetConnectivityDonorSide(Connectivity, Donors_);
  }
  void Wait();
  int NumMPIRequests() const { return MPIRequests_.size(); }
  MPI_Request *MPIRequests() { return MPIRequests_.data(); }
private:
  const exchange *Exchange_;
  const connectivity_d *Donors_;
  int Count_;
  int NumSends_;
  std::vector<std::vector<value_type>> Buffers_;
  std::vector<MPI_Request> MPIRequests_;
};

template <typename T> class send_impl {

public:

  using user_value_type = T;
  using value_type = no_bool<T>;
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

    NumSends_ = Exchange.Sends_.size();

  }

  void Send(const user_value_type * const *DonorValues, request &Request) {

    const exchange &Exchange = *Exchange_;
    const connectivity_d &Donors = *Donors_;

    long long NumDonors;
    GetConnectivityDonorSideCount(Donors, NumDonors);

    if (OVK_DEBUG) {
      if (NumDonors > 0) {
        for (int iCount = 0; iCount < Count_; ++iCount) {
          OVK_DEBUG_ASSERT(DonorValues[iCount], "Invalid donor values pointer.");
        }
      }
    }

    MPI_Datatype MPIDataType = core::GetMPIDataType<value_type>();

    std::vector<long long> NextBufferEntry(NumSends_, 0);

    std::vector<std::vector<value_type>> Buffers(NumSends_);

    for (int iSend = 0; iSend < NumSends_; ++iSend) {
      const exchange::send &Send = Exchange.Sends_[iSend];
      Buffers[iSend].resize(Count_*Send.Count);
    }

    for (long long iDonorOrder = 0; iDonorOrder < NumDonors; ++iDonorOrder) {
      long long iDonor = Exchange.DonorsSorted_[iDonorOrder];
      int iSend = Exchange.DonorSendIndices_[iDonor];
      if (iSend >= 0) {
        const exchange::send &Send = Exchange.Sends_[iSend];
        long long iBuffer = NextBufferEntry[iSend];
        for (int iCount = 0; iCount < Count_; ++iCount) {
          Buffers[iSend][iCount*Send.Count+iBuffer] = value_type(DonorValues[iCount][iDonor]);
        }
        ++NextBufferEntry[iSend];
      }
    }

    std::vector<MPI_Request> MPIRequests(NumSends_);

    for (int iSend = 0; iSend < NumSends_; ++iSend) {
      const exchange::send &Send = Exchange.Sends_[iSend];
      MPI_Isend(Buffers[iSend].data(), Count_*Send.Count, MPIDataType, Send.Rank, Tag_,
        Exchange.Comm_, &MPIRequests[iSend]);
    }

    Request = request_type(Exchange, Count_, NumSends_, std::move(Buffers), std::move(MPIRequests));

  }

private:

  const exchange *Exchange_;
  const connectivity_d *Donors_;
  int Count_;
  int Tag_;
  int NumSends_;

};

template <typename T> void send_request<T>::Wait() {

  MPI_Waitall(NumSends_, MPIRequests_.data(), MPI_STATUSES_IGNORE);

  Buffers_.clear();
  MPIRequests_.clear();

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
    using user_value_type = typename T::user_value_type;
    explicit model(T Recv):
      Recv_(std::move(Recv))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, int Tag) override {
      Recv_.Initialize(Exchange, Count, Tag);
      ReceiverValues_.resize(Count);
    }
    virtual void Recv(void **ReceiverValuesVoid, request &Request) override {
      OVK_DEBUG_ASSERT(ReceiverValuesVoid || ReceiverValues_.size() == 0, "Invalid receiver values "
        "pointer.");
      for (int iCount = 0; iCount < int(ReceiverValues_.size()); ++iCount) {
        ReceiverValues_[iCount] = static_cast<user_value_type *>(ReceiverValuesVoid[iCount]);
      }
      Recv_.Recv(ReceiverValues_.data(), Request);
    }
    T Recv_;
    std::vector<user_value_type *> ReceiverValues_;
  };

  std::unique_ptr<concept> Recv_;

};

template <typename T> class recv_request {
public:
  using user_value_type = T;
  using value_type = no_bool<T>;
  recv_request(const exchange &Exchange, int Count, int NumRecvs, 
    std::vector<std::vector<value_type>> Buffers, std::vector<MPI_Request> MPIRequests,
    std::vector<user_value_type *> ReceiverValues):
    Exchange_(&Exchange),
    Count_(Count),
    NumRecvs_(NumRecvs),
    Buffers_(std::move(Buffers)),
    MPIRequests_(std::move(MPIRequests)),
    ReceiverValues_(std::move(ReceiverValues))
  {
    const connectivity &Connectivity = *Exchange.Connectivity_;
    GetConnectivityReceiverSide(Connectivity, Receivers_);
  }
  void Wait();
  int NumMPIRequests() const { return MPIRequests_.size(); }
  MPI_Request *MPIRequests() { return MPIRequests_.data(); }
private:
  const exchange *Exchange_;
  const connectivity_r *Receivers_;
  int Count_;
  int NumRecvs_;
  std::vector<std::vector<value_type>> Buffers_;
  std::vector<MPI_Request> MPIRequests_;
  std::vector<user_value_type *> ReceiverValues_;
};

template <typename T> class recv_impl {

public:

  using user_value_type = T;
  using value_type = no_bool<T>;
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

    NumRecvs_ = Exchange.Recvs_.size();

  }

  void Recv(user_value_type **ReceiverValues, request &Request) {

    const exchange &Exchange = *Exchange_;
    const connectivity_r &Receivers = *Receivers_;

    long long NumReceivers;
    GetConnectivityReceiverSideCount(Receivers, NumReceivers);

    if (OVK_DEBUG) {
      if (NumReceivers > 0) {
        for (int iCount = 0; iCount < Count_; ++iCount) {
          OVK_DEBUG_ASSERT(ReceiverValues[iCount], "Invalid receiver data pointer.");
        }
      }
    }

    MPI_Datatype MPIDataType = core::GetMPIDataType<value_type>();

    std::vector<std::vector<value_type>> Buffers(NumRecvs_);

    for (int iRecv = 0; iRecv < NumRecvs_; ++iRecv) {
      const exchange::recv &Recv = Exchange.Recvs_[iRecv];
      Buffers[iRecv].resize(Count_*Recv.Count);
    }

    std::vector<MPI_Request> MPIRequests(NumRecvs_);

    for (int iRecv = 0; iRecv < NumRecvs_; ++iRecv) {
      const exchange::recv &Recv = Exchange.Recvs_[iRecv];
      MPI_Irecv(Buffers[iRecv].data(), Count_*Recv.Count, MPIDataType, Recv.Rank, Tag_,
        Exchange.Comm_, &MPIRequests[iRecv]);
    }

    std::vector<user_value_type *> ReceiverValuesSaved(ReceiverValues, ReceiverValues+Count_);

    Request = request_type(Exchange, Count_, NumRecvs_, std::move(Buffers), std::move(MPIRequests),
      std::move(ReceiverValuesSaved));

  }

private:

  const exchange *Exchange_;
  const connectivity_r *Receivers_;
  int Count_;
  int Tag_;
  int NumRecvs_;

};

template <typename T> void recv_request<T>::Wait() {

  const exchange &Exchange = *Exchange_;

  const connectivity_r &Receivers = *Receivers_;

  long long NumReceivers;
  GetConnectivityReceiverSideCount(Receivers, NumReceivers);

  std::vector<long long> NextBufferEntry(NumRecvs_, 0);

  MPI_Waitall(NumRecvs_, MPIRequests_.data(), MPI_STATUSES_IGNORE);

  for (long long iReceiverOrder = 0; iReceiverOrder < NumReceivers; ++iReceiverOrder) {
    long long iReceiver = Exchange.ReceiversSorted_[iReceiverOrder];
    int iRecv = Exchange.ReceiverRecvIndices_[iReceiver];
    if (iRecv >= 0) {
      const exchange::recv &Recv = Exchange.Recvs_[iRecv];
      long long iBuffer = NextBufferEntry[iRecv];
      for (int iCount = 0; iCount < Count_; ++iCount) {
        ReceiverValues_[iCount][iReceiver] = user_value_type(Buffers_[iRecv][iCount*Recv.Count+
          iBuffer]);
      }
      ++NextBufferEntry[iRecv];
    }
  }

  Buffers_.clear();
  MPIRequests_.clear();
  ReceiverValues_.clear();

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

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {
    Disperse_->Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);
  }

  void Disperse(const void * const *ReceiverValues, void **GridValues) {
    Disperse_->Disperse(ReceiverValues, GridValues);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange,
      array_layout GridValuesLayout) = 0;
    virtual void Disperse(const void * const *ReceiverValues, void **GridValues) = 0;
  };

  template <typename T> class model : public concept {
  public:
    using user_value_type = typename T::user_value_type;
    explicit model(T Disperse):
      Disperse_(std::move(Disperse))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange,
      array_layout GridValuesLayout) override {
      Disperse_.Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);
      ReceiverValues_.resize(Count);
      GridValues_.resize(Count);
    }
    virtual void Disperse(const void * const *ReceiverValuesVoid, void **GridValuesVoid) override {
      OVK_DEBUG_ASSERT(ReceiverValuesVoid || ReceiverValues_.size() == 0, "Invalid receiver values "
        "pointer.");
      OVK_DEBUG_ASSERT(GridValuesVoid || GridValues_.size() == 0, "Invalid grid values pointer.");
      for (int iCount = 0; iCount < int(ReceiverValues_.size()); ++iCount) {
        ReceiverValues_[iCount] = static_cast<const user_value_type *>(ReceiverValuesVoid[iCount]);
      }
      for (int iCount = 0; iCount < int(GridValues_.size()); ++iCount) {
        GridValues_[iCount] = static_cast<user_value_type *>(GridValuesVoid[iCount]);
      }
      Disperse_.Disperse(ReceiverValues_.data(), GridValues_.data());
    }
    T Disperse_;
    std::vector<const user_value_type *> ReceiverValues_;
    std::vector<user_value_type *> GridValues_;
  };

  std::unique_ptr<concept> Disperse_;

};

template <typename T> class disperse_base {

public:

  using user_value_type = T;
  using value_type = no_bool<T>;

  disperse_base() = default;
  disperse_base(const disperse_base &) = delete;
  disperse_base(disperse_base &&) = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {

    Count_ = Count;
    GridValuesRange_ = GridValuesRange;
    GridValuesLayout_ = GridValuesLayout;

    const connectivity &Connectivity = *Exchange.Connectivity_;

    GetConnectivityReceiverSide(Connectivity, Receivers_);
    const connectivity_r &Receivers = *Receivers_;

    GetConnectivityReceiverSideGrid(Receivers, Grid_);
    const grid &Grid = *Grid_;

    if (OVK_DEBUG) {
      range LocalRange;
      GetGridLocalRange(Grid, LocalRange);
      OVK_DEBUG_ASSERT(RangeIncludes(GridValuesRange, LocalRange), "Invalid grid values range.");
    }

  }

protected:

  const connectivity_r *Receivers_;
  const grid *Grid_;
  int Count_;
  range GridValuesRange_;
  array_layout GridValuesLayout_;

  void CheckValuesPointers(const user_value_type * const *ReceiverValues, const user_value_type *
    const *GridValues) {

    const connectivity_r &Receivers = *Receivers_;
    long long NumReceivers;
    GetConnectivityReceiverSideCount(Receivers, NumReceivers);

    if (NumReceivers > 0) {
      for (int iCount = 0; iCount < Count_; ++iCount) {
        OVK_DEBUG_ASSERT(ReceiverValues[iCount], "Invalid receiver values pointer.");
      }
      for (int iCount = 0; iCount < Count_; ++iCount) {
        OVK_DEBUG_ASSERT(GridValues[iCount], "Invalid grid values pointer.");
      }
    }

  }

};

template <typename T> class disperse_overwrite : public disperse_base<T> {

protected:

  using parent_type = disperse_base<T>;

  using parent_type::Receivers_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesLayout_;

public:

  using typename parent_type::user_value_type;
  using typename parent_type::value_type;

  disperse_overwrite() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {

    parent_type::Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);

  }

  void Disperse(const user_value_type * const *ReceiverValues, user_value_type **GridValues) {

    const connectivity_r &Receivers = *Receivers_;

    long long NumReceivers;
    GetConnectivityReceiverSideCount(Receivers, NumReceivers);

    if (OVK_DEBUG) parent_type::CheckValuesPointers(ReceiverValues, GridValues);

    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      int Point[MAX_DIMS] = {
        Receivers.Points_[0][iReceiver],
        Receivers.Points_[1][iReceiver],
        Receivers.Points_[2][iReceiver]
      };
      long long iPoint = RangeTupleToIndex(GridValuesRange_, GridValuesLayout_, Point);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        GridValues[iCount][iPoint] = ReceiverValues[iCount][iReceiver];
      }
    }

  }

};

template <typename T> class disperse_append : public disperse_base<T> {

protected:

  using parent_type = disperse_base<T>;

  using parent_type::Receivers_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesLayout_;

public:

  using typename parent_type::user_value_type;
  using typename parent_type::value_type;

  disperse_append() = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange, array_layout
    GridValuesLayout) {

    parent_type::Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);

  }

  void Disperse(const user_value_type * const *ReceiverValues, user_value_type **GridValues) {

    const connectivity_r &Receivers = *Receivers_;

    long long NumReceivers;
    GetConnectivityReceiverSideCount(Receivers, NumReceivers);

    if (OVK_DEBUG) parent_type::CheckValuesPointers(ReceiverValues, GridValues);

    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      int Point[MAX_DIMS] = {
        Receivers.Points_[0][iReceiver],
        Receivers.Points_[1][iReceiver],
        Receivers.Points_[2][iReceiver]
      };
      long long iPoint = RangeTupleToIndex(GridValuesRange_, GridValuesLayout_, Point);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        GridValues[iCount][iPoint] += ReceiverValues[iCount][iReceiver];
      }
    }

  }

};

}

namespace core {

void Disperse(const exchange &Exchange, data_type ValueType, int Count, disperse_op DisperseOp,
  const void * const *ReceiverValues, const range &GridValuesRange, array_layout GridValuesLayout,
  void **GridValues) {

  disperse DisperseWrapper;

  switch (DisperseOp) {
  case disperse_op::OVERWRITE:
    switch (ValueType) {
    case data_type::BOOL: DisperseWrapper = disperse_overwrite<bool>(); break;
    case data_type::BYTE: DisperseWrapper = disperse_overwrite<unsigned char>(); break;
    case data_type::INT: DisperseWrapper = disperse_overwrite<int>(); break;
    case data_type::LONG: DisperseWrapper = disperse_overwrite<long>(); break;
    case data_type::LONG_LONG: DisperseWrapper = disperse_overwrite<long long>(); break;
    case data_type::UNSIGNED_INT: DisperseWrapper = disperse_overwrite<unsigned int>(); break;
    case data_type::UNSIGNED_LONG: DisperseWrapper = disperse_overwrite<unsigned long>(); break;
    case data_type::UNSIGNED_LONG_LONG: DisperseWrapper = disperse_overwrite<unsigned long long>(); break;
    case data_type::FLOAT: DisperseWrapper = disperse_overwrite<float>(); break;
    case data_type::DOUBLE: DisperseWrapper = disperse_overwrite<double>(); break;
    }
    break;
  case disperse_op::APPEND:
    switch (ValueType) {
    case data_type::BYTE: DisperseWrapper = disperse_append<unsigned char>(); break;
    case data_type::INT: DisperseWrapper = disperse_append<int>(); break;
    case data_type::LONG: DisperseWrapper = disperse_append<long>(); break;
    case data_type::LONG_LONG: DisperseWrapper = disperse_append<long long>(); break;
    case data_type::UNSIGNED_INT: DisperseWrapper = disperse_append<unsigned int>(); break;
    case data_type::UNSIGNED_LONG: DisperseWrapper = disperse_append<unsigned long>(); break;
    case data_type::UNSIGNED_LONG_LONG: DisperseWrapper = disperse_append<unsigned long long>(); break;
    case data_type::FLOAT: DisperseWrapper = disperse_append<float>(); break;
    case data_type::DOUBLE: DisperseWrapper = disperse_append<double>(); break;
    default:
      OVK_DEBUG_ASSERT(false, "Invalid data type for append disperse operation.");
      break;
    }
    break;
  }

  DisperseWrapper.Initialize(Exchange, Count, GridValuesRange, GridValuesLayout);
  DisperseWrapper.Disperse(ReceiverValues, GridValues);

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
