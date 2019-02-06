// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/PartitionHash.hpp"

#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

namespace ovk {
namespace core {

namespace {

void CreatePartitionBin(partition_bin &Bin, int NumDims, int Index, const range &Range,
  int NumPartitions);
void DestroyPartitionBin(partition_bin &Bin);

void BinDecomp(int NumDims, int *GlobalSize, int MaxBins, int *NumBins);
inline void MapToUniformCell(const int *Origin, const int *CellSize, const int *Point,
  int *Cell);

}

void CreatePartitionHash(partition_hash &Hash, int NumDims, const comm &Comm, const range
  &GlobalRange, const range &LocalRange) {

  Hash.NumDims_ = NumDims;

  Hash.Comm_ = Comm;

  Hash.GlobalRange_ = GlobalRange;
  Hash.LocalRange_ = LocalRange;

  std::array<int,MAX_DIMS> GlobalSize = GlobalRange.Size();

  int NumBins[MAX_DIMS] = {1, 1, 1};
  BinDecomp(Hash.NumDims_, GlobalSize.data(), Hash.Comm_.Size(), NumBins);

  int TotalBins = 1;
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    TotalBins *= NumBins[iDim];
  }

  Hash.BinRange_ = range(NumDims, std::array<int,MAX_DIMS>({}), NumBins);

  Hash.RankHasBin_ = Hash.Comm_.Rank() < TotalBins;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Hash.BinSize_[iDim] = (GlobalSize[iDim]+NumBins[iDim]-1)/NumBins[iDim];
  }

  // Figure out which bins the local partition overlaps with
  int OverlappedBinBegin[MAX_DIMS] = {0, 0, 0};
  int OverlappedBinEnd[MAX_DIMS] = {0, 0, 1};
  if (!Hash.LocalRange_.Empty()) {
    int LowerCorner[MAX_DIMS], UpperCorner[MAX_DIMS];
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      LowerCorner[iDim] = Hash.LocalRange_.Begin(iDim);
      UpperCorner[iDim] = Hash.LocalRange_.End(iDim)-1;
    }
    int LowerBinLoc[MAX_DIMS], UpperBinLoc[MAX_DIMS];
    MapToUniformCell(Hash.GlobalRange_.BeginData(), Hash.BinSize_, LowerCorner, LowerBinLoc);
    MapToUniformCell(Hash.GlobalRange_.BeginData(), Hash.BinSize_, UpperCorner, UpperBinLoc);
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      OverlappedBinBegin[iDim] = LowerBinLoc[iDim];
      OverlappedBinEnd[iDim] = UpperBinLoc[iDim]+1;
    }
  }

  range OverlappedBinRange(NumDims, OverlappedBinBegin, OverlappedBinEnd);

  int NumOverlappedBins = OverlappedBinRange.Count<int>();

  std::vector<int> OverlappedBinIndices(NumOverlappedBins);

  int iBin = 0;
  for (int k = OverlappedBinRange.Begin(2); k < OverlappedBinRange.End(2); ++k) {
    for (int j = OverlappedBinRange.Begin(1); j < OverlappedBinRange.End(1); ++j) {
      for (int i = OverlappedBinRange.Begin(0); i < OverlappedBinRange.End(0); ++i) {
        int Bin[MAX_DIMS] = {i, j, k};
        OverlappedBinIndices[iBin] = RangeTupleToIndex<int>(Hash.BinRange_, array_layout::ROW_MAJOR,
          Bin);
        ++iBin;
      }
    }
  }

  // Collect partitions into bins

  std::set<int> PartitionRanks;

  core::DynamicHandshake(Hash.Comm_, NumOverlappedBins, OverlappedBinIndices.data(), PartitionRanks);

  int NumPartitions = PartitionRanks.size();

  std::vector<int> PartitionRangesFlat(2*MAX_DIMS*NumPartitions);
  std::vector<MPI_Request> Requests;

  Requests.reserve(NumOverlappedBins + NumPartitions);

  auto Isend = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int DestRank, int Tag,
    MPI_Comm Comm) {
    Requests.emplace_back();
    MPI_Isend(Buffer, Count, DataType, DestRank, Tag, Comm, &Requests.back());

  };
  auto Irecv = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int SourceRank, int Tag,
    MPI_Comm Comm) {
    Requests.emplace_back();
    MPI_Irecv(Buffer, Count, DataType, SourceRank, Tag, Comm, &Requests.back());
  };

  if (Hash.RankHasBin_) {
    int *PartitionRangeFlat = PartitionRangesFlat.data();
    for (int Rank : PartitionRanks) {
      Irecv(PartitionRangeFlat, 2*MAX_DIMS, MPI_INT, Rank, 0, Hash.Comm_);
      PartitionRangeFlat += 2*MAX_DIMS;
    }
  }

  int LocalRangeFlat[2*MAX_DIMS];
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    LocalRangeFlat[iDim] = Hash.LocalRange_.Begin(iDim);
    LocalRangeFlat[MAX_DIMS+iDim] = Hash.LocalRange_.End(iDim);
  }

  for (iBin = 0; iBin < NumOverlappedBins; ++iBin) {
    Isend(LocalRangeFlat, 2*MAX_DIMS, MPI_INT, OverlappedBinIndices[iBin], 0, Hash.Comm_);
  }

  MPI_Waitall(Requests.size(), Requests.data(), MPI_STATUSES_IGNORE);

  Requests.clear();
  OverlappedBinIndices.clear();

  if (Hash.RankHasBin_) {

    int BinIndex = Hash.Comm_.Rank();

    std::array<int,MAX_DIMS> BinLoc = RangeIndexToTuple(Hash.BinRange_, array_layout::ROW_MAJOR,
      BinIndex);

    int RangeBegin[MAX_DIMS], RangeEnd[MAX_DIMS];
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      RangeBegin[iDim] = Hash.GlobalRange_.Begin(iDim)+Hash.BinSize_[iDim]*BinLoc[iDim];
      RangeEnd[iDim] = Hash.GlobalRange_.Begin(iDim)+Hash.BinSize_[iDim]*(BinLoc[iDim]+1);
    }
    range Range(NumDims, RangeBegin, RangeEnd);
    Range = IntersectRanges(Range, Hash.GlobalRange_);

    Hash.Bin_.reset(new partition_bin());

    partition_bin &Bin = *Hash.Bin_;

    CreatePartitionBin(Bin, NumDims, BinIndex, Range, NumPartitions);

    int iPartition = 0;
    for (int Rank : PartitionRanks) {
      Bin.PartitionRanges_[iPartition] = range(NumDims, PartitionRangesFlat.data()+
        2*MAX_DIMS*iPartition, PartitionRangesFlat.data()+2*MAX_DIMS*iPartition+MAX_DIMS);
      Bin.PartitionRanks_[iPartition] = Rank;
      ++iPartition;
    }

  }

  MPI_Allreduce(&NumPartitions, &Hash.MaxBinPartitions_, 1, MPI_INT, MPI_MAX, Hash.Comm_);

}

void DestroyPartitionHash(partition_hash &Hash) {

  if (Hash.RankHasBin_) {
    DestroyPartitionBin(*Hash.Bin_);
    Hash.Bin_.reset();
  }

}

void MapToPartitionBins(const partition_hash &Hash, long long NumPoints, const int * const *Points,
  int *BinIndices) {

  for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
    int Point[MAX_DIMS] = {Points[0][iPoint], Points[1][iPoint], Points[2][iPoint]};
    int BinLoc[MAX_DIMS];
    MapToUniformCell(Hash.GlobalRange_.BeginData(), Hash.BinSize_, Point, BinLoc);
    ClampToRange(BinLoc, Hash.BinRange_);
    BinIndices[iPoint] = RangeTupleToIndex(Hash.BinRange_, array_layout::ROW_MAJOR, BinLoc);
  }

}

void RetrievePartitionBins(const partition_hash &Hash, std::map<int, partition_bin> &Bins) {

  int NumDims = Hash.NumDims_;

  std::vector<int> UnretrievedBinIndices;
  for (auto &Pair : Bins) {
    partition_bin &Bin = Pair.second;
    if (Bin.Index_ < 0) {
      int BinIndex = Pair.first;
      UnretrievedBinIndices.push_back(BinIndex);
    }
  }

  int NumUnretrievedBins = UnretrievedBinIndices.size();

  std::set<int> RetrieveRanks;

  core::DynamicHandshake(Hash.Comm_, NumUnretrievedBins, UnretrievedBinIndices.data(),
    RetrieveRanks);

  int NumRetrieves = RetrieveRanks.size();

  std::vector<int> NumPartitions(NumUnretrievedBins);
  std::vector<std::vector<int>> PartitionRangesFlat(NumUnretrievedBins);
  std::vector<std::vector<int>> PartitionRanks(NumUnretrievedBins);
  for (int iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    PartitionRangesFlat[iBin].resize(2*MAX_DIMS*Hash.MaxBinPartitions_);
    PartitionRanks[iBin].resize(Hash.MaxBinPartitions_);
  }

  std::vector<MPI_Request> Requests;

  Requests.reserve(3*NumUnretrievedBins + 3*NumRetrieves);

  auto Isend = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int DestRank, int Tag,
    MPI_Comm Comm) {
    Requests.emplace_back();
    MPI_Isend(Buffer, Count, DataType, DestRank, Tag, Comm, &Requests.back());

  };
  auto Irecv = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int SourceRank, int Tag,
    MPI_Comm Comm) {
    Requests.emplace_back();
    MPI_Irecv(Buffer, Count, DataType, SourceRank, Tag, Comm, &Requests.back());
  };

  for (int iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    Irecv(NumPartitions.data()+iBin, 1, MPI_INT, UnretrievedBinIndices[iBin], 0, Hash.Comm_);
    Irecv(PartitionRangesFlat[iBin].data(), 2*MAX_DIMS*Hash.MaxBinPartitions_, MPI_INT,
      UnretrievedBinIndices[iBin], 0, Hash.Comm_);
    Irecv(PartitionRanks[iBin].data(), Hash.MaxBinPartitions_, MPI_INT, UnretrievedBinIndices[iBin],
      0, Hash.Comm_);
  }

  if (Hash.RankHasBin_) {
    partition_bin &Bin = *Hash.Bin_;
    std::vector<int> LocalBinPartitionRangesFlat;
    LocalBinPartitionRangesFlat.resize(2*MAX_DIMS*Bin.NumPartitions_);
    for (int iPartition = 0; iPartition < Bin.NumPartitions_; ++iPartition) {
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        LocalBinPartitionRangesFlat[2*MAX_DIMS*iPartition+iDim] =
          Bin.PartitionRanges_[iPartition].Begin(iDim);
        LocalBinPartitionRangesFlat[2*MAX_DIMS*iPartition+MAX_DIMS+iDim] =
          Bin.PartitionRanges_[iPartition].End(iDim);
      }
    }
    auto RankIter = RetrieveRanks.begin();
    for (int iRetrieve = 0; iRetrieve < NumRetrieves; ++iRetrieve) {
      int RetrieveRank = *RankIter;
      Isend(&Bin.NumPartitions_, 1, MPI_INT, RetrieveRank, 0, Hash.Comm_);
      Isend(LocalBinPartitionRangesFlat.data(), 2*MAX_DIMS*Bin.NumPartitions_, MPI_INT,
        RetrieveRank, 0, Hash.Comm_);
      Isend(Bin.PartitionRanks_.data(), Bin.NumPartitions_, MPI_INT, RetrieveRank, 0, Hash.Comm_);
      ++RankIter;
    }
  }

  MPI_Waitall(Requests.size(), Requests.data(), MPI_STATUSES_IGNORE);

  Requests.clear();
  RetrieveRanks.clear();

  for (int iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    int BinIndex = UnretrievedBinIndices[iBin];
    std::array<int,MAX_DIMS> BinLoc = RangeIndexToTuple(Hash.BinRange_, array_layout::ROW_MAJOR,
      BinIndex);
    int RangeBegin[MAX_DIMS], RangeEnd[MAX_DIMS];
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      RangeBegin[iDim] = Hash.GlobalRange_.Begin(iDim)+Hash.BinSize_[iDim]*BinLoc[iDim];
      RangeEnd[iDim] = Hash.GlobalRange_.Begin(iDim)+Hash.BinSize_[iDim]*(BinLoc[iDim]+1);
    }
    range Range(NumDims, RangeBegin, RangeEnd);
    partition_bin &Bin = Bins[UnretrievedBinIndices[iBin]];
    CreatePartitionBin(Bin, NumDims, BinIndex, Range, NumPartitions[iBin]);
    for (int iPartition = 0; iPartition < NumPartitions[iBin]; ++iPartition) {
      Bin.PartitionRanges_[iPartition] = range(NumDims, PartitionRangesFlat[iBin].data() +
        2*MAX_DIMS*iPartition, PartitionRangesFlat[iBin].data()+2*MAX_DIMS*iPartition+MAX_DIMS);
      Bin.PartitionRanks_[iPartition] = PartitionRanks[iBin][iPartition];
    }
  }

}

void FindPartitions(const partition_hash &Hash, const std::map<int, partition_bin> &RetrievedBins,
  long long NumPoints, const int * const *Points, const int *BinIndices, int *PartitionRanks) {

  for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
    int Point[MAX_DIMS] = {Points[0][iPoint], Points[1][iPoint], Points[2][iPoint]};
    PartitionRanks[iPoint] = -1;
    auto BinIter = RetrievedBins.find(BinIndices[iPoint]);
    if (BinIter != RetrievedBins.end()) {
      const partition_bin &Bin = BinIter->second;
      for (int iPartition = 0; iPartition < Bin.NumPartitions_; ++iPartition) {
        if (RangeContains(Bin.PartitionRanges_[iPartition], Point)) {
          PartitionRanks[iPoint] = Bin.PartitionRanks_[iPartition];
          break;
        }
      }
    }
  }

}

namespace {

void CreatePartitionBin(partition_bin &Bin, int NumDims, int Index, const range &Range,
  int NumPartitions) {

  Bin.Index_ = Index;
  Bin.Range_ = Range;
  Bin.NumPartitions_ = NumPartitions;
  Bin.PartitionRanges_.resize(NumPartitions);
  Bin.PartitionRanks_.resize(NumPartitions);

}

void DestroyPartitionBin(partition_bin &Bin) {

  Bin.PartitionRanges_.clear();
  Bin.PartitionRanks_.clear();

}

void BinDecomp(int NumDims, int *GlobalSize, int MaxBins, int *NumBins) {

  if (NumDims == 1) {

    NumBins[0] = MaxBins;

  } else {

    double Length[MAX_DIMS];
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      Length[iDim] = (double)(GlobalSize[iDim]-1);
    }

    double Volume = 1.;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      Volume *= Length[iDim];
    }

    int iMinLengthDim = 0;
    double MinLength = std::numeric_limits<double>::max();
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (Length[iDim] <= MinLength) {
        iMinLengthDim = iDim;
        MinLength = Length[iDim];
      }
    }

    double Base = std::pow((double)MaxBins/Volume, 1./((double)NumDims));

    NumBins[iMinLengthDim] = std::max((int)(Length[iMinLengthDim]*Base),1);

    int GlobalSizeReduced[MAX_DIMS];

    int iReducedDim;

    iReducedDim = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        GlobalSizeReduced[iReducedDim] = GlobalSize[iDim];
        ++iReducedDim;
      }
    }

    int MaxBinsReduced = MaxBins/NumBins[iMinLengthDim];

    int NumBinsReduced[MAX_DIMS];
    BinDecomp(NumDims-1, GlobalSizeReduced, MaxBinsReduced, NumBinsReduced);

    iReducedDim = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        NumBins[iDim] = NumBinsReduced[iReducedDim];
        ++iReducedDim;
      }
    }

  }

}

inline void MapToUniformCell(const int *Origin, const int *CellSize, const int *Point, int *Cell) {

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    int Offset = Point[iDim] - Origin[iDim];
    // Division rounding down
    Cell[iDim] = Offset/CellSize[iDim] - (Offset % CellSize[iDim] < 0);
  }

}

}

}}
