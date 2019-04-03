// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/PartitionHash.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

#include <cmath>
#include <map>
#include <memory>
#include <numeric>

namespace ovk {
namespace core {

namespace {

void CreatePartitionBin(partition_bin &Bin, int NumDims, int Index, const range &Range,
  int NumPartitions);
void DestroyPartitionBin(partition_bin &Bin);

elem<int,MAX_DIMS> BinDecomp(int NumDims, const elem<int,MAX_DIMS> &GlobalSize, int MaxBins);
inline elem<int,MAX_DIMS> MapToUniformCell(const elem<int,MAX_DIMS> &Origin, const elem<int,
  MAX_DIMS> &CellSize, const elem<int,MAX_DIMS> &Point);

}

void CreatePartitionHash(partition_hash &Hash, int NumDims, const comm &Comm, const range
  &GlobalRange, const range &LocalRange) {

  Hash.NumDims_ = NumDims;

  Hash.Comm_ = Comm;

  Hash.GlobalRange_ = GlobalRange;
  Hash.LocalRange_ = LocalRange;

  elem<int,MAX_DIMS> NumBins = BinDecomp(Hash.NumDims_, Hash.GlobalRange_.Size(),
    Hash.Comm_.Size());

  int TotalBins = 1;
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    TotalBins *= NumBins[iDim];
  }

  Hash.BinRange_ = range(MakeUniformElem<int,MAX_DIMS>(0), NumBins);
  Hash.BinIndexer_ = partition_hash::bin_indexer(Hash.BinRange_);

  Hash.RankHasBin_ = Hash.Comm_.Rank() < TotalBins;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Hash.BinSize_[iDim] = (Hash.GlobalRange_.Size(iDim)+NumBins[iDim]-1)/NumBins[iDim];
  }

  // Figure out which bins the local partition overlaps with
  range OverlappedBinRange = MakeEmptyRange(NumDims);
  if (!Hash.LocalRange_.Empty()) {
    elem<int,MAX_DIMS> LowerCorner, UpperCorner;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      LowerCorner[iDim] = Hash.LocalRange_.Begin(iDim);
      UpperCorner[iDim] = Hash.LocalRange_.End(iDim)-1;
    }
    ExtendRange(OverlappedBinRange, MapToUniformCell(Hash.GlobalRange_.Begin(), Hash.BinSize_,
      LowerCorner));
    ExtendRange(OverlappedBinRange, MapToUniformCell(Hash.GlobalRange_.Begin(), Hash.BinSize_,
      UpperCorner));
  }

  int NumOverlappedBins = OverlappedBinRange.Count<int>();

  array<int> OverlappedBinIndices({NumOverlappedBins});

  int iBin = 0;
  for (int k = OverlappedBinRange.Begin(2); k < OverlappedBinRange.End(2); ++k) {
    for (int j = OverlappedBinRange.Begin(1); j < OverlappedBinRange.End(1); ++j) {
      for (int i = OverlappedBinRange.Begin(0); i < OverlappedBinRange.End(0); ++i) {
        elem<int,MAX_DIMS> Bin = {i,j,k};
        OverlappedBinIndices(iBin) = Hash.BinIndexer_.ToIndex(Bin);
        ++iBin;
      }
    }
  }

  // Collect partitions into bins

  array<int> PartitionRanks;

  core::DynamicHandshake(Hash.Comm_, OverlappedBinIndices, PartitionRanks);

  int NumPartitions = PartitionRanks.Size(0);

  array<int,3> PartitionRangeValues({{NumPartitions,2,MAX_DIMS}});
  array<MPI_Request> Requests;

  Requests.Reserve(NumOverlappedBins + NumPartitions);

  auto Isend = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int DestRank, int Tag,
    MPI_Comm Comm) {
    MPI_Request &Request = Requests.Append();
    MPI_Isend(Buffer, Count, DataType, DestRank, Tag, Comm, &Request);

  };
  auto Irecv = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int SourceRank, int Tag,
    MPI_Comm Comm) {
    MPI_Request &Request = Requests.Append();
    MPI_Irecv(Buffer, Count, DataType, SourceRank, Tag, Comm, &Request);
  };

  if (Hash.RankHasBin_) {
    for (int iPartition = 0; iPartition < NumPartitions; ++iPartition) {
      Irecv(PartitionRangeValues.Data(iPartition,0,0), 2*MAX_DIMS, MPI_INT,
        PartitionRanks(iPartition), 0, Hash.Comm_);
    }
  }

  array<int,2> LocalRangeValues({{2,MAX_DIMS}});
//   static_array<int,2*MAX_DIMS,2> LocalRangeValues({{2,MAX_DIMS}});
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    LocalRangeValues(0,iDim) = Hash.LocalRange_.Begin(iDim);
    LocalRangeValues(1,iDim) = Hash.LocalRange_.End(iDim);
  }

  for (iBin = 0; iBin < NumOverlappedBins; ++iBin) {
    Isend(LocalRangeValues.Data(), 2*MAX_DIMS, MPI_INT, OverlappedBinIndices(iBin), 0, Hash.Comm_);
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();
  OverlappedBinIndices.Clear();

  if (Hash.RankHasBin_) {

    int BinIndex = Hash.Comm_.Rank();

    elem<int,MAX_DIMS> BinLoc = Hash.BinIndexer_.ToTuple(BinIndex);

    range Range = MakeEmptyRange(NumDims);
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      Range.Begin(iDim) = Hash.GlobalRange_.Begin(iDim)+Hash.BinSize_[iDim]*BinLoc[iDim];
      Range.End(iDim) = Hash.GlobalRange_.Begin(iDim)+Hash.BinSize_[iDim]*(BinLoc[iDim]+1);
    }
    Range = IntersectRanges(Range, Hash.GlobalRange_);

    Hash.Bin_.reset(new partition_bin());

    partition_bin &Bin = *Hash.Bin_;

    CreatePartitionBin(Bin, NumDims, BinIndex, Range, NumPartitions);

    int iPartition = 0;
    for (int Rank : PartitionRanks) {
      Bin.PartitionRanges_(iPartition) = range(PartitionRangeValues.Data(iPartition,0,0),
        PartitionRangeValues.Data(iPartition,1,0));
      Bin.PartitionRanks_(iPartition) = Rank;
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

void MapToPartitionBins(const partition_hash &Hash, array_view<const int,2> Points, array_view<int>
  BinIndices) {

  long long NumPoints = (long long)(Points.Size(1));

  for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
    elem<int,MAX_DIMS> Point = {Points(0,iPoint), Points(1,iPoint), Points(2,iPoint)};
    elem<int,MAX_DIMS> BinLoc = MapToUniformCell(Hash.GlobalRange_.Begin(), Hash.BinSize_, Point);
    ClampToRange(Hash.BinRange_, BinLoc);
    BinIndices(iPoint) = Hash.BinIndexer_.ToIndex(BinLoc);
  }

}

void RetrievePartitionBins(const partition_hash &Hash, std::map<int, partition_bin> &Bins) {

  int NumDims = Hash.NumDims_;

  array<int> UnretrievedBinIndices;
  for (auto &Pair : Bins) {
    partition_bin &Bin = Pair.second;
    if (Bin.Index_ < 0) {
      int BinIndex = Pair.first;
      UnretrievedBinIndices.Append(BinIndex);
    }
  }

  int NumUnretrievedBins = UnretrievedBinIndices.Count();

  array<int> RetrieveRanks;

  core::DynamicHandshake(Hash.Comm_, UnretrievedBinIndices, RetrieveRanks);

  int NumRetrieves = RetrieveRanks.Size(0);

  array<int> NumPartitions({NumUnretrievedBins});
  array<int,4> PartitionRangeValues({{NumUnretrievedBins,Hash.MaxBinPartitions_,2,MAX_DIMS}});
  array<int,2> PartitionRanks({{NumUnretrievedBins,Hash.MaxBinPartitions_}});

  array<MPI_Request> Requests;

  Requests.Reserve(3*NumUnretrievedBins + 3*NumRetrieves);

  auto Isend = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int DestRank, int Tag,
    MPI_Comm Comm) {
    MPI_Request &Request = Requests.Append();
    MPI_Isend(Buffer, Count, DataType, DestRank, Tag, Comm, &Request);

  };
  auto Irecv = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int SourceRank, int Tag,
    MPI_Comm Comm) {
    MPI_Request &Request = Requests.Append();
    MPI_Irecv(Buffer, Count, DataType, SourceRank, Tag, Comm, &Request);
  };

  for (int iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    Irecv(NumPartitions.Data(iBin), 1, MPI_INT, UnretrievedBinIndices(iBin), 0, Hash.Comm_);
    Irecv(PartitionRangeValues.Data(iBin,0,0,0), 2*MAX_DIMS*Hash.MaxBinPartitions_, MPI_INT,
      UnretrievedBinIndices(iBin), 0, Hash.Comm_);
    Irecv(PartitionRanks.Data(iBin,0), Hash.MaxBinPartitions_, MPI_INT, UnretrievedBinIndices(iBin),
      0, Hash.Comm_);
  }

  if (Hash.RankHasBin_) {
    partition_bin &Bin = *Hash.Bin_;
    array<int,3> BinPartitionRangeValues({{Bin.NumPartitions_,2,MAX_DIMS}});
    for (int iPartition = 0; iPartition < Bin.NumPartitions_; ++iPartition) {
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        BinPartitionRangeValues(iPartition,0,iDim) = Bin.PartitionRanges_(iPartition).Begin(iDim);
        BinPartitionRangeValues(iPartition,1,iDim) = Bin.PartitionRanges_(iPartition).End(iDim);
      }
    }
    int iRetrieveRank = 0;
    for (int iRetrieve = 0; iRetrieve < NumRetrieves; ++iRetrieve) {
      int RetrieveRank = RetrieveRanks(iRetrieveRank);
      Isend(&Bin.NumPartitions_, 1, MPI_INT, RetrieveRank, 0, Hash.Comm_);
      Isend(BinPartitionRangeValues.Data(), 2*MAX_DIMS*Bin.NumPartitions_, MPI_INT, RetrieveRank, 0,
        Hash.Comm_);
      Isend(Bin.PartitionRanks_.Data(), Bin.NumPartitions_, MPI_INT, RetrieveRank, 0, Hash.Comm_);
      ++iRetrieveRank;
    }
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();
  RetrieveRanks.Clear();

  for (int iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    int BinIndex = UnretrievedBinIndices(iBin);
    elem<int,MAX_DIMS> BinLoc = Hash.BinIndexer_.ToTuple(BinIndex);
    range Range = MakeEmptyRange(NumDims);
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      Range.Begin(iDim) = Hash.GlobalRange_.Begin(iDim)+Hash.BinSize_[iDim]*BinLoc[iDim];
      Range.End(iDim) = Hash.GlobalRange_.Begin(iDim)+Hash.BinSize_[iDim]*(BinLoc[iDim]+1);
    }
    partition_bin &Bin = Bins[UnretrievedBinIndices(iBin)];
    CreatePartitionBin(Bin, NumDims, BinIndex, Range, NumPartitions(iBin));
    for (int iPartition = 0; iPartition < NumPartitions(iBin); ++iPartition) {
      Bin.PartitionRanges_(iPartition) = range(PartitionRangeValues.Data(iBin,iPartition,0,0),
        PartitionRangeValues.Data(iBin,iPartition,1,0));
      Bin.PartitionRanks_(iPartition) = PartitionRanks(iBin,iPartition);
    }
  }

}

void FindPartitions(const partition_hash &Hash, const std::map<int, partition_bin> &RetrievedBins,
  array_view<const int,2> Points, array_view<const int> BinIndices, array_view<int> PartitionRanks) {

  long long NumPoints = (long long)(Points.Size(1));

  for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
    elem<int,MAX_DIMS> Point = {Points(0,iPoint), Points(1,iPoint), Points(2,iPoint)};
    PartitionRanks(iPoint) = -1;
    auto BinIter = RetrievedBins.find(BinIndices(iPoint));
    if (BinIter != RetrievedBins.end()) {
      const partition_bin &Bin = BinIter->second;
      for (int iPartition = 0; iPartition < Bin.NumPartitions_; ++iPartition) {
        if (Bin.PartitionRanges_(iPartition).Contains(Point)) {
          PartitionRanks(iPoint) = Bin.PartitionRanks_(iPartition);
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
  Bin.PartitionRanges_.Resize({NumPartitions});
  Bin.PartitionRanks_.Resize({NumPartitions});

}

void DestroyPartitionBin(partition_bin &Bin) {

  Bin.PartitionRanges_.Clear();
  Bin.PartitionRanks_.Clear();

}

elem<int,MAX_DIMS> BinDecomp(int NumDims, const elem<int,MAX_DIMS> &GlobalSize, int MaxBins) {

  elem<int,MAX_DIMS> NumBins = {1,1,1};

  if (NumDims == 1) {

    NumBins[0] = MaxBins;

  } else {

    elem<double,MAX_DIMS> Length;
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

    elem<int,MAX_DIMS> GlobalSizeReduced;

    int iReducedDim;

    iReducedDim = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        GlobalSizeReduced[iReducedDim] = GlobalSize[iDim];
        ++iReducedDim;
      }
    }

    int MaxBinsReduced = MaxBins/NumBins[iMinLengthDim];

    elem<int,MAX_DIMS> NumBinsReduced = BinDecomp(NumDims-1, GlobalSizeReduced, MaxBinsReduced);

    iReducedDim = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        NumBins[iDim] = NumBinsReduced[iReducedDim];
        ++iReducedDim;
      }
    }

  }

  return NumBins;

}

inline elem<int,MAX_DIMS> MapToUniformCell(const elem<int,MAX_DIMS> &Origin, const elem<int,
  MAX_DIMS> &CellSize, const elem<int,MAX_DIMS> &Point) {

  elem<int,MAX_DIMS> Cell;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    int Offset = Point[iDim] - Origin[iDim];
    // Division rounding down
    Cell[iDim] = Offset/CellSize[iDim] - (Offset % CellSize[iDim] < 0);
  }

  return Cell;

}

}

}}
