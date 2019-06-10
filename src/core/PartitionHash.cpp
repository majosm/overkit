// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/PartitionHash.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Optional.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <cmath>
#include <map>
#include <memory>
#include <numeric>

namespace ovk {
namespace core {

namespace {

tuple<int> BinDecomp(int NumDims, const tuple<int> &GlobalSize, int MaxBins);
tuple<int> MapToUniformCell(const tuple<int> &Origin, const tuple<int> &CellSize, const tuple<int>
  &Point);

}

partition_hash::partition_hash(int NumDims, comm_view Comm):
  NumDims_(NumDims),
  Comm_(Comm),
  GlobalRange_(MakeEmptyRange(NumDims)),
  LocalRange_(MakeEmptyRange(NumDims)),
  BinSize_(MakeUniformTuple<int>(NumDims, 0, 1)),
  MaxBinPartitions_(0)
{}

partition_hash::partition_hash(int NumDims, comm_view Comm, const range &GlobalRange, const range
  &LocalRange):
  NumDims_(NumDims),
  Comm_(Comm),
  GlobalRange_(GlobalRange),
  LocalRange_(LocalRange)
{

  tuple<int> NumBins = BinDecomp(NumDims_, GlobalRange_.Size(), Comm_.Size());

  int TotalBins = 1;
  for (int iDim = 0; iDim < NumDims_; ++iDim) {
    TotalBins *= NumBins(iDim);
  }

  BinRange_ = range(NumBins);
  BinIndexer_ = partition_hash::bin_indexer(BinRange_);

  bool RankHasBin = Comm_.Rank() < TotalBins;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    BinSize_(iDim) = (GlobalRange_.Size(iDim)+NumBins(iDim)-1)/NumBins(iDim);
  }

  // Figure out which bins the local partition overlaps with
  range OverlappedBinRange = MakeEmptyRange(NumDims_);
  if (!LocalRange_.Empty()) {
    tuple<int> LowerCorner, UpperCorner;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      LowerCorner(iDim) = LocalRange_.Begin(iDim);
      UpperCorner(iDim) = LocalRange_.End(iDim)-1;
    }
    OverlappedBinRange = ExtendRange(OverlappedBinRange, MapToUniformCell(GlobalRange_.Begin(),
      BinSize_, LowerCorner));
    OverlappedBinRange = ExtendRange(OverlappedBinRange, MapToUniformCell(GlobalRange_.Begin(),
      BinSize_, UpperCorner));
  }

  int NumOverlappedBins = OverlappedBinRange.Count<int>();

  array<int> OverlappedBinIndices({NumOverlappedBins});

  int iBin = 0;
  for (int k = OverlappedBinRange.Begin(2); k < OverlappedBinRange.End(2); ++k) {
    for (int j = OverlappedBinRange.Begin(1); j < OverlappedBinRange.End(1); ++j) {
      for (int i = OverlappedBinRange.Begin(0); i < OverlappedBinRange.End(0); ++i) {
        OverlappedBinIndices(iBin) = BinIndexer_.ToIndex(i,j,k);
        ++iBin;
      }
    }
  }

  // Collect partitions into bins

  array<int> PartitionRanks = core::DynamicHandshake(Comm_, OverlappedBinIndices);

  int NumPartitions = PartitionRanks.Size(0);

  array<int,3> PartitionRangeValues({{NumPartitions,2,MAX_DIMS}});
  array<MPI_Request> Requests;

  Requests.Reserve(NumOverlappedBins + NumPartitions);

  auto Isend = [&Requests](const void *Buffer, int Count, MPI_Datatype DataType, int DestRank, int
    Tag, MPI_Comm SendComm) {
    MPI_Request &Request = Requests.Append();
    MPI_Isend(Buffer, Count, DataType, DestRank, Tag, SendComm, &Request);

  };
  auto Irecv = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int SourceRank, int Tag,
    MPI_Comm RecvComm) {
    MPI_Request &Request = Requests.Append();
    MPI_Irecv(Buffer, Count, DataType, SourceRank, Tag, RecvComm, &Request);
  };

  if (RankHasBin) {
    for (int iPartition = 0; iPartition < NumPartitions; ++iPartition) {
      Irecv(PartitionRangeValues.Data(iPartition,0,0), 2*MAX_DIMS, MPI_INT,
        PartitionRanks(iPartition), 0, Comm_);
    }
  }

  array<int,2> LocalRangeValues({{2,MAX_DIMS}});
//   static_array<int,2*MAX_DIMS,2> LocalRangeValues({{2,MAX_DIMS}});
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    LocalRangeValues(0,iDim) = LocalRange_.Begin(iDim);
    LocalRangeValues(1,iDim) = LocalRange_.End(iDim);
  }

  for (iBin = 0; iBin < NumOverlappedBins; ++iBin) {
    Isend(LocalRangeValues.Data(), 2*MAX_DIMS, MPI_INT, OverlappedBinIndices(iBin), 0, Comm_);
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();
  OverlappedBinIndices.Clear();

  if (RankHasBin) {

    int BinIndex = Comm_.Rank();

    tuple<int> BinLoc = BinIndexer_.ToTuple(BinIndex);

    range Range = MakeEmptyRange(NumDims_);
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      Range.Begin(iDim) = GlobalRange_.Begin(iDim)+BinSize_(iDim)*BinLoc(iDim);
      Range.End(iDim) = GlobalRange_.Begin(iDim)+BinSize_(iDim)*(BinLoc(iDim)+1);
    }
    Range = IntersectRanges(Range, GlobalRange_);

    Bin_ = bin(BinIndex, Range, NumPartitions);
    bin &Bin = *Bin_;

    int iPartition = 0;
    for (int Rank : PartitionRanks) {
      Bin.PartitionRanges_(iPartition) = range(PartitionRangeValues.Data(iPartition,0,0),
        PartitionRangeValues.Data(iPartition,1,0));
      Bin.PartitionRanks_(iPartition) = Rank;
      ++iPartition;
    }

  }

  MPI_Allreduce(&NumPartitions, &MaxBinPartitions_, 1, MPI_INT, MPI_MAX, Comm_);

}

void partition_hash::MapToBins(array_view<const int,2> Points, array_view<int> BinIndices) const {

  long long NumPoints = (long long)(Points.Size(1));

  for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
    tuple<int> Point = {Points(0,iPoint), Points(1,iPoint), Points(2,iPoint)};
    tuple<int> BinLoc = MapToUniformCell(GlobalRange_.Begin(), BinSize_, Point);
    BinLoc = ClampToRange(BinRange_, BinLoc);
    BinIndices(iPoint) = BinIndexer_.ToIndex(BinLoc);
  }

}

void partition_hash::RetrieveBins(std::map<int, bin> &Bins) const {

  array<int> UnretrievedBinIndices;
  for (auto &Pair : Bins) {
    bin &Bin = Pair.second;
    if (Bin.Index_ < 0) {
      int BinIndex = Pair.first;
      UnretrievedBinIndices.Append(BinIndex);
    }
  }

  int NumUnretrievedBins = UnretrievedBinIndices.Count();

  array<int> RetrieveRanks = core::DynamicHandshake(Comm_, UnretrievedBinIndices);

  int NumRetrieves = RetrieveRanks.Size(0);

  array<int> NumPartitions({NumUnretrievedBins});
  array<int,4> PartitionRangeValues({{NumUnretrievedBins,MaxBinPartitions_,2,MAX_DIMS}});
  array<int,2> PartitionRanks({{NumUnretrievedBins,MaxBinPartitions_}});

  array<MPI_Request> Requests;

  Requests.Reserve(3*NumUnretrievedBins + 3*NumRetrieves);

  auto Isend = [&Requests](const void *Buffer, int Count, MPI_Datatype DataType, int DestRank, int
    Tag, MPI_Comm Comm) {
    MPI_Request &Request = Requests.Append();
    MPI_Isend(Buffer, Count, DataType, DestRank, Tag, Comm, &Request);

  };
  auto Irecv = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int SourceRank, int Tag,
    MPI_Comm Comm) {
    MPI_Request &Request = Requests.Append();
    MPI_Irecv(Buffer, Count, DataType, SourceRank, Tag, Comm, &Request);
  };

  for (int iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    Irecv(NumPartitions.Data(iBin), 1, MPI_INT, UnretrievedBinIndices(iBin), 0, Comm_);
    Irecv(PartitionRangeValues.Data(iBin,0,0,0), 2*MAX_DIMS*MaxBinPartitions_, MPI_INT,
      UnretrievedBinIndices(iBin), 0, Comm_);
    Irecv(PartitionRanks.Data(iBin,0), MaxBinPartitions_, MPI_INT, UnretrievedBinIndices(iBin),
      0, Comm_);
  }

  if (Bin_.Present()) {
    const bin &Bin = *Bin_;
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
      Isend(&Bin.NumPartitions_, 1, MPI_INT, RetrieveRank, 0, Comm_);
      Isend(BinPartitionRangeValues.Data(), 2*MAX_DIMS*Bin.NumPartitions_, MPI_INT, RetrieveRank, 0,
        Comm_);
      Isend(Bin.PartitionRanks_.Data(), Bin.NumPartitions_, MPI_INT, RetrieveRank, 0, Comm_);
      ++iRetrieveRank;
    }
  }

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();
  RetrieveRanks.Clear();

  for (int iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    int BinIndex = UnretrievedBinIndices(iBin);
    tuple<int> BinLoc = BinIndexer_.ToTuple(BinIndex);
    range Range = MakeEmptyRange(NumDims_);
    for (int iDim = 0; iDim < NumDims_; ++iDim) {
      Range.Begin(iDim) = GlobalRange_.Begin(iDim)+BinSize_[iDim]*BinLoc[iDim];
      Range.End(iDim) = GlobalRange_.Begin(iDim)+BinSize_[iDim]*(BinLoc[iDim]+1);
    }
    bin &Bin = Bins[UnretrievedBinIndices(iBin)];
    Bin = bin(BinIndex, Range, NumPartitions(iBin));
    for (int iPartition = 0; iPartition < NumPartitions(iBin); ++iPartition) {
      Bin.PartitionRanges_(iPartition) = range(PartitionRangeValues.Data(iBin,iPartition,0,0),
        PartitionRangeValues.Data(iBin,iPartition,1,0));
      Bin.PartitionRanks_(iPartition) = PartitionRanks(iBin,iPartition);
    }
  }

}

void partition_hash::FindPartitions(const std::map<int, bin> &Bins, array_view<const int,2> Points,
  array_view<const int> BinIndices, array_view<int> PartitionRanks) const {

  long long NumPoints = (long long)(Points.Size(1));

  for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
    tuple<int> Point = {Points(0,iPoint), Points(1,iPoint), Points(2,iPoint)};
    PartitionRanks(iPoint) = -1;
    auto BinIter = Bins.find(BinIndices(iPoint));
    if (BinIter != Bins.end()) {
      const bin &Bin = BinIter->second;
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

tuple<int> BinDecomp(int NumDims, const tuple<int> &GlobalSize, int MaxBins) {

  tuple<int> NumBins = {1,1,1};

  if (NumDims == 1) {

    NumBins[0] = MaxBins;

  } else {

    tuple<double> Length;
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

    tuple<int> GlobalSizeReduced;

    int iReducedDim;

    iReducedDim = 0;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        GlobalSizeReduced[iReducedDim] = GlobalSize[iDim];
        ++iReducedDim;
      }
    }

    int MaxBinsReduced = MaxBins/NumBins[iMinLengthDim];

    tuple<int> NumBinsReduced = BinDecomp(NumDims-1, GlobalSizeReduced, MaxBinsReduced);

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

tuple<int> MapToUniformCell(const tuple<int> &Origin, const tuple<int> &CellSize, const tuple<int>
  &Point) {

  tuple<int> Cell;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    int Offset = Point[iDim] - Origin[iDim];
    // Division rounding down
    Cell[iDim] = Offset/CellSize[iDim] - (Offset % CellSize[iDim] < 0);
  }

  return Cell;

}

}

}}
