// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/PartitionHash.h"

#include "ovk/core/Global.h"
#include "ovk/core/MPIUtils.h"
#include "ovk/core/OrderedMap.h"
#include "ovk/core/Range.h"

#include <math.h>
#include <float.h>

static void CreatePartitionBin(t_partition_bin **Bin, int NumDims, int Index, ovk_range *Range,
  int NumPartitions);
static void DestroyPartitionBin(t_partition_bin **Bin);

static void BinDecomp(int NumDims, int *GlobalSize, int MaxBins, int *NumBins);
static inline void MapToUniformCell(const int *Origin, const int *CellSize, const int *Point,
  int *Cell);

void PRIVATE(CreatePartitionHash)(t_partition_hash **Hash_, int NumDims, MPI_Comm Comm,
  const ovk_range *GlobalRange, const ovk_range *LocalRange) {

  int iDim;
  int i, j, k;
  int iBin;
  int iPartition;

  *Hash_ = malloc(sizeof(t_partition_hash));
  t_partition_hash *Hash = *Hash_;

  Hash->num_dims = NumDims;

  Hash->comm = Comm;
  MPI_Comm_size(Hash->comm, &Hash->comm_size);
  MPI_Comm_rank(Hash->comm, &Hash->comm_rank);

  Hash->global_range = *GlobalRange;
  Hash->local_range = *LocalRange;

  int GlobalSize[MAX_DIMS];
  ovkRangeSize(GlobalRange, GlobalSize);

  int NumBins[MAX_DIMS] = {1, 1, 1};
  BinDecomp(Hash->num_dims, GlobalSize, Hash->comm_size, NumBins);

  int TotalBins = 1;
  for (iDim = 0; iDim < NumDims; ++iDim) {
    TotalBins *= NumBins[iDim];
  }

  ovkDefaultRange(&Hash->bin_range, NumDims);
  for (iDim = 0; iDim < NumDims; ++iDim) {
    Hash->bin_range.b[iDim] = 0;
    Hash->bin_range.e[iDim] = NumBins[iDim];
  }

  Hash->rank_has_bin = Hash->comm_rank < TotalBins;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Hash->bin_size[iDim] = (GlobalSize[iDim]+NumBins[iDim]-1)/NumBins[iDim];
  }

  // Figure out which bins the local partition overlaps with
  ovk_range OverlappedBinRange;
  ovkDefaultRange(&OverlappedBinRange, NumDims);
  if (!ovkRangeIsEmpty(&Hash->local_range)) {
    int LowerCorner[MAX_DIMS] = {0, 0, 0};
    int UpperCorner[MAX_DIMS] = {0, 0, 0};
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      LowerCorner[iDim] = Hash->local_range.b[iDim];
      UpperCorner[iDim] = Hash->local_range.e[iDim]-1;
    }
    int LowerBinLoc[MAX_DIMS], UpperBinLoc[MAX_DIMS];
    MapToUniformCell(Hash->global_range.b, Hash->bin_size, LowerCorner, LowerBinLoc);
    MapToUniformCell(Hash->global_range.b, Hash->bin_size, UpperCorner, UpperBinLoc);
    for (iDim = 0; iDim < NumDims; ++iDim) {
      OverlappedBinRange.b[iDim] = LowerBinLoc[iDim];
      OverlappedBinRange.e[iDim] = UpperBinLoc[iDim]+1;
    }
  }

  int NumOverlappedBins;
  ovkRangeCountSmall(&OverlappedBinRange, &NumOverlappedBins);

  int *OverlappedBinIndices = malloc(NumOverlappedBins*sizeof(int));

  iBin = 0;
  for (k = OverlappedBinRange.b[2]; k < OverlappedBinRange.e[2]; ++k) {
    for (j = OverlappedBinRange.b[1]; j < OverlappedBinRange.e[1]; ++j) {
      for (i = OverlappedBinRange.b[0]; i < OverlappedBinRange.e[0]; ++i) {
        int Bin[MAX_DIMS] = {i, j, k};
        ovkRangeTupleToIndexSmall(&Hash->bin_range, OVK_ROW_MAJOR, Bin, OverlappedBinIndices+iBin);
        ++iBin;
      }
    }
  }

  // Collect partitions into bins

  int LocalRangeFlat[2*MAX_DIMS];
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    LocalRangeFlat[iDim] = Hash->local_range.b[iDim];
    LocalRangeFlat[MAX_DIMS+iDim] = Hash->local_range.e[iDim];
  }

  t_ordered_map *Partitions;
  OMCreate(&Partitions);

  DynamicHandshake(Hash->comm, NumOverlappedBins, OverlappedBinIndices, Partitions);

  int NumPartitions = OMSize(Partitions);

  MPI_Request *SendRequests = malloc(NumOverlappedBins*sizeof(MPI_Request));
  MPI_Request *RecvRequests = malloc(NumPartitions*sizeof(MPI_Request));

  if (Hash->rank_has_bin) {
    t_ordered_map_entry *Entry = OMBegin(Partitions);
    iPartition = 0;
    while (Entry != OMEnd(Partitions)) {
      int PartitionRank = OMKey(Entry);
      int *PartitionRangeFlat = malloc(2*MAX_DIMS*sizeof(int));
      OMSetData(Entry, PartitionRangeFlat);
      MPI_Irecv(PartitionRangeFlat, 2*MAX_DIMS, MPI_INT, PartitionRank, 0, Hash->comm,
        RecvRequests+iPartition);
      Entry = OMNext(Entry);
      ++iPartition;
    }
  }

  for (iBin = 0; iBin < NumOverlappedBins; ++iBin) {
    MPI_Isend(LocalRangeFlat, 2*MAX_DIMS, MPI_INT, OverlappedBinIndices[iBin], 0, Hash->comm,
      SendRequests+iBin);
  }

  MPI_Waitall(NumPartitions, RecvRequests, MPI_STATUSES_IGNORE);
  MPI_Waitall(NumOverlappedBins, SendRequests, MPI_STATUSES_IGNORE);

  free(SendRequests);
  free(RecvRequests);
  free(OverlappedBinIndices);

  if (Hash->rank_has_bin) {

    int BinIndex = Hash->comm_rank;

    int BinLoc[MAX_DIMS];
    ovkRangeIndexToTuple(&Hash->bin_range, OVK_ROW_MAJOR, BinIndex, BinLoc);

    ovk_range Range;
    ovkDefaultRange(&Range, NumDims);
    for (iDim = 0; iDim < NumDims; ++iDim) {
      Range.b[iDim] = Hash->global_range.b[iDim]+Hash->bin_size[iDim]*BinLoc[iDim];
      Range.e[iDim] = Hash->global_range.b[iDim]+Hash->bin_size[iDim]*(BinLoc[iDim]+1);
    }
    ovkRangeIntersect(&Range, &Hash->global_range, &Range);

    CreatePartitionBin(&Hash->bin, NumDims, BinIndex, &Range, NumPartitions);

    t_ordered_map_entry *Entry = OMBegin(Partitions);
    iPartition = 0;
    while (Entry != OMEnd(Partitions)) {
      int *PartitionRangeFlat = OMData(Entry);
      ovkSetRange(Hash->bin->partition_ranges+iPartition, NumDims, PartitionRangeFlat,
        PartitionRangeFlat+MAX_DIMS);
      Hash->bin->partition_ranks[iPartition] = OMKey(Entry);
      free(PartitionRangeFlat);
      Entry = OMNext(Entry);
      ++iPartition;
    }

  }

  OMDestroy(&Partitions);

  MPI_Allreduce(&NumPartitions, &Hash->max_bin_partitions, 1, MPI_INT, MPI_MAX, Hash->comm);

}

void PRIVATE(DestroyPartitionHash)(t_partition_hash **Hash_) {

  t_partition_hash *Hash = *Hash_;

  if (Hash->rank_has_bin) {
    DestroyPartitionBin(&Hash->bin);
  }

  free_null(Hash_);

}

void PRIVATE(MapToPartitionBins)(const t_partition_hash *Hash, long long NumPoints, const int
  **Points, int *BinIndices) {

  long long iPoint;

  for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
    int Point[MAX_DIMS] = {Points[0][iPoint], Points[1][iPoint], Points[2][iPoint]};
    int BinLoc[MAX_DIMS];
    MapToUniformCell(Hash->global_range.b, Hash->bin_size, Point, BinLoc);
    ovkRangeClamp(&Hash->bin_range, BinLoc);
    ovkRangeTupleToIndexSmall(&Hash->bin_range, OVK_ROW_MAJOR, BinLoc, BinIndices+iPoint);
  }

}

void PRIVATE(RetrievePartitionBins)(const t_partition_hash *Hash, t_ordered_map *Bins) {

  int iDim;
  int iBin;
  int iPartition;
  t_ordered_map_entry *Entry;

  int NumDims = Hash->num_dims;

  int NumUnretrievedBins = 0;
  Entry = OMBegin(Bins);
  while (Entry != OMEnd(Bins)) {
    if (!OMData(Entry)) ++NumUnretrievedBins;
    Entry = OMNext(Entry);
  }

  int *UnretrievedBinIndices = malloc(NumUnretrievedBins*sizeof(int));
  Entry = OMBegin(Bins);
  iBin = 0;
  while (Entry != OMEnd(Bins)) {
    if (!OMData(Entry)) {
      UnretrievedBinIndices[iBin] = OMKey(Entry);
      ++iBin;
    }
    Entry = OMNext(Entry);
  }

  t_ordered_map *RetrieveRanks;
  OMCreate(&RetrieveRanks);

  DynamicHandshake(Hash->comm, NumUnretrievedBins, UnretrievedBinIndices, RetrieveRanks);

  int *NumPartitions = malloc(NumUnretrievedBins*sizeof(int));
  int **PartitionRangesFlat = malloc(NumUnretrievedBins*sizeof(int *));
  int **PartitionRanks = malloc(NumUnretrievedBins*sizeof(int *));
  for (iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    PartitionRangesFlat[iBin] = malloc(2*MAX_DIMS*Hash->max_bin_partitions*sizeof(int));
    PartitionRanks[iBin] = malloc(Hash->max_bin_partitions*sizeof(int));
  }

  MPI_Request *RecvRequests = malloc(3*NumUnretrievedBins*sizeof(MPI_Request));
  for (iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    MPI_Irecv(NumPartitions+iBin, 1, MPI_INT, UnretrievedBinIndices[iBin], 0, Hash->comm,
      RecvRequests+3*iBin);
    MPI_Irecv(PartitionRangesFlat[iBin], 2*MAX_DIMS*Hash->max_bin_partitions, MPI_INT,
      UnretrievedBinIndices[iBin], 0, Hash->comm, RecvRequests+3*iBin+1);
    MPI_Irecv(PartitionRanks[iBin], Hash->max_bin_partitions, MPI_INT,
      UnretrievedBinIndices[iBin], 0, Hash->comm, RecvRequests+3*iBin+2);
  }

  t_partition_bin *Bin = NULL;
  int *LocalBinPartitionRangesFlat = NULL;
  if (Hash->rank_has_bin) {
    Bin = Hash->bin;
    LocalBinPartitionRangesFlat = malloc(2*MAX_DIMS*Bin->num_partitions*sizeof(int));
    for (iPartition = 0; iPartition < Bin->num_partitions; ++iPartition) {
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        LocalBinPartitionRangesFlat[2*MAX_DIMS*iPartition+iDim] =
          Bin->partition_ranges[iPartition].b[iDim];
        LocalBinPartitionRangesFlat[2*MAX_DIMS*iPartition+MAX_DIMS+iDim] =
          Bin->partition_ranges[iPartition].e[iDim];
      }
    }
  }

  int NumRetrieves = OMSize(RetrieveRanks);
  MPI_Request *SendRequests = malloc(3*NumRetrieves*sizeof(MPI_Request));

  Entry = OMBegin(RetrieveRanks);
  int iRetrieve = 0;
  while (Entry != OMEnd(RetrieveRanks)) {
    int RetrieveRank = OMKey(Entry);
    MPI_Isend(&Bin->num_partitions, 1, MPI_INT, RetrieveRank, 0, Hash->comm,
      SendRequests+3*iRetrieve);
    MPI_Isend(LocalBinPartitionRangesFlat, 2*MAX_DIMS*Bin->num_partitions, MPI_INT, RetrieveRank,
      0, Hash->comm, SendRequests+3*iRetrieve+1);
    MPI_Isend(Bin->partition_ranks, Bin->num_partitions, MPI_INT, RetrieveRank, 0, Hash->comm,
      SendRequests+3*iRetrieve+2);
    Entry = OMNext(Entry);
    ++iRetrieve;
  }

  MPI_Waitall(3*NumRetrieves, SendRequests, MPI_STATUSES_IGNORE);
  MPI_Waitall(3*NumUnretrievedBins, RecvRequests, MPI_STATUSES_IGNORE);

  free(SendRequests);
  free(RecvRequests);

  OMDestroy(&RetrieveRanks);

  for (iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    int BinIndex = UnretrievedBinIndices[iBin];
    int BinLoc[MAX_DIMS];
    ovkRangeIndexToTuple(&Hash->bin_range, OVK_ROW_MAJOR, BinIndex, BinLoc);
    ovk_range Range;
    ovkDefaultRange(&Range, NumDims);
    for (iDim = 0; iDim < NumDims; ++iDim) {
      Range.b[iDim] = Hash->global_range.b[iDim]+Hash->bin_size[iDim]*BinLoc[iDim];
      Range.e[iDim] = Hash->global_range.b[iDim]+Hash->bin_size[iDim]*(BinLoc[iDim]+1);
    }
    t_partition_bin *Bin;
    CreatePartitionBin(&Bin, NumDims, BinIndex, &Range, NumPartitions[iBin]);
    for (iPartition = 0; iPartition < NumPartitions[iBin]; ++iPartition) {
      ovkSetRange(Bin->partition_ranges+iPartition, NumDims, PartitionRangesFlat[iBin] +
        2*MAX_DIMS*iPartition, PartitionRangesFlat[iBin]+2*MAX_DIMS*iPartition+MAX_DIMS);
      Bin->partition_ranks[iPartition] = PartitionRanks[iBin][iPartition];
    }
    Entry = OMFind(Bins, UnretrievedBinIndices[iBin]);
    OMSetData(Entry, Bin);
  }

  free(NumPartitions);
  for (iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    free(PartitionRangesFlat[iBin]);
    free(PartitionRanks[iBin]);
  }
  free(PartitionRangesFlat);
  free(PartitionRanks);

  if (Hash->rank_has_bin) {
    free(LocalBinPartitionRangesFlat);
  }

  free(UnretrievedBinIndices);

}

void PRIVATE(ClearPartitionBins)(t_ordered_map *Bins) {

  t_ordered_map_entry *Entry = OMBegin(Bins);
  while (Entry != OMEnd(Bins)) {
    t_partition_bin *Bin = OMData(Entry);
    DestroyPartitionBin(&Bin);
    Entry = OMNext(Entry);
  }

  OMClear(Bins);

}

void PRIVATE(FindPartitions)(const t_partition_hash *Hash, const t_ordered_map *RetrievedBins,
  long long NumPoints, const int **Points, const int *BinIndices, int *PartitionRanks) {

  long long iPoint;
  int iPartition;

  for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
    int Point[MAX_DIMS] = {Points[0][iPoint], Points[1][iPoint], Points[2][iPoint]};
    PartitionRanks[iPoint] = -1;
    const t_ordered_map_entry *Entry = OMFindC(RetrievedBins, BinIndices[iPoint]);
    if (Entry != OMEndC(RetrievedBins)) {
      const t_partition_bin *Bin = OMDataC(Entry);
      for (iPartition = 0; iPartition < Bin->num_partitions; ++iPartition) {
        if (ovkRangeContains(Bin->partition_ranges+iPartition, Point)) {
          PartitionRanks[iPoint] = Bin->partition_ranks[iPartition];
          break;
        }
      }
    }
  }

}

static void CreatePartitionBin(t_partition_bin **Bin_, int NumDims, int Index, ovk_range *Range,
  int NumPartitions) {

  *Bin_ = malloc(sizeof(t_partition_bin));
  t_partition_bin *Bin = *Bin_;

  Bin->index = Index;
  Bin->range = *Range;
  Bin->num_partitions = NumPartitions;
  Bin->partition_ranges = malloc(NumPartitions*sizeof(ovk_range));
  Bin->partition_ranks = malloc(NumPartitions*sizeof(int));

}

static void DestroyPartitionBin(t_partition_bin **Bin_) {

  t_partition_bin *Bin = *Bin_;

  free(Bin->partition_ranges);
  free(Bin->partition_ranks);

  free_null(Bin_);

}

static void BinDecomp(int NumDims, int *GlobalSize, int MaxBins, int *NumBins) {

  int iDim;
  int iReducedDim;

  if (NumDims == 1) {

    NumBins[0] = MaxBins;

  } else {

    double Length[MAX_DIMS];
    for (iDim = 0; iDim < NumDims; ++iDim) {
      Length[iDim] = (double)(GlobalSize[iDim]-1);
    }

    double Volume = 1.;
    for (iDim = 0; iDim < NumDims; ++iDim) {
      Volume *= Length[iDim];
    }

    int iMinLengthDim = 0;
    double MinLength = DBL_MAX;
    for (iDim = 0; iDim < NumDims; ++iDim) {
      if (Length[iDim] <= MinLength) {
        iMinLengthDim = iDim;
        MinLength = Length[iDim];
      }
    }

    double Base = pow((double)MaxBins/Volume, 1./((double)NumDims));

    NumBins[iMinLengthDim] = max((int)(Length[iMinLengthDim]*Base),1);

    int GlobalSizeReduced[MAX_DIMS];

    iReducedDim = 0;
    for (iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        GlobalSizeReduced[iReducedDim] = GlobalSize[iDim];
        ++iReducedDim;
      }
    }

    int MaxBinsReduced = MaxBins/NumBins[iMinLengthDim];

    int NumBinsReduced[MAX_DIMS];
    BinDecomp(NumDims-1, GlobalSizeReduced, MaxBinsReduced, NumBinsReduced);

    iReducedDim = 0;
    for (iDim = 0; iDim < NumDims; ++iDim) {
      if (iDim != iMinLengthDim) {
        NumBins[iDim] = NumBinsReduced[iReducedDim];
        ++iReducedDim;
      }
    }

  }

}

static inline void MapToUniformCell(const int *Origin, const int *CellSize, const int *Point,
  int *Cell) {

  int iDim;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    int Offset = Point[iDim] - Origin[iDim];
    // Division rounding down
    Cell[iDim] = Offset/CellSize[iDim] - (Offset % CellSize[iDim] < 0);
  }

}
