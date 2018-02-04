// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/PartitionHash.h"

#include "ovk/core/Global.h"
#include "ovk/core/MPIUtils.h"
#include "ovk/core/OrderedMap.h"
#include "ovk/core/Range.h"

#include <math.h>

static void CreatePartitionBin(t_partition_bin **Bin, int NumDims, int Index, ovk_range *Range,
  int NumPartitions);
static void DestroyPartitionBin(t_partition_bin **Bin);

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

  double Aspect[MAX_DIMS];
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Aspect[iDim] = (double)(GlobalSize[iDim]-1)/((double)GlobalSize[0]-1);
  }

  double Volume = 1.;
  for (iDim = 0; iDim < NumDims; ++iDim) {
    Volume *= Aspect[iDim];
  }
  double BaseSize = pow((double)Hash->comm_size/Volume, 1./((double)NumDims));

  int NumBinsBase[MAX_DIMS];
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    NumBinsBase[iDim] = clamp((int)(Aspect[iDim]*BaseSize), 1, GlobalSize[iDim]);
  }

  // Find largest number of bins that can be distributed to ranks (at most 1 to any given rank)
  int TotalBins = NumBinsBase[0]*NumBinsBase[1]*NumBinsBase[2];
  int NumBins[MAX_DIMS] = {NumBinsBase[0], NumBinsBase[1], NumBinsBase[2]};
  int SearchEnd[MAX_DIMS] = {1, 1, 1};
  for (iDim = 0; iDim < NumDims; ++iDim) {
    ++SearchEnd[iDim];
  }
  for (k = 0; k < SearchEnd[2]; ++k) {
    for (j = 0; j < SearchEnd[1]; ++j) {
      for (i = 0; i < SearchEnd[0]; ++i) {
        int NumBinsPerturbed[MAX_DIMS] = {
          min(NumBins[0]+i, GlobalSize[0]),
          min(NumBins[1]+j, GlobalSize[1]),
          min(NumBins[2]+k, GlobalSize[2])};
        int TotalBinsPerturbed = NumBinsPerturbed[0]*NumBinsPerturbed[1]*NumBinsPerturbed[2];
        if (TotalBinsPerturbed > TotalBins && TotalBinsPerturbed <= Hash->comm_size) {
          TotalBins = TotalBinsPerturbed;
          for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
            NumBins[iDim] = NumBinsPerturbed[iDim];
          }
        }
      }
    }
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

  char SendBuffer[1], RecvBuffer[1];

  MPI_Request *SendRequests = malloc(NumOverlappedBins*sizeof(MPI_Request));
  for (iBin = 0; iBin < NumOverlappedBins; ++iBin) {
    MPI_Issend(SendBuffer, 1, MPI_CHAR, OverlappedBinIndices[iBin], 0, Hash->comm,
      SendRequests+iBin);
  }

  t_ordered_map *Partitions;
  if (Hash->rank_has_bin) {
    OMCreate(&Partitions);
  }

  t_signal *AllSendsDoneSignal;
  CreateSignal(&AllSendsDoneSignal, Hash->comm);

  bool Done = false;
  int SendsDone = false;
  while (!Done) {
    while (true) {
      int IncomingMessage;
      MPI_Status Status;
      MPI_Iprobe(MPI_ANY_SOURCE, 0, Hash->comm, &IncomingMessage, &Status);
      if (!IncomingMessage) break;
      int PartitionRank = Status.MPI_SOURCE;
      MPI_Recv(RecvBuffer, 1, MPI_CHAR, PartitionRank, 0, Hash->comm, MPI_STATUS_IGNORE);
      OMInsert(Partitions, PartitionRank, NULL);
    }
    if (SendsDone) {
      CheckSignal(AllSendsDoneSignal, &Done);
    } else {
      MPI_Testall(NumOverlappedBins, SendRequests, &SendsDone, MPI_STATUSES_IGNORE);
      if (SendsDone) {
        StartSignal(AllSendsDoneSignal);
      }
    }
  }

  MPI_Barrier(Hash->comm);

  DestroySignal(&AllSendsDoneSignal);

  int NumPartitions = 0;
  if (Hash->rank_has_bin) {
    NumPartitions = OMSize(Partitions);
  }

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

    OMDestroy(&Partitions);

  }

  OMCreate(&Hash->retrieved_bins);
  if (Hash->rank_has_bin) {
    OMInsert(Hash->retrieved_bins, Hash->bin->index, Hash->bin);
  }

  MPI_Allreduce(&NumPartitions, &Hash->max_bin_partitions, 1, MPI_INT, MPI_MAX, Hash->comm);

}

void PRIVATE(DestroyPartitionHash)(t_partition_hash **Hash_) {

  t_partition_hash *Hash = *Hash_;

  ClearRetrievedPartitionBins(Hash);

  OMDestroy(&Hash->retrieved_bins);

  if (Hash->rank_has_bin) {
    DestroyPartitionBin(&Hash->bin);
  }

  free_null(Hash_);

}

void PRIVATE(MapToPartitionBins)(const t_partition_hash *Hash, size_t NumPoints, const int **Points,
  int *BinIndices) {

  size_t iPoint;

  for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
    int Point[MAX_DIMS] = {Points[0][iPoint], Points[1][iPoint], Points[2][iPoint]};
    int BinLoc[MAX_DIMS];
    MapToUniformCell(Hash->global_range.b, Hash->bin_size, Point, BinLoc);
    ovkRangeClamp(&Hash->bin_range, BinLoc);
    ovkRangeTupleToIndexSmall(&Hash->bin_range, OVK_ROW_MAJOR, BinLoc, BinIndices+iPoint);
  }

}

void PRIVATE(RetrievePartitionBins)(t_partition_hash *Hash, int NumBins, int *BinIndices) {

  int iDim;
  int iBin;
  int iPartition;

  int NumDims = Hash->num_dims;

  int *UnretrievedBinIndices = malloc(NumBins*sizeof(int));
  int NumUnretrievedBins = 0;
  for (iBin = 0; iBin < NumBins; ++iBin) {
    t_ordered_map_entry *Entry = OMFind(Hash->retrieved_bins, BinIndices[iBin]);
    if (Entry == OMEnd(Hash->retrieved_bins)) {
      UnretrievedBinIndices[NumUnretrievedBins] = BinIndices[iBin];
      ++NumUnretrievedBins;
    }
  }

  char SendBuffer[1], RecvBuffer[1];

  MPI_Request *HandshakeSendRequests = malloc(NumUnretrievedBins*sizeof(MPI_Request));
  for (iBin = 0; iBin < NumUnretrievedBins; ++iBin) {
    MPI_Issend(SendBuffer, 1, MPI_CHAR, UnretrievedBinIndices[iBin], 0, Hash->comm,
      HandshakeSendRequests+iBin);
  }

  t_ordered_map *RetrieveRanks;
  OMCreate(&RetrieveRanks);

  t_signal *AllSendsDoneSignal;
  CreateSignal(&AllSendsDoneSignal, Hash->comm);

  bool Done = false;
  int SendsDone = false;
  while (!Done) {
    while (true) {
      int IncomingMessage;
      MPI_Status Status;
      MPI_Iprobe(MPI_ANY_SOURCE, 0, Hash->comm, &IncomingMessage, &Status);
      if (!IncomingMessage) break;
      int RetrieveRank = Status.MPI_SOURCE;
      MPI_Recv(&RecvBuffer, 1, MPI_CHAR, RetrieveRank, 0, Hash->comm, MPI_STATUS_IGNORE);
      OMInsert(RetrieveRanks, RetrieveRank, NULL);
    }
    if (SendsDone) {
      CheckSignal(AllSendsDoneSignal, &Done);
    } else {
      MPI_Testall(NumUnretrievedBins, HandshakeSendRequests, &SendsDone, MPI_STATUSES_IGNORE);
      if (SendsDone) {
        StartSignal(AllSendsDoneSignal);
      }
    }
  }

  MPI_Barrier(Hash->comm);

  DestroySignal(&AllSendsDoneSignal);

  free(HandshakeSendRequests);

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

  t_ordered_map_entry *Entry = OMBegin(RetrieveRanks);
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
    OMInsert(Hash->retrieved_bins, UnretrievedBinIndices[iBin], Bin);
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

void PRIVATE(ClearRetrievedPartitionBins)(t_partition_hash *Hash) {

  t_ordered_map_entry *Entry = OMBegin(Hash->retrieved_bins);
  while (Entry != OMEnd(Hash->retrieved_bins)) {
    t_partition_bin *Bin = OMData(Entry);
    if (!Hash->rank_has_bin || Bin != Hash->bin) {
      DestroyPartitionBin(&Bin);
    }
    Entry = OMNext(Entry);
  }

  OMClear(Hash->retrieved_bins);

  if (Hash->rank_has_bin) {
    OMInsert(Hash->retrieved_bins, Hash->bin->index, Hash->bin);
  }

}

void PRIVATE(FindPartitions)(const t_partition_hash *Hash, size_t NumPoints, const int **Points,
  const int *BinIndices, int *PartitionRanks) {

  size_t iPoint;
  int iPartition;

  for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
    int Point[MAX_DIMS] = {Points[0][iPoint], Points[1][iPoint], Points[2][iPoint]};
    PartitionRanks[iPoint] = -1;
    const t_ordered_map_entry *Entry = OMFindC(Hash->retrieved_bins, BinIndices[iPoint]);
    if (Entry != OMEndC(Hash->retrieved_bins)) {
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

static inline void MapToUniformCell(const int *Origin, const int *CellSize, const int *Point,
  int *Cell) {

  int iDim;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    int Offset = Point[iDim] - Origin[iDim];
    // Division rounding down
    Cell[iDim] = Offset/CellSize[iDim] - (Offset % CellSize[iDim] < 0);
  }

}
