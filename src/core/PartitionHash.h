// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PARTITION_HASH_INCLUDED
#define OVK_CORE_PARTITION_HASH_INCLUDED

#include "Global.h"
#include "OrderedMap.h"
#include "Range.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int index;
  ovk_range range;
  int num_partitions;
  ovk_range *partition_ranges;
  int *partition_ranks;
} t_partition_bin;

typedef struct {
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  ovk_range global_range;
  ovk_range local_range;
  ovk_range bin_range;
  int bin_size[MAX_DIMS];
  int max_bin_partitions;
  bool rank_has_bin;
  t_partition_bin *bin;
  t_ordered_map *retrieved_bins;
} t_partition_hash;

void PRIVATE(CreatePartitionHash)(t_partition_hash **Hash, int NumDims, MPI_Comm Comm,
  const ovk_range *GlobalRange, const ovk_range *LocalRange);
#define CreatePartitionHash(...) PRIVATE(CreatePartitionHash)(__VA_ARGS__)
void PRIVATE(DestroyPartitionHash)(t_partition_hash **Hash);
#define DestroyPartitionHash(...) PRIVATE(DestroyPartitionHash)(__VA_ARGS__)

void PRIVATE(MapToPartitionBins)(const t_partition_hash *Hash, size_t NumPoints, const int **Points,
  int *BinIndices);
#define MapToPartitionBins(...) PRIVATE(MapToPartitionBins)(__VA_ARGS__)

void PRIVATE(RetrievePartitionBins)(t_partition_hash *Hash, int NumBins, int *BinIndices);
#define RetrievePartitionBins(...) PRIVATE(RetrievePartitionBins)(__VA_ARGS__)
void PRIVATE(ClearRetrievedPartitionBins)(t_partition_hash *Hash);
#define ClearRetrievedPartitionBins(...) PRIVATE(ClearRetrievedPartitionBins)(__VA_ARGS__)

void PRIVATE(FindPartitions)(const t_partition_hash *Hash, size_t NumPoints, const int **Points,
  const int *BinIndices, int *PartitionRanks);
#define FindPartitions(...) PRIVATE(FindPartitions)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
