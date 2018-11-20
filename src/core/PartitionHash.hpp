// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PARTITION_HASH_HPP_INCLUDED
#define OVK_CORE_PARTITION_HASH_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

#include <map>
#include <memory>
#include <vector>

namespace ovk {
namespace core {

struct partition_bin {
  int Index_;
  ovk_range Range_;
  int NumPartitions_;
  std::vector<range> PartitionRanges_;
  std::vector<int> PartitionRanks_;
  // Need this for detecting unretrieved bins
  partition_bin():
    Index_(-1)
  {}
};

struct partition_hash {
  int NumDims_;
  MPI_Comm Comm_;
  int CommSize_;
  int CommRank_;
  range GlobalRange_;
  range LocalRange_;
  range BinRange_;
  int BinSize_[MAX_DIMS];
  int MaxBinPartitions_;
  bool RankHasBin_;
  std::unique_ptr<partition_bin> Bin_;
//   optional<partition_bin> Bin_;
};

void CreatePartitionHash(partition_hash &Hash, int NumDims, MPI_Comm Comm, const range &GlobalRange,
  const range &LocalRange);
void DestroyPartitionHash(partition_hash &Hash);

void MapToPartitionBins(const partition_hash &Hash, long long NumPoints, const int * const *Points,
  int *BinIndices);

void RetrievePartitionBins(const partition_hash &Hash, std::map<int, partition_bin> &Bins);

void FindPartitions(const partition_hash &Hash, const std::map<int, partition_bin> &RetrievedBins,
  long long NumPoints, const int * const *Points, const int *BinIndices, int *PartitionRanks);

}}

#endif
