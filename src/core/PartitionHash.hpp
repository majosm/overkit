// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PARTITION_HASH_HPP_INCLUDED
#define OVK_CORE_PARTITION_HASH_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <map>
#include <memory>

namespace ovk {
namespace core {

struct partition_bin {
  int Index_;
  range Range_;
  int NumPartitions_;
  array<range> PartitionRanges_;
  array<int> PartitionRanks_;
  // Need this for detecting unretrieved bins
  partition_bin():
    Index_(-1)
  {}
};

struct partition_hash {

  using bin_indexer = indexer<int, int, MAX_DIMS>;

  int NumDims_;
  comm_view Comm_;
  range GlobalRange_;
  range LocalRange_;
  range BinRange_;
  bin_indexer BinIndexer_;
  tuple<int> BinSize_;
  int MaxBinPartitions_;
  bool RankHasBin_;
  std::unique_ptr<partition_bin> Bin_;
//   optional<partition_bin> Bin_;

};

void CreatePartitionHash(partition_hash &Hash, int NumDims, comm_view Comm, const range
  &GlobalRange, const range &LocalRange);
void DestroyPartitionHash(partition_hash &Hash);

void MapToPartitionBins(const partition_hash &Hash, array_view<const int,2> Points, array_view<int>
  BinIndices);

void RetrievePartitionBins(const partition_hash &Hash, std::map<int, partition_bin> &Bins);

void FindPartitions(const partition_hash &Hash, const std::map<int, partition_bin> &RetrievedBins,
  array_view<const int,2> Points, array_view<const int> BinIndices, array_view<int> PartitionRanks);

}}

#endif
