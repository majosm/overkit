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
#include <ovk/core/Optional.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <map>
#include <memory>

namespace ovk {
namespace core {

class partition_hash {

public:

  class bin {
  public:
    bin():
      Index_(-1),
      Range_(MakeEmptyRange(3)),
      NumPartitions_(0)
    {}
  private:
    int Index_;
    range Range_;
    int NumPartitions_;
    array<range> PartitionRanges_;
    array<int> PartitionRanks_;
    bin(int Index, const range &Range, int NumPartitions):
      Index_(Index),
      Range_(Range),
      NumPartitions_(NumPartitions),
      PartitionRanges_({NumPartitions}),
      PartitionRanks_({NumPartitions})
    {}
    friend class partition_hash;
  };

  partition_hash(int NumDims, comm_view Comm);
  partition_hash(int NumDims, comm_view Comm, const range &GlobalRange, const range &LocalRange);

  partition_hash(const partition_hash &Other) = delete;
  partition_hash(partition_hash &&Other) noexcept = default;

  partition_hash &operator=(const partition_hash &Other) = delete;
  partition_hash &operator=(partition_hash &&Other) noexcept = default;

  void MapToBins(array_view<const int,2> Points, array_view<int> BinIndices) const;

  void RetrieveBins(std::map<int, bin> &Bins) const;

  void FindPartitions(const std::map<int, bin> &Bins, array_view<const int,2> Points, array_view<
    const int> BinIndices, array_view<int> PartitionRanks) const;

private:

  using bin_indexer = indexer<int, int, MAX_DIMS>;

  int NumDims_;
  comm_view Comm_;
  range GlobalRange_;
  range LocalRange_;
  range BinRange_;
  bin_indexer BinIndexer_;
  tuple<int> BinSize_;
  int MaxBinPartitions_;
  optional<bin> Bin_;

};

}}

#endif
