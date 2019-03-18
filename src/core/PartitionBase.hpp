// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PARTITION_BASE_HPP_INCLUDED
#define OVK_CORE_PARTITION_BASE_HPP_INCLUDED

#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/PartitionHash.hpp>
#include <ovk/core/Range.hpp>

namespace ovk {
namespace core {

struct partition_info {
  int Rank;
  range LocalRange;
  range ExtendedRange;
};

range ExtendLocalRange(const cart &Cart, const range &LocalRange, int ExtendAmount);

array<int> DetectNeighbors(const cart &Cart, comm_view Comm, const range &LocalRange, const
  partition_hash &Hash);

array<partition_info> RetrievePartitionInfo(comm_view Comm, const array<int> &Ranks, const range
  &LocalRange, const range &ExtendedRange);

}}

#endif
