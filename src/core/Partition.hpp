// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PARTITION_HPP_INCLUDED
#define OVK_CORE_PARTITION_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Halo.hpp>
#include <ovk/core/PartitionBase.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Request.hpp>
#include <ovk/core/Tuple.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <map>
#include <memory>
#include <utility>

namespace ovk {
namespace core {

class partition {

public:

  partition(const cart &Cart, comm_view Comm, const range &LocalRange, int ExtendAmount, int
    NumSubregions, const array<int> &NeighborRanks, profiler &Profiler);

  partition(const partition &Other) = delete;
  partition(partition &&Other) noexcept = default;

  partition &operator=(const partition &Other) = delete;
  partition &operator=(partition &&Other) noexcept = default;

  const cart &Cart() const { return Cart_; }
  comm_view Comm() const { return Comm_; }
  const range &LocalRange() const { return LocalRange_; }
  const range &ExtendedRange() const { return ExtendedRange_; }
  const array<range> &LocalSubregions() const { return LocalSubregions_; }
  const array<range> &ExtendedSubregions() const { return ExtendedSubregions_; }
  const array<partition_info> &Neighbors() const { return Neighbors_; }
  const halo &Halo() const { return Halo_; }

private:

  cart Cart_;
  comm_view Comm_;
  range LocalRange_;
  range ExtendedRange_;
  array<range> LocalSubregions_;
  array<range> ExtendedSubregions_;
  array<partition_info> Neighbors_;
  halo Halo_;

};

}}

#endif
