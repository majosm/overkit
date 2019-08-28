// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PARTITION_HPP_INCLUDED
#define OVK_CORE_PARTITION_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Decomp.hpp>
#include <ovk/core/DistributedRegionHash.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Halo.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Request.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/ScopeGuard.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/Tuple.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <utility>
#include <type_traits>

namespace ovk {

class partition {

public:

  using neighbor_info = core::decomp_info;

  partition(std::shared_ptr<context> Context, const cart &Cart, comm_view Comm, const range
    &LocalRange, const range &ExtendedRange, int NumSubregions, array_view<const int>
    NeighborRanks);

  partition(const partition &Other) = delete;
  partition(partition &&Other) noexcept = default;

  partition &operator=(const partition &Other) = delete;
  partition &operator=(partition &&Other) noexcept = default;

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const cart &Cart() const { return Cart_; }

  comm_view Comm() const { return Comm_; }

  const range &GlobalRange() const { return Cart_.Range(); }
  const range &LocalRange() const { return LocalRange_; }
  const range &ExtendedRange() const { return ExtendedRange_; }

  int SubregionCount() const { return NumSubregions_; }

  const array<range> &LocalSubregions() const { return LocalSubregions_; }
  const array<range> &ExtendedSubregions() const { return ExtendedSubregions_; }

  const set<int> &NeighborRanks() const { return Neighbors_.Keys(); }
  const map<int,neighbor_info> &Neighbors() const { return Neighbors_; }

  template <typename FieldType, OVK_FUNCTION_REQUIRES(core::IsField<FieldType>())> request
    Exchange(FieldType &Field) const {
    return Halo_.Exchange(Field);
  }

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_const<T>::value)> request
    Exchange(field_view<T> View) const {
    return Halo_.Exchange(View);
  }

private:

  std::shared_ptr<context> Context_;

  cart Cart_;
  comm_view Comm_;
  range LocalRange_;
  range ExtendedRange_;
  int NumSubregions_;
  array<range> LocalSubregions_;
  array<range> ExtendedSubregions_;
  map<int,neighbor_info> Neighbors_;
  core::halo Halo_;

};

namespace core {

class partition_pool {

public:

  partition_pool(std::shared_ptr<context> Context, comm_view Comm, set<int> NeighborRanks);

  partition_pool(const partition_pool &Other) = delete;
  partition_pool(partition_pool &&Other) noexcept = default;

  partition_pool &operator=(const partition_pool &Other) = delete;
  partition_pool &operator=(partition_pool &&Other) noexcept = default;

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  comm_view Comm() const { return Comm_; }

  const set<int> &NeighborRanks() const { return NeighborRanks_; }

  const std::shared_ptr<const partition> &Insert(std::shared_ptr<const partition> Partition);

  const std::shared_ptr<const partition> &Fetch(const cart &Cart, const range &LocalRange, const
    range &ExtendedRange, int NumSubregions) const;

  void Erase(const cart &Cart, const range &LocalRange, const range &ExtendedRange, int
    NumSubregions);


private:

  using cart_less = bool(*)(const cart &, const cart &);

  std::shared_ptr<context> Context_;
  comm_view Comm_;
  set<int> NeighborRanks_;

  mutable map<cart,array<std::shared_ptr<const partition>>,cart_less> Partitions_;

  static bool CartLess_(const cart &Left, const cart &Right);

};

}

}

#endif
