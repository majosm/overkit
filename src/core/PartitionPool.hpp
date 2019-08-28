// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PARTITION_POOL_HPP_INCLUDED
#define OVK_CORE_PARTITION_POOL_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Partition.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
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

}}

#endif
