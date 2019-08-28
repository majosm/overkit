// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/PartitionPool.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Set.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

partition_pool::partition_pool(std::shared_ptr<context> Context, comm_view Comm, set<int>
  NeighborRanks):
  Context_(std::move(Context)),
  Comm_(Comm),
  NeighborRanks_(std::move(NeighborRanks)),
  Partitions_(&CartLess_)
{}

const std::shared_ptr<const partition> &partition_pool::Insert(std::shared_ptr<const partition>
  Partition) {

  auto &GloballyMatching = Partitions_.Fetch(Partition->Cart());

  auto Iter = GloballyMatching.Begin();

  while (Iter != GloballyMatching.End()) {
    int Matches =
      (*Iter)->LocalRange() == Partition->LocalRange() &&
      (*Iter)->ExtendedRange() == Partition->ExtendedRange() &&
      (*Iter)->SubregionCount() == Partition->SubregionCount();
    MPI_Allreduce(MPI_IN_PLACE, &Matches, 1, MPI_INT, MPI_LAND, Comm_);
    if (Matches) {
      *Iter = std::move(Partition);
      break;
    }
    ++Iter;
  }

  if (Iter == GloballyMatching.End()) {
    Iter = GloballyMatching.Insert(Iter, std::move(Partition));
  }

  return *Iter;

}

const std::shared_ptr<const partition> &partition_pool::Fetch(const cart &Cart, const range
  &LocalRange, const range &ExtendedRange, int NumSubregions) const {

  auto &GloballyMatching = Partitions_.Fetch(Cart);

  auto Iter = GloballyMatching.Begin();

  while (Iter != GloballyMatching.End()) {
    int Matches =
      (*Iter)->LocalRange() == LocalRange &&
      (*Iter)->ExtendedRange() == ExtendedRange &&
      (*Iter)->SubregionCount() == NumSubregions;
    MPI_Allreduce(MPI_IN_PLACE, &Matches, 1, MPI_INT, MPI_LAND, Comm_);
    if (Matches) break;
    ++Iter;
  }

  if (Iter == GloballyMatching.End()) {
    Iter = GloballyMatching.Insert(Iter, std::make_shared<partition>(Context_, Cart, Comm_,
      LocalRange, ExtendedRange, NumSubregions, NeighborRanks_));
  }

  return *Iter;

}

void partition_pool::Erase(const cart &Cart, const range &LocalRange, const range &ExtendedRange,
  int NumSubregions) {

  auto CartMapIter = Partitions_.Find(Cart);

  if (CartMapIter == Partitions_.End()) return;

  auto &GloballyMatching = CartMapIter->Value();

  auto Iter = GloballyMatching.Begin();

  while (Iter != GloballyMatching.End()) {
    int Matches =
      (*Iter)->LocalRange() == LocalRange &&
      (*Iter)->ExtendedRange() == ExtendedRange &&
      (*Iter)->SubregionCount() == NumSubregions;
    MPI_Allreduce(MPI_IN_PLACE, &Matches, 1, MPI_INT, MPI_LAND, Comm_);
    if (Matches) {
      GloballyMatching.Erase(Iter);
      break;
    }
    ++Iter;
  }

  if (GloballyMatching.Empty()) {
    Partitions_.Erase(CartMapIter);
  }

}

bool partition_pool::CartLess_(const cart &Left, const cart &Right) {

  if (Left.Dimension() < Right.Dimension()) return true;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    if (Left.Range().Begin(iDim) < Right.Range().Begin(iDim)) {
      return true;
    }
  }
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    if (Left.Range().End(iDim) < Right.Range().End(iDim)) {
      return true;
    }
  }

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    if (Left.Periodic(iDim) < Right.Periodic(iDim)) {
      return true;
    }
  }

  using underlying_type = typename std::underlying_type<periodic_storage>::type;

  if (underlying_type(Left.PeriodicStorage()) < underlying_type(Right.PeriodicStorage())) {
    return true;
  }

  return false;

}

}}
