// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Partition.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Decomp.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Halo.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Set.hpp"

#include <mpi.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {

namespace {

array<range> CreateSubregions(int NumDims, const range &Range, int NumSubregions) {

  int iMaxStrideDim = NumDims-1;

  int Size = Range.Size(iMaxStrideDim);
  int SubregionSize = Size/NumSubregions;
  int Remainder = Size - NumSubregions*SubregionSize;

  array<range> Subregions({NumSubregions});
  for (int iSubregion = 0; iSubregion < NumSubregions; ++iSubregion) {
    Subregions(iSubregion) = Range;
    Subregions(iSubregion).Begin(iMaxStrideDim) = Range.Begin(iMaxStrideDim) +
      SubregionSize*iSubregion + Min(Remainder, iSubregion);
    Subregions(iSubregion).End(iMaxStrideDim) = Range.Begin(iMaxStrideDim) +
      SubregionSize*(iSubregion+1) + Min(Remainder, iSubregion+1);
  }

  return Subregions;

}

}

partition::partition(std::shared_ptr<context> Context, const cart &Cart, comm_view Comm, const range
  &LocalRange, const range &ExtendedRange, int NumSubregions, array_view<const int> NeighborRanks):
  Context_(std::move(Context)),
  Cart_(Cart),
  Comm_(Comm),
  LocalRange_(LocalRange),
  ExtendedRange_(ExtendedRange),
  NumSubregions_(NumSubregions),
  LocalSubregions_(CreateSubregions(Cart.Dimension(), LocalRange, NumSubregions)),
  ExtendedSubregions_(CreateSubregions(Cart.Dimension(), ExtendedRange_, NumSubregions)),
  Neighbors_(core::RetrieveDecompInfo(Comm_, NeighborRanks, LocalRange, ExtendedRange_)),
  Halo_(Context_, Cart, DuplicateComm(Comm_), LocalRange, ExtendedRange_, Neighbors_)
{}

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

}

}
