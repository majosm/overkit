// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/DistributedFieldOps.hpp"

#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/DistributedField.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Partition.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/ScalarOps.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <cmath>
#include <memory>
#include <utility>

namespace ovk {
namespace core {

namespace {

void DilateErode(distributed_field<bool> &Mask, int Amount, mask_bc BoundaryCondition);

}

long long CountDistributedMask(const distributed_field<bool> &Mask) {

  long long Count = 0;

  for (int k = Mask.LocalRange().Begin(2); k < Mask.LocalRange().End(2); ++k) {
    for (int j = Mask.LocalRange().Begin(1); j < Mask.LocalRange().End(1); ++j) {
      for (int i = Mask.LocalRange().Begin(0); i < Mask.LocalRange().End(0); ++i) {
        Count += (long long)(Mask(i,j,k));
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &Count, 1, MPI_LONG_LONG, MPI_SUM, Mask.Comm());

  return Count;

}

void DetectEdge(const distributed_field<bool> &Mask, edge_type EdgeType, mask_bc BoundaryCondition,
  bool IncludeExteriorPoint, distributed_field<bool> &EdgeMask, const partition_pool
  *MaybePartitionPool) {

  const std::shared_ptr<const partition> &Partition = Mask.SharedPartition();
  const cart &Cart = Partition->Cart();
  int NumDims = Cart.Dimension();

  std::shared_ptr<const partition> EdgePartition;
  if (IncludeExteriorPoint) {
    const range &LocalRange = Mask.LocalRange();
    const range &ExtendedRange = Mask.ExtendedRange();
    cart EdgeCart = CartIncludeExteriorPoint(Cart);
    range EdgeLocalRange = RangeIncludeExteriorPoint(Cart, LocalRange);
    range EdgeExtendedRange = RangeIncludeExteriorPoint(Cart, ExtendedRange);
    if (MaybePartitionPool) {
      const partition_pool &PartitionPool = *MaybePartitionPool;
      int NumSubregions = Partition->SubregionCount();
      EdgePartition = PartitionPool.Fetch(EdgeCart, EdgeLocalRange, EdgeExtendedRange,
        NumSubregions);
    } else {
      const std::shared_ptr<context> &Context = Partition->SharedContext();
      comm_view Comm = Partition->Comm();
      int NumSubregions = Partition->SubregionCount();
      const set<int> &NeighborRanks = Partition->NeighborRanks();
      EdgePartition = std::make_shared<partition>(Context, EdgeCart, Comm, EdgeLocalRange,
        EdgeExtendedRange, NumSubregions, NeighborRanks);
    }
  } else {
    EdgePartition = Partition;
  }

  EdgeMask.Assign(std::move(EdgePartition), false);

  const range &EdgeLocalRange = EdgeMask.LocalRange();

  tuple<int> GlobalLowerCorner = MakeUniformTuple<int>(NumDims, 0);
  tuple<int> GlobalUpperCorner = MakeUniformTuple<int>(NumDims, 0);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    GlobalLowerCorner(iDim) = Cart.Range().Begin(iDim);
    GlobalUpperCorner(iDim) = Cart.Range().End(iDim)-1;
  }

  auto GetMaskValue = [&](const tuple<int> &Point) -> bool {
    static constexpr tuple<int> Zero = {0,0,0};
    bool Value;
    if (Mask.ExtendedRange().Contains(Point)) {
      Value = Mask(Point);
    } else {
      switch (BoundaryCondition) {
      case mask_bc::TRUE:
        Value = true;
        break;
      case mask_bc::FALSE:
        Value = false;
        break;
      case mask_bc::MIRROR: {
        tuple<int> MirrorPoint = Point + Max(GlobalLowerCorner-Point, Zero) +
          Min(GlobalUpperCorner-Point, Zero);
        Value = Mask(MirrorPoint);
        break;
      }}
    }
    return Value;
  };

  bool EdgeValue = EdgeType == edge_type::INNER;

  for (int k = EdgeLocalRange.Begin(2); k < EdgeLocalRange.End(2); ++k) {
    for (int j = EdgeLocalRange.Begin(1); j < EdgeLocalRange.End(1); ++j) {
      for (int i = EdgeLocalRange.Begin(0); i < EdgeLocalRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        bool Value = GetMaskValue(Point);
        if (Value == EdgeValue) {
          range NeighborRange;
          for (int iDim = 0; iDim < NumDims; ++iDim) {
            NeighborRange.Begin(iDim) = Point(iDim) - 1;
            NeighborRange.End(iDim) = Point(iDim) + 2;
          }
          for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
            NeighborRange.Begin(iDim) = 0;
            NeighborRange.End(iDim) = 1;
          }
          for (int o = NeighborRange.Begin(2); o < NeighborRange.End(2); ++o) {
            for (int n = NeighborRange.Begin(1); n < NeighborRange.End(1); ++n) {
              for (int m = NeighborRange.Begin(0); m < NeighborRange.End(0); ++m) {
                tuple<int> Neighbor = {m,n,o};
                bool NeighborValue = GetMaskValue(Neighbor);
                if (NeighborValue != Value) {
                  EdgeMask(Point) = true;
                }
              }
            }
          }
        }
      }
    }
  }

  EdgeMask.Exchange();

}

void DilateMask(distributed_field<bool> &Mask, int Amount, mask_bc BoundaryCondition) {

  DilateErode(Mask, Amount, BoundaryCondition);

}

void ErodeMask(distributed_field<bool> &Mask, int Amount, mask_bc BoundaryCondition) {

  DilateErode(Mask, -Amount, BoundaryCondition);

}

namespace {

void DilateErode(distributed_field<bool> &Mask, int Amount, mask_bc BoundaryCondition) {

  if (Amount == 0) return;

  bool FillValue;
  edge_type EdgeType;
  if (Amount > 0) {
    FillValue = true;
    EdgeType = edge_type::OUTER;
  } else {
    FillValue = false;
    EdgeType = edge_type::INNER;
  }

  distributed_field<bool> EdgeMask;

  for (int iFill = 0; iFill < std::abs(Amount); ++iFill) {
    DetectEdge(Mask, EdgeType, BoundaryCondition, false, EdgeMask);
    for (long long l = 0; l < Mask.Count(); ++l) {
      if (EdgeMask[l]) {
        Mask[l] = FillValue;
      }
    }
  }

}

}

}}
