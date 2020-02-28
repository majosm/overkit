// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/DistributedFieldOps.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayOps.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/DistributedField.hpp"
#include "ovk/core/ElemSet.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Partition.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/ScalarOps.hpp"
#include "ovk/core/Tuple.hpp"
#include "ovk/core/UnionFind.hpp"

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
      }
      default:
        OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
        Value = false;
        break;
      }
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

void ConnectedComponents(const distributed_field<bool> &Mask, int &NumComponents,
  distributed_field<int> &ComponentLabels) {

  const std::shared_ptr<const partition> &Partition = Mask.SharedPartition();
  const cart &Cart = Partition->Cart();
  int NumDims = Cart.Dimension();
  comm_view Comm = Partition->Comm();
  const range &LocalRange = Partition->LocalRange();
  const range &ExtendedRange = Partition->ExtendedRange();
  long long NumExtended = ExtendedRange.Count();

  ComponentLabels.Assign(Partition, -1);

  int NumLocalComponents = 0;
  union_find LocalComponentSets;

  for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
    for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        bool Value = Mask(Point);
        int &Label = ComponentLabels(Point);
        range NeighborRange;
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          NeighborRange.Begin(iDim) = Point(iDim) - 1;
          NeighborRange.End(iDim) = Point(iDim) + 2;
        }
        for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
          NeighborRange.Begin(iDim) = 0;
          NeighborRange.End(iDim) = 1;
        }
        NeighborRange = IntersectRanges(NeighborRange, LocalRange);
        for (int o = NeighborRange.Begin(2); o < NeighborRange.End(2); ++o) {
          for (int n = NeighborRange.Begin(1); n < NeighborRange.End(1); ++n) {
            for (int m = NeighborRange.Begin(0); m < NeighborRange.End(0); ++m) {
              tuple<int> Neighbor = {m,n,o};
              bool NeighborValue = Mask(Neighbor);
              if (NeighborValue == Value) {
                int NeighborLabel = ComponentLabels(Neighbor);
                if (NeighborLabel >= 0) {
                  Label = NeighborLabel;
                  goto done_looping_over_neighbors;
                }
              }
            }
          }
          done_looping_over_neighbors:;
          if (Label < 0) {
            Label = NumLocalComponents;
            LocalComponentSets.Insert(Label);
            ++NumLocalComponents;
          }
        }
      }
    }
  }

  for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
    for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        bool Value = Mask(Point);
        int Label = ComponentLabels(Point);
        range NeighborRange;
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          NeighborRange.Begin(iDim) = Point(iDim) - 1;
          NeighborRange.End(iDim) = Point(iDim) + 2;
        }
        for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
          NeighborRange.Begin(iDim) = 0;
          NeighborRange.End(iDim) = 1;
        }
        NeighborRange = IntersectRanges(NeighborRange, LocalRange);
        for (int o = NeighborRange.Begin(2); o < NeighborRange.End(2); ++o) {
          for (int n = NeighborRange.Begin(1); n < NeighborRange.End(1); ++n) {
            for (int m = NeighborRange.Begin(0); m < NeighborRange.End(0); ++m) {
              tuple<int> Neighbor = {m,n,o};
              bool NeighborValue = Mask(Neighbor);
              if (NeighborValue == Value) {
                int NeighborLabel = ComponentLabels(Neighbor);
                if (NeighborLabel >= 0) {
                  LocalComponentSets.Union(Label, NeighborLabel);
                }
              }
            }
          }
        }
      }
    }
  }

  LocalComponentSets.Relabel();

  map<int,int> LocalRootToContiguous;

  for (int Label : LocalComponentSets) {
    int RootLabel = LocalComponentSets.Find(Label);
    LocalRootToContiguous.Fetch(RootLabel, -1);
  }

  for (int iLabel = 0; iLabel < LocalRootToContiguous.Count(); ++iLabel) {
    LocalRootToContiguous[iLabel].Value() = iLabel;
  }

  for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
    for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        int &Label = ComponentLabels(Point);
        int RootLabel = LocalComponentSets.Find(Label);
        Label = LocalRootToContiguous(RootLabel);
      }
    }
  }

  union_find GlobalComponentSets;

  int NumComponentsBeforeRank;
  MPI_Scan(&NumLocalComponents, &NumComponentsBeforeRank, 1, MPI_INT, MPI_SUM, Comm);
  NumComponentsBeforeRank -= NumLocalComponents;

  for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
    for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        int &Label = ComponentLabels(Point);
        Label += NumComponentsBeforeRank;
        GlobalComponentSets.Insert(Label);
      }
    }
  }

  ComponentLabels.Exchange();

  // TODO: Reimplement this part using graph contraction algorithm from Iverson, Kamath,
  // Karypsis '15

  elem_set<int,2> ExtendedUnions;

  for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
    for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        bool Value = Mask(Point);
        int Label = ComponentLabels(Point);
        range NeighborRange;
        for (int iDim = 0; iDim < NumDims; ++iDim) {
          NeighborRange.Begin(iDim) = Point(iDim) - 1;
          NeighborRange.End(iDim) = Point(iDim) + 2;
        }
        for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
          NeighborRange.Begin(iDim) = 0;
          NeighborRange.End(iDim) = 1;
        }
        NeighborRange = IntersectRanges(NeighborRange, ExtendedRange);
        for (int o = NeighborRange.Begin(2); o < NeighborRange.End(2); ++o) {
          for (int n = NeighborRange.Begin(1); n < NeighborRange.End(1); ++n) {
            for (int m = NeighborRange.Begin(0); m < NeighborRange.End(0); ++m) {
              tuple<int> Neighbor = {m,n,o};
              if (!LocalRange.Contains(Neighbor)) {
                bool NeighborValue = Mask(Neighbor);
                if (NeighborValue == Value) {
                  int NeighborLabel = ComponentLabels(Neighbor);
                  ExtendedUnions.Insert({Label,NeighborLabel});
                }
              }
            }
          }
        }
      }
    }
  }

  int NumExtendedUnions = ExtendedUnions.Count();

  array<int> NumExtendedUnionsFromRank({Comm.Size()});
  MPI_Allgather(&NumExtendedUnions, 1, MPI_INT, NumExtendedUnionsFromRank.Data(), 1, MPI_INT, Comm);

  int NumGlobalUnions = ArraySum(NumExtendedUnionsFromRank);

  array<int,2> GlobalUnions({{NumGlobalUnions,2}});

  array<int,2> ExtendedUnionsFlat({{NumExtendedUnions,2}});
  for (int iUnion = 0; iUnion < NumExtendedUnions; ++iUnion) {
    const elem<int,2> &Union = ExtendedUnions[iUnion];
    ExtendedUnionsFlat(iUnion,0) = Union(0);
    ExtendedUnionsFlat(iUnion,1) = Union(1);
  }

  array<int> GatherCounts({Comm.Size()});
  array<int> GatherOffsets({Comm.Size()});
  GatherCounts(0) = 2*NumExtendedUnionsFromRank(0);
  GatherOffsets(0) = 0;
  for (int OtherRank = 1; OtherRank < Comm.Size(); ++OtherRank) {
    GatherCounts(OtherRank) = 2*NumExtendedUnionsFromRank(OtherRank);
    GatherOffsets(OtherRank) = GatherOffsets(OtherRank-1) + GatherCounts(OtherRank-1);
  }

  MPI_Allgatherv(ExtendedUnionsFlat.Data(), 2*NumExtendedUnions, MPI_INT, GlobalUnions.Data(),
    GatherCounts.Data(), GatherOffsets.Data(), MPI_INT, Comm);

  for (int iUnion = 0; iUnion < NumGlobalUnions; ++iUnion) {
    GlobalComponentSets.Insert(GlobalUnions(iUnion,0));
    GlobalComponentSets.Insert(GlobalUnions(iUnion,1));
    GlobalComponentSets.Union(GlobalUnions(iUnion,0), GlobalUnions(iUnion,1));
  }

  GlobalComponentSets.Relabel();

  for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
    for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        int &Label = ComponentLabels(Point);
        Label = GlobalComponentSets.Find(Label);
      }
    }
  }

  ComponentLabels.Exchange();

  int MaxLabel = ArrayMax(ComponentLabels);
  MPI_Allreduce(MPI_IN_PLACE, &MaxLabel, 1, MPI_INT, MPI_MAX, Comm);

  array<int> LabelExists({MaxLabel+1}, 0);

  for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
    for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        int Label = ComponentLabels(Point);
        LabelExists(Label) = true;
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, LabelExists.Data(), LabelExists.Count(), MPI_INT, MPI_LOR, Comm);

  array<int> GlobalRootToContiguous({MaxLabel+1});

  NumComponents = 0;
  for (int Label = 0; Label <= MaxLabel; ++Label) {
    if (LabelExists(Label)) {
      GlobalRootToContiguous(Label) = NumComponents;
      ++NumComponents;
    }
  }

  for (long long l = 0; l < NumExtended; ++l) {
    int &Label = ComponentLabels[l];
    Label = GlobalRootToContiguous(Label);
  }

}

void FloodMask(distributed_field<bool> &Mask, const distributed_field<bool> &BarrierMask) {

  const std::shared_ptr<const partition> &Partition = Mask.SharedPartition();
  comm_view Comm = Partition->Comm();
  const range &ExtendedRange = Partition->ExtendedRange();

  int NumComponents;
  distributed_field<int> ComponentLabels(Partition);

  ConnectedComponents(BarrierMask, NumComponents, ComponentLabels);

  array<int> IsFloodComponent({NumComponents}, 0);

  for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
    for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
      for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        if (Mask(Point) && !BarrierMask(Point)) {
          int Label = ComponentLabels(Point);
          IsFloodComponent(Label) = 1;
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, IsFloodComponent.Data(), NumComponents, MPI_INT, MPI_MAX, Comm);

  for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
    for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
      for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
        tuple<int> Point = {i,j,k};
        if (!BarrierMask(Point)) {
          int Label = ComponentLabels(Point);
          if (IsFloodComponent(Label)) {
            Mask(Point) = true;
          }
        }
      }
    }
  }

}

}}
