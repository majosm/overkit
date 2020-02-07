// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RANGE_HPP_INCLUDED
#define OVK_CORE_RANGE_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/HashableRegionTraits.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/Math.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/Tuple.hpp>

namespace ovk {

using range = interval<int,MAX_DIMS>;

inline range MakeEmptyRange(int NumDims);
inline range ExtendRange(const range &Range, const tuple<int> &Point);
inline bool RangesOverlap(const range &LeftRange, const range &RightRange);
inline range UnionRanges(const range &LeftRange, const range &RightRange);
inline range IntersectRanges(const range &LeftRange, const range &RightRange);
inline tuple<int> ClampToRange(const range &Range, const tuple<int> &Point);

template <typename IndexType, array_layout Layout=array_layout::ROW_MAJOR> using range_indexer =
  indexer<IndexType, int, MAX_DIMS, Layout>;
template <typename IndexType> using range_indexer_r = range_indexer<IndexType,
  array_layout::ROW_MAJOR>;
template <typename IndexType> using range_indexer_c = range_indexer<IndexType,
  array_layout::COLUMN_MAJOR>;

namespace core {
template <> struct hashable_region_traits<range> {
  using coord_type = int;
  static range ComputeExtents(int, const range &Region) { return Region; }
  template <typename IndexerType> static set<typename IndexerType::index_type> MapToBins(int
    NumDims, const range &BinRange, const IndexerType &BinIndexer, const tuple<int> &LowerCorner,
    const tuple<int> &BinSize, const range &Region) {
    using index_type = typename IndexerType::index_type;
    set<index_type> Bins;
    tuple<int> RegionLower, RegionUpper;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      RegionLower(iDim) = Region.Begin(iDim);
      RegionUpper(iDim) = Region.End(iDim)-1;
    }
    tuple<int> BinLocLower = ClampToRange(BinRange, MapToUniformGridCell(NumDims, LowerCorner,
      BinSize, RegionLower));
    tuple<int> BinLocUpper = ClampToRange(BinRange, MapToUniformGridCell(NumDims, LowerCorner,
      BinSize, RegionUpper));
    range OverlappedBinRange = MakeEmptyRange(NumDims);
    OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocLower);
    OverlappedBinRange = ExtendRange(OverlappedBinRange, BinLocUpper);
    Bins.Reserve(OverlappedBinRange.Count());
    for (int k = OverlappedBinRange.Begin(2); k < OverlappedBinRange.End(2); ++k) {
      for (int j = OverlappedBinRange.Begin(1); j < OverlappedBinRange.End(1); ++j) {
        for (int i = OverlappedBinRange.Begin(0); i < OverlappedBinRange.End(0); ++i) {
          Bins.Insert(BinIndexer.ToIndex(i,j,k));
        }
      }
    }
    return Bins;
  }
};
}

}

#include <ovk/core/Range.inl>

#endif
