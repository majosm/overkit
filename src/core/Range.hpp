// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RANGE_HPP_INCLUDED
#define OVK_CORE_RANGE_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/RegionTraits.hpp>
#include <ovk/core/ScalarOps.hpp>
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
template <> struct region_traits<range> {
  using coord_type = int;
  static range MakeEmptyRegion(int NumDims) {
    return MakeEmptyRange(NumDims);
  }
  static range UnionRegions(const range &Left, const range &Right) {
    return UnionRanges(Left, Right);
  }
  static range IntersectRegions(const range &Left, const range &Right) {
    return IntersectRanges(Left, Right);
  }
  static tuple<int> GetRegionLowerCorner(const range &Region) {
    return Region.Begin();
  }
  static tuple<int> GetRegionUpperCorner(const range &Region) {
    tuple<int> UpperCorner;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      UpperCorner(iDim) = Region.End(iDim)-1;
    }
    return UpperCorner;
  }
};
}

}

#include <ovk/core/Range.inl>

#endif
