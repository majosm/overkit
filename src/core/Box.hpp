// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_BOX_HPP_INCLUDED
#define OVK_CORE_BOX_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/HashableRegionTraits.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/Math.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/Tuple.hpp>

namespace ovk {

using box = interval<double,MAX_DIMS>;

inline box MakeEmptyBox(int NumDims);
inline box ExtendBox(const box &Box, const tuple<double> &Point);
inline bool BoxesOverlap(const box &LeftBox, const box &RightBox);
inline box UnionBoxes(const box &LeftBox, const box &RightBox);
inline box IntersectBoxes(const box &LeftBox, const box &RightBox);
inline tuple<double> ClampToBox(const box &Box, const tuple<double> &Point);
inline box MoveBox(const box &Box, const tuple<double> &Amount);
inline box GrowBox(const box &Box, const tuple<double> &Amount);
inline box ScaleBox(const box &Box, const tuple<double> &Factor);

namespace core {
template <> struct hashable_region_traits<box> {
  using coord_type = double;
  static constexpr hashable_region_maps_to MapsTo() { return hashable_region_maps_to::RANGE; }
  static box ComputeExtents(int, const box &Region) { return Region; }
  static range MapToBins(int NumDims, const range &BinRange, const tuple<double> &LowerCorner,
    const tuple<double> &BinSize, const box &Region) {
    tuple<int> BinLocLower = ClampToRange(BinRange, MapToUniformGridCell(NumDims, LowerCorner,
      BinSize, Region.Begin()));
    tuple<int> BinLocUpper = ClampToRange(BinRange, MapToUniformGridCell(NumDims, LowerCorner,
      BinSize, Region.End()));
    range Bins = MakeEmptyRange(NumDims);
    Bins = ExtendRange(Bins, BinLocLower);
    Bins = ExtendRange(Bins, BinLocUpper);
    return Bins;
  }
};
}

}

#include <ovk/core/Box.inl>

#endif
