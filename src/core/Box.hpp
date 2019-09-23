// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_BOX_HPP_INCLUDED
#define OVK_CORE_BOX_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/RegionTraits.hpp>
#include <ovk/core/ScalarOps.hpp>
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
template <> struct region_traits<box> {
  using coord_type = double;
  static box MakeEmptyRegion(int NumDims) {
    return MakeEmptyBox(NumDims);
  }
  static box UnionRegions(const box &Left, const box &Right) {
    return UnionBoxes(Left, Right);
  }
  static box IntersectRegions(const box &Left, const box &Right) {
    return IntersectBoxes(Left, Right);
  }
  static tuple<double> GetRegionLowerCorner(const box &Region) {
    return Region.Begin();
  }
  static tuple<double> GetRegionUpperCorner(const box &Region) {
    return Region.End();
  }
};
}

}

#include <ovk/core/Box.inl>

#endif
