// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RANGE_HPP_INCLUDED
#define OVK_CORE_RANGE_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Interval.hpp>
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

}

#include <ovk/core/Range.inl>

#endif
