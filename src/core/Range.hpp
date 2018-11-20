// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RANGE_HPP_INCLUDED
#define OVK_CORE_RANGE_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/RangeBase.h>

namespace ovk {

using range = ovk_range;

inline void DefaultRange(range &Range, int NumDims) { ovkDefaultRange(&Range, NumDims); }
inline void SetRange(range &Range, int NumDims, const int *Begin, const int *End) { ovkSetRange(&Range, NumDims, Begin, End); }
inline void RangeSize(const range &Range, int *Size) { ovkRangeSize(&Range, Size); }
inline void RangeCount(const range &Range, long long &Count) { ovkRangeCount(&Range, &Count); }
inline void RangeCount(const range &Range, int &Count) { ovkRangeCountSmall(&Range, &Count); }
inline void RangeTupleToIndex(const range &Range, array_layout Layout, const int *Tuple, long long &Index) { ovkRangeTupleToIndex(&Range, ovk_array_layout(Layout), Tuple, &Index); }
inline void RangeTupleToIndex(const range &Range, array_layout Layout, const int *Tuple, int &Index) { ovkRangeTupleToIndexSmall(&Range, ovk_array_layout(Layout), Tuple, &Index); }
inline void RangeIndexToTuple(const range &Range, array_layout Layout, long long Index, int *Tuple) { ovkRangeIndexToTuple(&Range, ovk_array_layout(Layout), Index, Tuple); }
inline bool RangeIsEmpty(const range &Range) { return ovkRangeIsEmpty(&Range); }
inline bool RangeEquals(const range &LeftRange, const range &RightRange) { return ovkRangeEquals(&LeftRange, &RightRange); }
inline bool RangeContains(const range &Range, const int *Point) { return ovkRangeContains(&Range, Point); }
inline bool RangeIncludes(const range &Range, const range &OtherRange) { return ovkRangeIncludes(&Range, &OtherRange); }
inline bool RangeOverlaps(const range &LeftRange, const range &RightRange) { return ovkRangeOverlaps(&LeftRange, &RightRange); }
inline void RangeUnion(const range &LeftRange, const range &RightRange, range &UnionRange) { ovkRangeUnion(&LeftRange, &RightRange, &UnionRange); }
inline void RangeIntersect(const range &LeftRange, const range &RightRange, range &IntersectRange) { ovkRangeIntersect(&LeftRange, &RightRange, &IntersectRange); }
inline void RangeClamp(const range &Range, int *Point) { ovkRangeClamp(&Range, Point); }
inline void RangeExtend(const range &Range, const int *Point, range &ExtendRange) { ovkRangeExtend(&Range, Point, &ExtendRange); }

}

#endif
