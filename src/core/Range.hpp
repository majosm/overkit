// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RANGE_HPP_INCLUDED
#define OVK_CORE_RANGE_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/RangeBase.h>

namespace ovk {

using range = ovk_range;

inline void DefaultRange(range &Range, int NumDims);
template <typename IntArrayType1, typename IntArrayType2> inline void SetRange(range &Range,
  int NumDims, const IntArrayType1 &Begin, const IntArrayType2 &End);
template <typename IntArrayType> inline void RangeSize(const range &Range, IntArrayType &Size);
template <typename IntegerType> inline void RangeCount(const range &Range, IntegerType &Count);
inline bool RangeIsEmpty(const range &Range);
template <typename IntArrayType, typename IntegerType> inline void RangeTupleToIndex(const range
  &Range, array_layout Layout, const IntArrayType &Tuple, IntegerType &Index);
template <typename IntegerType, typename IntArrayType> inline void RangeIndexToTuple(const range
  &Range, array_layout Layout, IntegerType Index, IntArrayType &Tuple);
template <typename IntArrayType> inline bool RangeContains(const range &Range, const IntArrayType
  &Point);
inline bool RangeIncludes(const range &Range, const range &OtherRange);
inline bool RangeOverlaps(const range &LeftRange, const range &RightRange);
inline void RangeUnion(const range &LeftRange, const range &RightRange, range &UnionRange);
inline void RangeIntersect(const range &LeftRange, const range &RightRange, range &IntersectRange);
inline void RangeClamp(const range &Range, int *Point);
template <typename IntArrayType> inline void RangeExtend(const range &Range, const IntArrayType
  &Point, range &ExtendRange);

}

// Can't put these in the ovk namespace (ADL won't work since ovk_range isn't defined there)
inline bool operator==(const ovk::range &LeftRange, const ovk::range &RightRange);
inline bool operator!=(const ovk::range &LeftRange, const ovk::range &RightRange);

#include <ovk/core/Range.inl>

#endif
