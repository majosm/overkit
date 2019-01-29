// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline void DefaultRange(range &Range, int NumDims) {

  ovkDefaultRange(&Range, NumDims);

}

template <typename IntArrayType1, typename IntArrayType2> inline void SetRange(range &Range,
  int NumDims, const IntArrayType1 &Begin, const IntArrayType2 &End) {

  ovkSetRange(&Range, NumDims, &Begin[0], &End[0]);

}

template <typename IntArrayType> inline void RangeSize(const range &Range, IntArrayType &Size) {

  ovkRangeSize(&Range, &Size[0]);

}

template <typename IntegerType> inline void RangeCount(const range &Range, IntegerType &Count) {

  long long CountLongLong;
  ovkRangeCount(&Range, &CountLongLong);

  Count = IntegerType(CountLongLong);

}

inline bool RangeIsEmpty(const range &Range) {

  return ovkRangeIsEmpty(&Range);

}

template <typename IntArrayType, typename IntegerType> inline void RangeTupleToIndex(const range
  &Range, array_layout Layout, const IntArrayType &Tuple, IntegerType &Index) {

  long long IndexLongLong;
  ovkRangeTupleToIndex(&Range, ovk_array_layout(Layout), &Tuple[0], &IndexLongLong);

  Index = IntegerType(IndexLongLong);

}

template <typename IntegerType, typename IntArrayType> inline void RangeIndexToTuple(const range
  &Range, array_layout Layout, IntegerType Index, IntArrayType &Tuple) {

  long long IndexLongLong = (long long)(Index);
  ovkRangeIndexToTuple(&Range, ovk_array_layout(Layout), IndexLongLong, &Tuple[0]);

}

template <typename IntArrayType> inline bool RangeContains(const range &Range, const IntArrayType
  &Point) {

  return ovkRangeContains(&Range, &Point[0]);

}

inline bool RangeIncludes(const range &Range, const range &OtherRange) {

  return ovkRangeIncludes(&Range, &OtherRange);

}

inline bool RangeOverlaps(const range &LeftRange, const range &RightRange) {

  return ovkRangeOverlaps(&LeftRange, &RightRange);

}

inline void RangeUnion(const range &LeftRange, const range &RightRange, range &UnionRange) {

  ovkRangeUnion(&LeftRange, &RightRange, &UnionRange);

}

inline void RangeIntersect(const range &LeftRange, const range &RightRange, range &IntersectRange) {

  ovkRangeIntersect(&LeftRange, &RightRange, &IntersectRange);

}

inline void RangeClamp(const range &Range, int *Point) {

  ovkRangeClamp(&Range, Point);

}

template <typename IntArrayType> inline void RangeExtend(const range &Range, const IntArrayType
  &Point, range &ExtendRange) {

  ovkRangeExtend(&Range, &Point[0], &ExtendRange);

}

}

inline bool operator==(const ovk::range &LeftRange, const ovk::range &RightRange) {

  return ovkRangeEquals(&LeftRange, &RightRange);

}

inline bool operator!=(const ovk::range &LeftRange, const ovk::range &RightRange) {

  return !ovkRangeEquals(&LeftRange, &RightRange);

}
