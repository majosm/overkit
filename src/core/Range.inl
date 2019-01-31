// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline range::range(int NumDims) {

  *this = static_cast<range &&>(ovkDefaultRange(NumDims));

}

template <typename BeginArrayType, typename EndArrayType> inline range::range(int NumDims,
  const BeginArrayType &Begin, const EndArrayType &End) {

  *this = static_cast<range &&>(ovkRange(NumDims, &Begin[0], &End[0]));

}

inline int range::Dimension() const {

  return ovkRangeDimension(this);

}

template <int N> inline std::array<int,N> range::Begin() const {

  std::array<int,N> Result;

  ovkRangeBegin(this, N, Result.data());

  return Result;

}

inline int range::Begin(int iDim) const {

  return ovkRangeBeginDim(this, iDim);

}

inline const int *range::BeginData() const {

  return ovkRangeBeginData(this);

}

template <int N> inline std::array<int,N> range::End() const {

  std::array<int,N> Result;

  ovkRangeEnd(this, N, Result.data());

  return Result;

}

inline int range::End(int iDim) const {

  return ovkRangeEndDim(this, iDim);

}
inline const int *range::EndData() const {

  return ovkRangeEndData(this);

}

template <int N> inline std::array<int,N> range::Size() const {

  std::array<int,N> Result;

  ovkRangeSize(this, N, Result.data());

  return Result;

}

inline int range::Size(int iDim) const {

  return ovkRangeSizeDim(this, iDim);

}

template <typename CountType> inline CountType range::Count() const {

  return CountType(ovkRangeCount(this));

}

inline bool range::Empty() const {

  return ovkRangeIsEmpty(this);

}

inline bool operator==(const ovk::range &LeftRange, const ovk::range &RightRange) {

  return ovkRangesAreEqual(&LeftRange, &RightRange);

}

inline bool operator!=(const ovk::range &LeftRange, const ovk::range &RightRange) {

  return !ovkRangesAreEqual(&LeftRange, &RightRange);

}

template <typename IndexType, typename TupleArrayType> inline IndexType RangeTupleToIndex(const
  range &Range, array_layout Layout, const TupleArrayType &Tuple) {

  return IndexType(ovkRangeTupleToIndex(&Range, ovk_array_layout(Layout), &Tuple[0]));

}

template <int N, typename IndexType> inline std::array<int,N> RangeIndexToTuple(const range &Range,
  array_layout Layout, IndexType Index) {

  std::array<int,N> Tuple;

  ovkRangeIndexToTuple(&Range, ovk_array_layout(Layout), (long long)(Index), N, Tuple.data());

  return Tuple;

}

template <typename TupleArrayType> inline bool RangeContains(const range &Range, const
  TupleArrayType &Tuple) {

  return ovkRangeContains(&Range, &Tuple[0]);

}

inline bool RangeIncludes(const range &Range, const range &OtherRange) {

  return ovkRangeIncludes(&Range, &OtherRange);

}

template <typename TupleArrayType> inline void ExtendRange(range &Range, const TupleArrayType
  &Tuple) {

  ovkExtendRange(&Range, &Tuple[0]);

}

inline bool RangesOverlap(const range &LeftRange, const range &RightRange) {

  return ovkRangesOverlap(&LeftRange, &RightRange);

}

inline range UnionRanges(const range &LeftRange, const range &RightRange) {

  return static_cast<range &&>(ovkUnionRanges(&LeftRange, &RightRange));

}

inline range IntersectRanges(const range &LeftRange, const range &RightRange) {

  return static_cast<range &&>(ovkIntersectRanges(&LeftRange, &RightRange));

}

template <typename TupleArrayType> inline void ClampToRange(TupleArrayType &Tuple, const range
  &Range) {

  ovkClampToRange(&Tuple[0], &Range);

}

}
