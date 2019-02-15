// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline range::range(int NumDims) {

  *this = static_cast<range &&>(ovkDefaultRange(NumDims));

}

template <typename BeginTupleType, typename EndTupleType, OVK_FUNCDEF_REQUIRES(
  range_internal::is_compatible_tuple_type<BeginTupleType>() &&
  range_internal::is_compatible_tuple_type<EndTupleType>())>
  inline range::range(int NumDims, const BeginTupleType &Begin, const EndTupleType &End) {

  *this = static_cast<range &&>(ovkRange(NumDims, &Begin[0], &End[0]));

}

inline int range::Dimension() const {

  return ovkRangeDimension(this);

}

inline range::tuple_type range::Begin() const {

  tuple_type Result;

  ovkRangeBegin(this, MAX_DIMS, Result.Data());

  return Result;

}

inline int range::Begin(int iDim) const {

  return ovkRangeBeginDim(this, iDim);

}

inline const int *range::BeginData() const {

  return ovkRangeBeginData(this);

}

inline range::tuple_type range::End() const {

  tuple_type Result;

  ovkRangeEnd(this, MAX_DIMS, Result.Data());

  return Result;

}

inline int range::End(int iDim) const {

  return ovkRangeEndDim(this, iDim);

}
inline const int *range::EndData() const {

  return ovkRangeEndData(this);

}

inline range::tuple_type range::Size() const {

  tuple_type Result;

  ovkRangeSize(this, MAX_DIMS, Result.Data());

  return Result;

}

inline int range::Size(int iDim) const {

  return ovkRangeSizeDim(this, iDim);

}

template <typename IntegerType, OVK_FUNCDEF_REQUIRES(std::is_integral<IntegerType>::value)> inline
  IntegerType range::Count() const {

  return IntegerType(ovkRangeCount(this));

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

template <typename TupleType, OVK_FUNCDEF_REQUIRES(range_internal::is_compatible_tuple_type<
  TupleType>())> inline bool RangeContains(const range &Range, const TupleType &Point) {

  return ovkRangeContains(&Range, &Point[0]);

}

inline bool RangeIncludes(const range &Range, const range &OtherRange) {

  return ovkRangeIncludes(&Range, &OtherRange);

}

template <typename TupleType, OVK_FUNCDEF_REQUIRES(range_internal::is_compatible_tuple_type<
  TupleType>())> inline void ExtendRange(range &Range, const TupleType &Point) {

  ovkExtendRange(&Range, &Point[0]);

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

template <typename TupleType, OVK_FUNCDEF_REQUIRES(range_internal::is_compatible_tuple_type<
  typename std::decay<TupleType>::type>() && !std::is_const<typename std::decay<TupleType>::type
  >::value)> inline void ClampToRange(TupleType &&Point, const range &Range) {

  ovkClampToRange(&Point[0], &Range);

}

}
