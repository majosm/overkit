// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RANGE_HPP_INCLUDED
#define OVK_CORE_RANGE_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/RangeBase.h>
#include <ovk/core/Requires.hpp>

#include <type_traits>

namespace ovk {

namespace range_internal {

template <typename TupleType> static constexpr bool is_compatible_tuple_type() {
  return (core::IsArray<TupleType>() && std::is_same<core::array_value_type<TupleType>, int>::value
    && core::ArrayRank<TupleType>() == 1) || std::is_same<TupleType, int *>::value;
}

}

class range : public ovk_range {

public:

  using tuple_type = elem<int,MAX_DIMS>;

  range() = default;
  explicit range(int NumDims);
  template <typename BeginTupleType, typename EndTupleType, OVK_FUNCDECL_REQUIRES(
    range_internal::is_compatible_tuple_type<BeginTupleType>() &&
    range_internal::is_compatible_tuple_type<EndTupleType>())> range(int NumDims, const
    BeginTupleType &Begin, const EndTupleType &End);

  int Dimension() const;

  tuple_type Begin() const;
  int Begin(int iDim) const;
  const int *BeginData() const;

  tuple_type End() const;
  int End(int iDim) const;
  const int *EndData() const;

  tuple_type Size() const;
  int Size(int iDim) const;

  template <typename IntegerType=long long, OVK_FUNCDECL_REQUIRES(std::is_integral<
    IntegerType>::value)> IntegerType Count() const;

  bool Empty() const;

};

inline bool operator==(const range &LeftRange, const range &RightRange);
inline bool operator!=(const range &LeftRange, const range &RightRange);

template <typename TupleType, OVK_FUNCDECL_REQUIRES(range_internal::is_compatible_tuple_type<
  TupleType>())> inline bool RangeContains(const range &Range, const TupleType &Point);

inline bool RangeIncludes(const range &Range, const range &OtherRange);

template <typename TupleType, OVK_FUNCDECL_REQUIRES(range_internal::is_compatible_tuple_type<
  TupleType>())> inline void ExtendRange(range &Range, const TupleType &Point);

inline bool RangesOverlap(const range &LeftRange, const range &RightRange);

inline range UnionRanges(const range &LeftRange, const range &RightRange);
inline range IntersectRanges(const range &LeftRange, const range &RightRange);

// && in order to bind to temporary pointers (i.e., const pointer to non-const data); not sure if
// there is a nicer way to do this
template <typename TupleType, OVK_FUNCDECL_REQUIRES(range_internal::is_compatible_tuple_type<
  typename std::decay<TupleType>::type>() && !std::is_const<typename std::decay<TupleType>::type
  >::value)> inline void ClampToRange(TupleType &&Point, const range &Range);

}

#include <ovk/core/Range.inl>

#endif
