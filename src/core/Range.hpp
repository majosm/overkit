// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RANGE_HPP_INCLUDED
#define OVK_CORE_RANGE_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/RangeBase.h>

#include <array>

namespace ovk {

class range : public ovk_range {

public:

  range() = default;
  explicit range(int NumDims);
  template <typename BeginArrayType, typename EndArrayType> range(int NumDims, const BeginArrayType
    &Begin, const EndArrayType &End);
  range(const range &Other) = default;

  range &operator=(const range &Other) = default;

  int Dimension() const;

  template <int N=MAX_DIMS> std::array<int,N> Begin() const;
  int Begin(int iDim) const;
  const int *BeginData() const;
  template <int N=MAX_DIMS> std::array<int,N> End() const;
  int End(int iDim) const;
  const int *EndData() const;

  template <int N=MAX_DIMS> std::array<int,N> Size() const;
  int Size(int iDim) const;

  template <typename CountType=long long> CountType Count() const;

  bool Empty() const;

};

inline bool operator==(const range &LeftRange, const range &RightRange);
inline bool operator!=(const range &LeftRange, const range &RightRange);

template <typename IndexType=long long, typename TupleArrayType> IndexType RangeTupleToIndex(
  const range &Range, array_layout Layout, const TupleArrayType &Tuple);
template <int N=MAX_DIMS, typename IndexType=long long> std::array<int,N>
  RangeIndexToTuple(const range &Range, array_layout Layout, IndexType Index);

template <typename TupleArrayType> bool RangeContains(const range &Range, const TupleArrayType
  &Tuple);

bool RangeIncludes(const range &Range, const range &OtherRange);

template <typename TupleArrayType> void ExtendRange(range &Range, const TupleArrayType &Tuple);

inline bool RangesOverlap(const range &LeftRange, const range &RightRange);

inline range UnionRanges(const range &LeftRange, const range &RightRange);
inline range IntersectRanges(const range &LeftRange, const range &RightRange);

template <typename TupleArrayType> inline void ClampToRange(TupleArrayType &Tuple, const range
  &Range);

}

#include <ovk/core/Range.inl>

#endif
