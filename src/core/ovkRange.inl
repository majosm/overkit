// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_RANGE_INL_INCLUDED
#define OVK_CORE_PUBLIC_RANGE_INL_INCLUDED

#include <ovkRange.h>

#include <ovkGlobal.h>

struct ovk_range {
  int nd;
  int b[OVK_MAX_DIMS];
  int e[OVK_MAX_DIMS];
};

static inline void ovkRangeDefault(ovk_range *Range, int NumDims) {

  int d;

  Range->nd = NumDims;
  for (d = 0; d < NumDims; ++d) {
    Range->b[d] = 0;
    Range->e[d] = 0;
  }
  for (d = NumDims; d < OVK_MAX_DIMS; ++d) {
    Range->b[d] = 0;
    Range->e[d] = 1;
  }

}

static inline void ovkRangeSize(const ovk_range *Range, int *Size) {

  int d;

  for (d = 0; d < Range->nd; ++d) {
    Size[d] = Range->e[d] - Range->b[d];
  }

}

static inline void ovkRangeCount(const ovk_range *Range, size_t *Count) {

  int d;

  *Count = 1;

  for (d = 0; d < Range->nd; ++d) {
    *Count *= (size_t)(Range->e[d] - Range->b[d]);
  }

}

static inline void ovkRangeTupleToIndex(const ovk_range *Range, const int *Tuple, size_t *Index) {

  int d;

  *Index = 0;
  size_t Stride = 1;
  for (d = 0; d < Range->nd; ++d) {
    *Index += Stride * (size_t)(Tuple[d] - Range->b[d]);
    Stride *= (size_t)(Range->e[d] - Range->b[d]);
  }

}

static inline void ovkRangeIndexToTuple(const ovk_range *Range, size_t Index, int *Tuple) {

  int d;

  size_t Stride = 1;
  for (d = 0; d < Range->nd; ++d) {
    int N = Range->e[d] - Range->b[d];
    Tuple[d] = Range->b[d] + (int)((Index/Stride) % N);
    Stride *= (size_t)N;
  }

}

static inline bool ovkRangeIsEmpty(const ovk_range *Range) {

  return Range->e[0] <= Range->b[0] || Range->e[1] <= Range->b[1] || Range->e[2] <= Range->b[2];

}

static inline bool ovkRangeEquals(const ovk_range *LeftRange, const ovk_range *RightRange) {

  return LeftRange->nd == RightRange->nd &&
    LeftRange->b[0] == RightRange->b[0] &&
    LeftRange->b[1] == RightRange->b[1] &&
    LeftRange->b[2] == RightRange->b[2] &&
    LeftRange->e[0] == RightRange->e[0] &&
    LeftRange->e[1] == RightRange->e[1] &&
    LeftRange->e[2] == RightRange->e[2];

}

static inline bool ovkRangeContains(const ovk_range *Range, const int *Point) {

  int d;

  for (d = 0; d < Range->nd; ++d) {
    if (Point[d] < Range->b[d] || Point[d] >= Range->e[d]) return false;
  }

  return true;

}

static inline bool ovkRangeOverlaps(const ovk_range *LeftRange, const ovk_range *RightRange) {

  return !ovkRangeIsEmpty(LeftRange) && !ovkRangeIsEmpty(RightRange) &&
    RightRange->e[0] > LeftRange->b[0] && LeftRange->e[0] > RightRange->b[0] &&
    RightRange->e[1] > LeftRange->b[1] && LeftRange->e[1] > RightRange->b[1] &&
    RightRange->e[2] > LeftRange->b[2] && LeftRange->e[2] > RightRange->b[2];

}

static inline void ovkRangeUnion(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *UnionRange) {

  int d;

  if (ovkRangeIsEmpty(LeftRange)) {
    *UnionRange = *RightRange;
  } else if (ovkRangeIsEmpty(RightRange)) {
    *UnionRange = *LeftRange;
  } else {
    UnionRange->nd = LeftRange->nd;
    for (d = 0; d < UnionRange->nd; ++d) {
      UnionRange->b[d] = ovk_min(LeftRange->b[d], RightRange->b[d]);
      UnionRange->e[d] = ovk_max(LeftRange->e[d], RightRange->e[d]);
    }
    for (d = UnionRange->nd; d < OVK_MAX_DIMS; ++d) {
      UnionRange->b[d] = 0;
      UnionRange->e[d] = 1;
    }
  }

}

static inline void ovkRangeIntersect(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *IntersectRange) {

  int d;

  IntersectRange->nd = LeftRange->nd;
  for (d = 0; d < IntersectRange->nd; ++d) {
    IntersectRange->b[d] = ovk_max(LeftRange->b[d], RightRange->b[d]);
    IntersectRange->e[d] = ovk_min(LeftRange->e[d], RightRange->e[d]);
  }
  for (d = IntersectRange->nd; d < OVK_MAX_DIMS; ++d) {
    IntersectRange->b[d] = 0;
    IntersectRange->e[d] = 1;
  }

}

static inline void ovkRangeClamp(const ovk_range *Range, const int *Point, int *ClampedPoint) {

  int d;

  for (d = 0; d < Range->nd; ++d) {
    if (Point[d] < Range->b[d]) {
      ClampedPoint[d] = Range->b[d];
    } else if (Point[d] >= Range->e[d]) {
      ClampedPoint[d] = Range->e[d]-1;
    } else {
      ClampedPoint[d] = Point[d];
    }
  }

}

#endif
