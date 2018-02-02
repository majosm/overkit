// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_RANGE_INL_INCLUDED
#define OVK_CORE_PUBLIC_RANGE_INL_INCLUDED

#include <ovkRange.h>

#include <ovkGlobal.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline void ovkDefaultRange(ovk_range *Range, int NumDims) {

  int iDim;

  Range->nd = NumDims;
  for (iDim = 0; iDim < Range->nd; ++iDim) {
    Range->b[iDim] = 0;
    Range->e[iDim] = 0;
  }
  for (iDim = Range->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    Range->b[iDim] = 0;
    Range->e[iDim] = 1;
  }

}

static inline void ovkSetRange(ovk_range *Range, int NumDims, const int *Begin, const int *End) {

  int iDim;

  Range->nd = NumDims;
  for (iDim = 0; iDim < Range->nd; ++iDim) {
    Range->b[iDim] = Begin[iDim];
    Range->e[iDim] = End[iDim];
  }
  for (iDim = Range->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    Range->b[iDim] = 0;
    Range->e[iDim] = 1;
  }

}

static inline void ovkRangeSize(const ovk_range *Range, int *Size) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    Size[iDim] = Range->e[iDim] - Range->b[iDim];
  }

}

static inline void ovkRangeCount(const ovk_range *Range, size_t *Count) {

  int iDim;

  *Count = 1;

  for (iDim = 0; iDim < Range->nd; ++iDim) {
    *Count *= (size_t)(Range->e[iDim] - Range->b[iDim]);
  }

}

static inline void ovkRangeCountSmall(const ovk_range *Range, int *Count) {

  size_t CountLarge;

  ovkRangeCount(Range, &CountLarge);

  *Count = (int)CountLarge;

}

static inline void ovkRangeTupleToIndex(const ovk_range *Range, ovk_array_layout Layout,
  const int *Tuple, size_t *Index) {

  int iDim;
  size_t Stride;

  switch (Layout) {
  case OVK_COLUMN_MAJOR:
    *Index = 0;
    Stride = 1;
    for (iDim = 0; iDim < Range->nd; ++iDim) {
      *Index += Stride * (size_t)(Tuple[iDim] - Range->b[iDim]);
      Stride *= (size_t)(Range->e[iDim] - Range->b[iDim]);
    }
    break;
  case OVK_ROW_MAJOR:
    *Index = 0;
    Stride = 1;
    for (iDim = Range->nd-1; iDim >= 0; --iDim) {
      *Index += Stride * (size_t)(Tuple[iDim] - Range->b[iDim]);
      Stride *= (size_t)(Range->e[iDim] - Range->b[iDim]);
    }
    break;
  }

}

static inline void ovkRangeTupleToIndexSmall(const ovk_range *Range, ovk_array_layout Layout,
  const int *Tuple, int *Index) {

  size_t IndexLarge;

  ovkRangeTupleToIndex(Range, Layout, Tuple, &IndexLarge);

  *Index = (int)IndexLarge;

}

static inline void ovkRangeIndexToTuple(const ovk_range *Range, ovk_array_layout Layout,
  size_t Index, int *Tuple) {

  int iDim;
  size_t Stride;

  switch (Layout) {
  case OVK_COLUMN_MAJOR:
    Stride = 1;
    for (iDim = 0; iDim < Range->nd; ++iDim) {
      size_t N = Range->e[iDim] - Range->b[iDim];
      Tuple[iDim] = Range->b[iDim] + (int)((Index/Stride) % N);
      Stride *= N;
    }
    break;
  case OVK_ROW_MAJOR:
    Stride = 1;
    for (iDim = Range->nd-1; iDim >= 0; --iDim) {
      size_t N = Range->e[iDim] - Range->b[iDim];
      Tuple[iDim] = Range->b[iDim] + (int)((Index/Stride) % N);
      Stride *= N;
    }
    break;
  }

  for (iDim = Range->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    Tuple[iDim] = 0;
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

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Point[iDim] < Range->b[iDim] || Point[iDim] >= Range->e[iDim]) return false;
  }

  return true;

}

static inline bool ovkRangeIncludes(const ovk_range *Range, const ovk_range *OtherRange) {

  return ovkRangeIsEmpty(OtherRange) ||
    (Range->b[0] <= OtherRange->b[0] && Range->e[0] >= OtherRange->e[0] &&
    Range->b[1] <= OtherRange->b[1] && Range->e[1] >= OtherRange->e[1] &&
    Range->b[2] <= OtherRange->b[2] && Range->e[2] >= OtherRange->e[2]);

}

static inline bool ovkRangeOverlaps(const ovk_range *LeftRange, const ovk_range *RightRange) {

  return !ovkRangeIsEmpty(LeftRange) && !ovkRangeIsEmpty(RightRange) &&
    RightRange->e[0] > LeftRange->b[0] && LeftRange->e[0] > RightRange->b[0] &&
    RightRange->e[1] > LeftRange->b[1] && LeftRange->e[1] > RightRange->b[1] &&
    RightRange->e[2] > LeftRange->b[2] && LeftRange->e[2] > RightRange->b[2];

}

static inline void ovkRangeUnion(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *UnionRange) {

  int iDim;

  if (ovkRangeIsEmpty(LeftRange)) {
    *UnionRange = *RightRange;
  } else if (ovkRangeIsEmpty(RightRange)) {
    *UnionRange = *LeftRange;
  } else {
    UnionRange->nd = LeftRange->nd;
    for (iDim = 0; iDim < UnionRange->nd; ++iDim) {
      UnionRange->b[iDim] = ovk_min(LeftRange->b[iDim], RightRange->b[iDim]);
      UnionRange->e[iDim] = ovk_max(LeftRange->e[iDim], RightRange->e[iDim]);
    }
    for (iDim = UnionRange->nd; iDim < OVK_MAX_DIMS; ++iDim) {
      UnionRange->b[iDim] = 0;
      UnionRange->e[iDim] = 1;
    }
  }

}

static inline void ovkRangeIntersect(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *IntersectRange) {

  int iDim;

  IntersectRange->nd = LeftRange->nd;
  for (iDim = 0; iDim < IntersectRange->nd; ++iDim) {
    IntersectRange->b[iDim] = ovk_max(LeftRange->b[iDim], RightRange->b[iDim]);
    IntersectRange->e[iDim] = ovk_min(LeftRange->e[iDim], RightRange->e[iDim]);
  }
  for (iDim = IntersectRange->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    IntersectRange->b[iDim] = 0;
    IntersectRange->e[iDim] = 1;
  }

}

static inline void ovkRangeClamp(const ovk_range *Range, int *Point) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Point[iDim] < Range->b[iDim]) {
      Point[iDim] = Range->b[iDim];
    } else if (Point[iDim] >= Range->e[iDim]) {
      Point[iDim] = Range->e[iDim]-1;
    }
  }

}

static inline void ovkRangeExtend(const ovk_range *Range, const int *Point, ovk_range *ExtendRange)
  {

  int iDim;

  ExtendRange->nd = Range->nd;
  if (ovkRangeIsEmpty(Range)) {
    for (iDim = 0; iDim < Range->nd; ++iDim) {
      ExtendRange->b[iDim] = Point[iDim];
      ExtendRange->e[iDim] = Point[iDim];
    }
  } else {
    for (iDim = 0; iDim < Range->nd; ++iDim) {
      ExtendRange->b[iDim] = ovk_min(Range->b[iDim], Point[iDim]);
      ExtendRange->e[iDim] = ovk_max(Range->e[iDim], Point[iDim]);
    }
  }
  for (iDim = Range->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    ExtendRange->b[iDim] = 0;
    ExtendRange->e[iDim] = 1;
  }

}

#ifdef __cplusplus
}
#endif

#endif
