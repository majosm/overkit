// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_RANGE_INCLUDED
#define OVK_CORE_PUBLIC_RANGE_INCLUDED

#include <ovkGlobal.h>

struct ovk_range;
typedef struct ovk_range ovk_range;

static inline void ovkRangeDefault(ovk_range *Range, int NumDims);
static inline void ovkRangeSize(const ovk_range *Range, int *Size);
static inline void ovkRangeCount(const ovk_range *Range, size_t *Count);
static inline void ovkRangeTupleToIndex(const ovk_range *Range, const int *Tuple, size_t *Index);
static inline void ovkRangeIndexToTuple(const ovk_range *Range, size_t Index, int *Tuple);
static inline bool ovkRangeIsEmpty(const ovk_range *Range);
static inline bool ovkRangeEquals(const ovk_range *LeftRange, const ovk_range *RightRange);
static inline bool ovkRangeContains(const ovk_range *Range, const int *Point);
static inline bool ovkRangeOverlaps(const ovk_range *LeftRange, const ovk_range *RightRange);
static inline void ovkRangeUnion(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *UnionRange);
static inline void ovkRangeIntersect(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *IntersectRange);
static inline void ovkRangeClamp(const ovk_range *Range, const int *Point, int *ClampedPoint);

#endif
