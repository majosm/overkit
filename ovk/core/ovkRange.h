// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_RANGE_INCLUDED
#define OVK_CORE_PUBLIC_RANGE_INCLUDED

#include <ovk/core/ovkGlobal.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int nd;
  int b[OVK_MAX_DIMS];
  int e[OVK_MAX_DIMS];
} ovk_range;

static inline void ovkDefaultRange(ovk_range *Range, int NumDims);
static inline void ovkSetRange(ovk_range *Range, int NumDims, const int *Begin, const int *End);
static inline void ovkRangeSize(const ovk_range *Range, int *Size);
static inline void ovkRangeCount(const ovk_range *Range, size_t *Count);
static inline void ovkRangeCountSmall(const ovk_range *Range, int *Count);
static inline void ovkRangeTupleToIndex(const ovk_range *Range, ovk_array_layout Layout,
  const int *Tuple, size_t *Index);
static inline void ovkRangeTupleToIndexSmall(const ovk_range *Range, ovk_array_layout Layout,
  const int *Tuple, int *Index);
static inline void ovkRangeIndexToTuple(const ovk_range *Range, ovk_array_layout Layout,
  size_t Index, int *Tuple);
static inline bool ovkRangeIsEmpty(const ovk_range *Range);
static inline bool ovkRangeEquals(const ovk_range *LeftRange, const ovk_range *RightRange);
static inline bool ovkRangeContains(const ovk_range *Range, const int *Point);
static inline bool ovkRangeIncludes(const ovk_range *LeftRange, const ovk_range *RightRange);
static inline bool ovkRangeOverlaps(const ovk_range *LeftRange, const ovk_range *RightRange);
static inline void ovkRangeUnion(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *UnionRange);
static inline void ovkRangeIntersect(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *IntersectRange);
static inline void ovkRangeClamp(const ovk_range *Range, int *Point);
static inline void ovkRangeExtend(const ovk_range *Range, const int *Point, ovk_range *ExtendRange);

#ifdef __cplusplus
}
#endif

#endif
