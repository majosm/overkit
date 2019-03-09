// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RANGE_BASE_H_INCLUDED
#define OVK_CORE_RANGE_BASE_H_INCLUDED

#include <ovk/core/ConstantsBase.h>
#include <ovk/core/GlobalBase.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int NumDims_;
  int Begin_[OVK_MAX_DIMS];
  int End_[OVK_MAX_DIMS];
} ovk_range;

static inline ovk_range ovkDefaultRange(int NumDims);
static inline ovk_range ovkRange(int NumDims, const int *Begin, const int *End);
static inline int ovkRangeDimension(const ovk_range *Range);
static inline void ovkRangeBegin(const ovk_range *Range, int N, int *Begin);
static inline int ovkRangeBeginDim(const ovk_range *Range, int iDim);
static inline const int *ovkRangeBeginData(const ovk_range *Range);
static inline void ovkRangeEnd(const ovk_range *Range, int N, int *End);
static inline int ovkRangeEndDim(const ovk_range *Range, int iDim);
static inline const int *ovkRangeEndData(const ovk_range *Range);
static inline void ovkRangeSize(const ovk_range *Range, int N, int *Size);
static inline int ovkRangeSizeDim(const ovk_range *Range, int iDim);
static inline long long ovkRangeCount(const ovk_range *Range);
static inline bool ovkRangeIsEmpty(const ovk_range *Range);
static inline bool ovkRangesAreEqual(const ovk_range *LeftRange, const ovk_range *RightRange);
static inline long long ovkRangeTupleToIndex(const ovk_range *Range, ovk_array_layout Layout,
  const int *Tuple);
static inline void ovkRangeIndexToTuple(const ovk_range *Range, ovk_array_layout Layout,
  long long Index, int N, int *Tuple);
static inline bool ovkRangeContains(const ovk_range *Range, const int *Tuple);
static inline bool ovkRangeIncludes(const ovk_range *Range, const ovk_range *OtherRange);
static inline bool ovkRangesOverlap(const ovk_range *LeftRange, const ovk_range *RightRange);
static inline ovk_range ovkUnionRanges(const ovk_range *LeftRange, const ovk_range *RightRange);
static inline ovk_range ovkIntersectRanges(const ovk_range *LeftRange, const ovk_range *RightRange);
static inline void ovkClampToRange(int *Tuple, const ovk_range *Range);
static inline void ovkExtendRange(ovk_range *Range, const int *Tuple);

#ifdef __cplusplus
}
#endif

#include <ovk/core/RangeBase.inl>

#endif
