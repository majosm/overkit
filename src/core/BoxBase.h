// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_BOX_BASE_H_INCLUDED
#define OVK_CORE_BOX_BASE_H_INCLUDED

#include <ovk/core/ConstantsBase.h>
#include <ovk/core/GlobalBase.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int NumDims;
  double Begin[OVK_MAX_DIMS];
  double End[OVK_MAX_DIMS];
} ovk_box;

static inline void ovkDefaultBox(ovk_box *Box, int NumDims);
static inline void ovkSetBox(ovk_box *Box, int NumDims, const double *Begin, const double *End);
static inline bool ovkBoxEquals(const ovk_box *LeftBox, const ovk_box *RightBox);
static inline void ovkBoxSize(const ovk_box *Box, double *Size);
static inline void ovkBoxVolume(const ovk_box *Box, double *Volume);
static inline bool ovkBoxIsEmpty(const ovk_box *Box);
static inline bool ovkBoxContains(const ovk_box *Box, const double *Point);
static inline bool ovkBoxOverlaps(const ovk_box *LeftBox, const ovk_box *RightBox);
static inline void ovkBoxUnion(const ovk_box *LeftBox, const ovk_box *RightBox, ovk_box *UnionBox);
static inline void ovkBoxIntersect(const ovk_box *LeftBox, const ovk_box *RightBox,
  ovk_box *IntersectBox);
static inline void ovkBoxClamp(const ovk_box *Box, double *Point);
static inline void ovkBoxMove(const ovk_box *Box, const double *Amount, ovk_box *MoveBox);
static inline void ovkBoxGrow(const ovk_box *Box, const double *Amount, ovk_box *GrowBox);
static inline void ovkBoxScale(const ovk_box *Box, const double *Factor, ovk_box *ScaleBox);
static inline void ovkBoxExtend(const ovk_box *Box, const double *Point, ovk_box *ExtendBox);

#ifdef __cplusplus
}
#endif

#include <ovk/core/BoxBase.inl>

#endif
