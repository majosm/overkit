// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_BOX_INCLUDED
#define OVK_CORE_PUBLIC_BOX_INCLUDED

#include <ovkGlobal.h>

struct ovk_box;
typedef struct ovk_box ovk_box;

static inline void ovkBoxDefault(ovk_box *Box, int NumDims);
static inline void ovkBoxSize(const ovk_box *Box, double *Size);
static inline void ovkBoxVolume(const ovk_box *Box, double *Volume);
static inline bool ovkBoxIsEmpty(const ovk_box *Box);
static inline bool ovkBoxEquals(const ovk_box *LeftBox, const ovk_box *RightBox);
static inline bool ovkBoxContains(const ovk_box *Box, const double *Point);
static inline bool ovkBoxOverlaps(const ovk_box *LeftBox, const ovk_box *RightBox);
static inline void ovkBoxUnion(const ovk_box *LeftBox, const ovk_box *RightBox, ovk_box *UnionBox);
static inline void ovkBoxIntersect(const ovk_box *LeftBox, const ovk_box *RightBox,
  ovk_box *IntersectBox);
static inline void ovkBoxClamp(const ovk_box *Box, const double *Point, double *ClampedPoint);
static inline void ovkBoxMove(const ovk_box *Box, const double *Amount, ovk_box *MoveBox);
static inline void ovkBoxGrow(const ovk_box *Box, const double *Amount, ovk_box *GrowBox);
static inline void ovkBoxScale(const ovk_box *Box, const double *Factor, ovk_box *ScaleBox);
static inline void ovkBoxExtend(const ovk_box *Box, const double *Point, ovk_box *ExtendBox);

#endif
