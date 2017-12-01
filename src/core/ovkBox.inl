// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_BOX_INL_INCLUDED
#define OVK_CORE_PUBLIC_BOX_INL_INCLUDED

#include <ovkBox.h>

#include <ovkGlobal.h>

struct ovk_box {
  int nd;
  double b[OVK_MAX_DIMS];
  double e[OVK_MAX_DIMS];
};

static inline void ovkBoxDefault(ovk_box *Box, int NumDims) {

  Box->nd = NumDims;
  for (int d = 0; d < NumDims; ++d) {
    Box->b[d] = 0.;
    Box->e[d] = -1.;
  }
  for (int d = NumDims; d < OVK_MAX_DIMS; ++d) {
    Box->b[d] = 0.;
    Box->e[d] = 0.;
  }

}

static inline void ovkBoxSize(const ovk_box *Box, double *Size) {

  for (int d = 0; d < Box->nd; ++d) {
    Size[d] = Box->e[d] - Box->b[d];
  }

}

static inline void ovkBoxVolume(const ovk_box *Box, double *Volume) {

  *Volume = 1.;

  for (int d = 0; d < Box->nd; ++d) {
    *Volume *= Box->e[d] - Box->b[d];
  }

}

static inline bool ovkBoxIsEmpty(const ovk_box *Box) {

  return Box->e[0] < Box->b[0] || Box->e[1] < Box->b[1] || Box->e[2] < Box->b[2];

}

static inline bool ovkBoxEquals(const ovk_box *LeftBox, const ovk_box *RightBox) {

  return LeftBox->nd == RightBox->nd &&
    LeftBox->b[0] == RightBox->b[0] &&
    LeftBox->b[1] == RightBox->b[1] &&
    LeftBox->b[2] == RightBox->b[2] &&
    LeftBox->e[0] == RightBox->e[0] &&
    LeftBox->e[1] == RightBox->e[1] &&
    LeftBox->e[2] == RightBox->e[2];

}

static inline bool ovkBoxContains(const ovk_box *Box, const double *Point) {

  for (int d = 0; d < Box->nd; ++d) {
    if (Point[d] < Box->b[d] || Point[d] > Box->e[d]) return false;
  }

  return true;

}

static inline bool ovkBoxOverlaps(const ovk_box *LeftBox, const ovk_box *RightBox) {

  return !ovkBoxIsEmpty(LeftBox) && !ovkBoxIsEmpty(RightBox) &&
    RightBox->e[0] >= LeftBox->b[0] && LeftBox->e[0] >= RightBox->b[0] &&
    RightBox->e[1] >= LeftBox->b[1] && LeftBox->e[1] >= RightBox->b[1] &&
    RightBox->e[2] >= LeftBox->b[2] && LeftBox->e[2] >= RightBox->b[2];

}

static inline void ovkBoxUnion(const ovk_box *LeftBox, const ovk_box *RightBox, ovk_box *UnionBox) {

  if (ovkBoxIsEmpty(LeftBox)) {
    *UnionBox = *RightBox;
  } else if (ovkBoxIsEmpty(RightBox)) {
    *UnionBox = *LeftBox;
  } else {
    UnionBox->nd = LeftBox->nd;
    for (int d = 0; d < UnionBox->nd; ++d) {
      UnionBox->b[d] = ovk_min(LeftBox->b[d], RightBox->b[d]);
      UnionBox->e[d] = ovk_max(LeftBox->e[d], RightBox->e[d]);
    }
    for (int d = UnionBox->nd; d < OVK_MAX_DIMS; ++d) {
      UnionBox->b[d] = 0.;
      UnionBox->e[d] = 0.;
    }
  }

}

static inline void ovkBoxIntersect(const ovk_box *LeftBox, const ovk_box *RightBox,
  ovk_box *IntersectBox) {

  IntersectBox->nd = LeftBox->nd;
  for (int d = 0; d < IntersectBox->nd; ++d) {
    IntersectBox->b[d] = ovk_max(LeftBox->b[d], RightBox->b[d]);
    IntersectBox->e[d] = ovk_min(LeftBox->e[d], RightBox->e[d]);
  }
  for (int d = IntersectBox->nd; d < OVK_MAX_DIMS; ++d) {
    IntersectBox->b[d] = 0.;
    IntersectBox->e[d] = 0.;
  }

}

static inline void ovkBoxClamp(const ovk_box *Box, const double *Point, double *ClampedPoint) {

  for (int d = 0; d < Box->nd; ++d) {
    if (Point[d] < Box->b[d]) {
      ClampedPoint[d] = Box->b[d];
    } else if (Point[d] > Box->e[d]) {
      ClampedPoint[d] = Box->e[d];
    } else {
      ClampedPoint[d] = Point[d];
    }
  }

}

static inline void ovkBoxMove(const ovk_box *Box, const double *Amount, ovk_box *MoveBox) {

  MoveBox->nd = Box->nd;
  for (int d = 0; d < Box->nd; ++d) {
    MoveBox->b[d] = Box->b[d] + Amount[d];
    MoveBox->e[d] = Box->e[d] + Amount[d];
  }
  for (int d = Box->nd; d < OVK_MAX_DIMS; ++d) {
    MoveBox->b[d] = 0.;
    MoveBox->e[d] = 0.;
  }

}

static inline void ovkBoxGrow(const ovk_box *Box, const double *Amount, ovk_box *GrowBox) {

  if (ovkBoxIsEmpty(Box)) {
    *GrowBox = *Box;
  } else {
    GrowBox->nd = Box->nd;
    for (int d = 0; d < Box->nd; ++d) {
      GrowBox->b[d] = Box->b[d] - Amount[d];
      GrowBox->e[d] = Box->e[d] + Amount[d];
    }
    for (int d = Box->nd; d < OVK_MAX_DIMS; ++d) {
      GrowBox->b[d] = 0.;
      GrowBox->e[d] = 0.;
    }
  }

}

static inline void ovkBoxScale(const ovk_box *Box, const double *Factor, ovk_box *ScaleBox) {

  ScaleBox->nd = Box->nd;
  for (int d = 0; d < Box->nd; ++d) {
    double Center = 0.5 * (Box->e[d] + Box->b[d]);
    double HalfSize = 0.5 * Factor[d] * (Box->e[d] - Box->b[d]);
    ScaleBox->b[d] = Center - HalfSize;
    ScaleBox->e[d] = Center + HalfSize;
  }
  for (int d = Box->nd; d < OVK_MAX_DIMS; ++d) {
    ScaleBox->b[d] = 0.;
    ScaleBox->e[d] = 0.;
  }

}

static inline void ovkBoxExtend(const ovk_box *Box, const double *Point, ovk_box *ExtendBox) {

  ExtendBox->nd = Box->nd;
  for (int d = 0; d < Box->nd; ++d) {
    ExtendBox->b[d] = ovk_min(ExtendBox->b[d], Point[d]);
    ExtendBox->e[d] = ovk_max(ExtendBox->e[d], Point[d]);
  }
  for (int d = Box->nd; d < OVK_MAX_DIMS; ++d) {
    ExtendBox->b[d] = 0.;
    ExtendBox->e[d] = 0.;
  }

}

#endif
