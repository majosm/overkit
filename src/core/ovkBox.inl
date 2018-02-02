// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_BOX_INL_INCLUDED
#define OVK_CORE_PUBLIC_BOX_INL_INCLUDED

#include <ovkBox.h>

#include <ovkGlobal.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline void ovkDefaultBox(ovk_box *Box, int NumDims) {

  int iDim;

  Box->nd = NumDims;
  for (iDim = 0; iDim < Box->nd; ++iDim) {
    Box->b[iDim] = 0.;
    Box->e[iDim] = -1.;
  }
  for (iDim = Box->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    Box->b[iDim] = 0.;
    Box->e[iDim] = 0.;
  }

}

static inline void ovkSetBox(ovk_box *Box, int NumDims, const double *Begin, const double *End) {

  int iDim;

  Box->nd = NumDims;
  for (iDim = 0; iDim < Box->nd; ++iDim) {
    Box->b[iDim] = Begin[iDim];
    Box->e[iDim] = End[iDim];
  }
  for (iDim = Box->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    Box->b[iDim] = 0.;
    Box->e[iDim] = 0.;
  }

}

static inline void ovkBoxSize(const ovk_box *Box, double *Size) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    Size[iDim] = Box->e[iDim] - Box->b[iDim];
  }

}

static inline void ovkBoxVolume(const ovk_box *Box, double *Volume) {

  int iDim;

  *Volume = 1.;

  for (iDim = 0; iDim < Box->nd; ++iDim) {
    *Volume *= Box->e[iDim] - Box->b[iDim];
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

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Point[iDim] < Box->b[iDim] || Point[iDim] > Box->e[iDim]) return false;
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

  int iDim;

  if (ovkBoxIsEmpty(LeftBox)) {
    *UnionBox = *RightBox;
  } else if (ovkBoxIsEmpty(RightBox)) {
    *UnionBox = *LeftBox;
  } else {
    UnionBox->nd = LeftBox->nd;
    for (iDim = 0; iDim < UnionBox->nd; ++iDim) {
      UnionBox->b[iDim] = ovk_min(LeftBox->b[iDim], RightBox->b[iDim]);
      UnionBox->e[iDim] = ovk_max(LeftBox->e[iDim], RightBox->e[iDim]);
    }
    for (iDim = UnionBox->nd; iDim < OVK_MAX_DIMS; ++iDim) {
      UnionBox->b[iDim] = 0.;
      UnionBox->e[iDim] = 0.;
    }
  }

}

static inline void ovkBoxIntersect(const ovk_box *LeftBox, const ovk_box *RightBox,
  ovk_box *IntersectBox) {

  int iDim;

  IntersectBox->nd = LeftBox->nd;
  for (iDim = 0; iDim < IntersectBox->nd; ++iDim) {
    IntersectBox->b[iDim] = ovk_max(LeftBox->b[iDim], RightBox->b[iDim]);
    IntersectBox->e[iDim] = ovk_min(LeftBox->e[iDim], RightBox->e[iDim]);
  }
  for (iDim = IntersectBox->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    IntersectBox->b[iDim] = 0.;
    IntersectBox->e[iDim] = 0.;
  }

}

static inline void ovkBoxClamp(const ovk_box *Box, double *Point) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Point[iDim] < Box->b[iDim]) {
      Point[iDim] = Box->b[iDim];
    } else if (Point[iDim] > Box->e[iDim]) {
      Point[iDim] = Box->e[iDim];
    }
  }

}

static inline void ovkBoxMove(const ovk_box *Box, const double *Amount, ovk_box *MoveBox) {

  int iDim;

  MoveBox->nd = Box->nd;
  for (iDim = 0; iDim < Box->nd; ++iDim) {
    MoveBox->b[iDim] = Box->b[iDim] + Amount[iDim];
    MoveBox->e[iDim] = Box->e[iDim] + Amount[iDim];
  }
  for (iDim = Box->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    MoveBox->b[iDim] = 0.;
    MoveBox->e[iDim] = 0.;
  }

}

static inline void ovkBoxGrow(const ovk_box *Box, const double *Amount, ovk_box *GrowBox) {

  int iDim;

  if (ovkBoxIsEmpty(Box)) {
    *GrowBox = *Box;
  } else {
    GrowBox->nd = Box->nd;
    for (iDim = 0; iDim < Box->nd; ++iDim) {
      GrowBox->b[iDim] = Box->b[iDim] - Amount[iDim];
      GrowBox->e[iDim] = Box->e[iDim] + Amount[iDim];
    }
    for (iDim = Box->nd; iDim < OVK_MAX_DIMS; ++iDim) {
      GrowBox->b[iDim] = 0.;
      GrowBox->e[iDim] = 0.;
    }
  }

}

static inline void ovkBoxScale(const ovk_box *Box, const double *Factor, ovk_box *ScaleBox) {

  int iDim;

  ScaleBox->nd = Box->nd;
  for (iDim = 0; iDim < Box->nd; ++iDim) {
    double Center = 0.5 * (Box->e[iDim] + Box->b[iDim]);
    double HalfSize = 0.5 * Factor[iDim] * (Box->e[iDim] - Box->b[iDim]);
    ScaleBox->b[iDim] = Center - HalfSize;
    ScaleBox->e[iDim] = Center + HalfSize;
  }
  for (iDim = Box->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    ScaleBox->b[iDim] = 0.;
    ScaleBox->e[iDim] = 0.;
  }

}

static inline void ovkBoxExtend(const ovk_box *Box, const double *Point, ovk_box *ExtendBox) {

  int iDim;

  ExtendBox->nd = Box->nd;
  if (ovkBoxIsEmpty(Box)) {
    for (iDim = 0; iDim < Box->nd; ++iDim) {
      ExtendBox->b[iDim] = Point[iDim];
      ExtendBox->e[iDim] = Point[iDim];
    }
  } else {
    for (iDim = 0; iDim < Box->nd; ++iDim) {
      ExtendBox->b[iDim] = ovk_min(Box->b[iDim], Point[iDim]);
      ExtendBox->e[iDim] = ovk_max(Box->e[iDim], Point[iDim]);
    }
  }
  for (iDim = Box->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    ExtendBox->b[iDim] = 0.;
    ExtendBox->e[iDim] = 0.;
  }

}

#ifdef __cplusplus
}
#endif

#endif
