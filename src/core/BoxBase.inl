// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifdef __cplusplus
extern "C" {
#endif

static inline void ovkDefaultBox(ovk_box *Box, int NumDims) {

  int iDim;

  Box->NumDims = NumDims;
  for (iDim = 0; iDim < Box->NumDims; ++iDim) {
    Box->Begin[iDim] = 0.;
    Box->End[iDim] = -1.;
  }
  for (iDim = Box->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    Box->Begin[iDim] = 0.;
    Box->End[iDim] = 0.;
  }

}

static inline void ovkSetBox(ovk_box *Box, int NumDims, const double *Begin, const double *End) {

  int iDim;

  Box->NumDims = NumDims;
  for (iDim = 0; iDim < Box->NumDims; ++iDim) {
    Box->Begin[iDim] = Begin[iDim];
    Box->End[iDim] = End[iDim];
  }
  for (iDim = Box->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    Box->Begin[iDim] = 0.;
    Box->End[iDim] = 0.;
  }

}

static inline void ovkBoxSize(const ovk_box *Box, double *Size) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    Size[iDim] = Box->End[iDim] - Box->Begin[iDim];
  }

}

static inline void ovkBoxVolume(const ovk_box *Box, double *Volume) {

  int iDim;

  *Volume = 1.;

  for (iDim = 0; iDim < Box->NumDims; ++iDim) {
    *Volume *= Box->End[iDim] - Box->Begin[iDim];
  }

}

static inline bool ovkBoxIsEmpty(const ovk_box *Box) {

  return Box->End[0] < Box->Begin[0] || Box->End[1] < Box->Begin[1] || Box->End[2] < Box->Begin[2];

}

static inline bool ovkBoxEquals(const ovk_box *LeftBox, const ovk_box *RightBox) {

  return LeftBox->NumDims == RightBox->NumDims &&
    LeftBox->Begin[0] == RightBox->Begin[0] &&
    LeftBox->Begin[1] == RightBox->Begin[1] &&
    LeftBox->Begin[2] == RightBox->Begin[2] &&
    LeftBox->End[0] == RightBox->End[0] &&
    LeftBox->End[1] == RightBox->End[1] &&
    LeftBox->End[2] == RightBox->End[2];

}

static inline bool ovkBoxContains(const ovk_box *Box, const double *Point) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Point[iDim] < Box->Begin[iDim] || Point[iDim] > Box->End[iDim]) return false;
  }

  return true;

}

static inline bool ovkBoxOverlaps(const ovk_box *LeftBox, const ovk_box *RightBox) {

  return !ovkBoxIsEmpty(LeftBox) && !ovkBoxIsEmpty(RightBox) &&
    RightBox->End[0] >= LeftBox->Begin[0] && LeftBox->End[0] >= RightBox->Begin[0] &&
    RightBox->End[1] >= LeftBox->Begin[1] && LeftBox->End[1] >= RightBox->Begin[1] &&
    RightBox->End[2] >= LeftBox->Begin[2] && LeftBox->End[2] >= RightBox->Begin[2];

}

static inline void ovkBoxUnion(const ovk_box *LeftBox, const ovk_box *RightBox, ovk_box *UnionBox) {

  int iDim;

  if (ovkBoxIsEmpty(LeftBox)) {
    *UnionBox = *RightBox;
  } else if (ovkBoxIsEmpty(RightBox)) {
    *UnionBox = *LeftBox;
  } else {
    UnionBox->NumDims = LeftBox->NumDims;
    for (iDim = 0; iDim < UnionBox->NumDims; ++iDim) {
      UnionBox->Begin[iDim] = ovk_min(LeftBox->Begin[iDim], RightBox->Begin[iDim]);
      UnionBox->End[iDim] = ovk_max(LeftBox->End[iDim], RightBox->End[iDim]);
    }
    for (iDim = UnionBox->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
      UnionBox->Begin[iDim] = 0.;
      UnionBox->End[iDim] = 0.;
    }
  }

}

static inline void ovkBoxIntersect(const ovk_box *LeftBox, const ovk_box *RightBox, ovk_box
  *IntersectBox) {

  int iDim;

  IntersectBox->NumDims = LeftBox->NumDims;
  for (iDim = 0; iDim < IntersectBox->NumDims; ++iDim) {
    IntersectBox->Begin[iDim] = ovk_max(LeftBox->Begin[iDim], RightBox->Begin[iDim]);
    IntersectBox->End[iDim] = ovk_min(LeftBox->End[iDim], RightBox->End[iDim]);
  }
  for (iDim = IntersectBox->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    IntersectBox->Begin[iDim] = 0.;
    IntersectBox->End[iDim] = 0.;
  }

}

static inline void ovkBoxClamp(const ovk_box *Box, double *Point) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Point[iDim] < Box->Begin[iDim]) {
      Point[iDim] = Box->Begin[iDim];
    } else if (Point[iDim] > Box->End[iDim]) {
      Point[iDim] = Box->End[iDim];
    }
  }

}

static inline void ovkBoxMove(const ovk_box *Box, const double *Amount, ovk_box *MoveBox) {

  int iDim;

  MoveBox->NumDims = Box->NumDims;
  for (iDim = 0; iDim < Box->NumDims; ++iDim) {
    MoveBox->Begin[iDim] = Box->Begin[iDim] + Amount[iDim];
    MoveBox->End[iDim] = Box->End[iDim] + Amount[iDim];
  }
  for (iDim = Box->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    MoveBox->Begin[iDim] = 0.;
    MoveBox->End[iDim] = 0.;
  }

}

static inline void ovkBoxGrow(const ovk_box *Box, const double *Amount, ovk_box *GrowBox) {

  int iDim;

  if (ovkBoxIsEmpty(Box)) {
    *GrowBox = *Box;
  } else {
    GrowBox->NumDims = Box->NumDims;
    for (iDim = 0; iDim < Box->NumDims; ++iDim) {
      GrowBox->Begin[iDim] = Box->Begin[iDim] - Amount[iDim];
      GrowBox->End[iDim] = Box->End[iDim] + Amount[iDim];
    }
    for (iDim = Box->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
      GrowBox->Begin[iDim] = 0.;
      GrowBox->End[iDim] = 0.;
    }
  }

}

static inline void ovkBoxScale(const ovk_box *Box, const double *Factor, ovk_box *ScaleBox) {

  int iDim;

  ScaleBox->NumDims = Box->NumDims;
  for (iDim = 0; iDim < Box->NumDims; ++iDim) {
    double Center = 0.5 * (Box->End[iDim] + Box->Begin[iDim]);
    double HalfSize = 0.5 * Factor[iDim] * (Box->End[iDim] - Box->Begin[iDim]);
    ScaleBox->Begin[iDim] = Center - HalfSize;
    ScaleBox->End[iDim] = Center + HalfSize;
  }
  for (iDim = Box->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    ScaleBox->Begin[iDim] = 0.;
    ScaleBox->End[iDim] = 0.;
  }

}

static inline void ovkBoxExtend(const ovk_box *Box, const double *Point, ovk_box *ExtendBox) {

  int iDim;

  ExtendBox->NumDims = Box->NumDims;
  if (ovkBoxIsEmpty(Box)) {
    for (iDim = 0; iDim < Box->NumDims; ++iDim) {
      ExtendBox->Begin[iDim] = Point[iDim];
      ExtendBox->End[iDim] = Point[iDim];
    }
  } else {
    for (iDim = 0; iDim < Box->NumDims; ++iDim) {
      ExtendBox->Begin[iDim] = ovk_min(Box->Begin[iDim], Point[iDim]);
      ExtendBox->End[iDim] = ovk_max(Box->End[iDim], Point[iDim]);
    }
  }
  for (iDim = Box->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    ExtendBox->Begin[iDim] = 0.;
    ExtendBox->End[iDim] = 0.;
  }

}

#ifdef __cplusplus
}
#endif
