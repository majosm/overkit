// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifdef __cplusplus
extern "C" {
#endif

static inline void ovkDefaultRange(ovk_range *Range, int NumDims) {

  int iDim;

  Range->NumDims = NumDims;
  for (iDim = 0; iDim < Range->NumDims; ++iDim) {
    Range->Begin[iDim] = 0;
    Range->End[iDim] = 0;
  }
  for (iDim = Range->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    Range->Begin[iDim] = 0;
    Range->End[iDim] = 1;
  }

}

static inline void ovkSetRange(ovk_range *Range, int NumDims, const int *Begin, const int *End) {

  int iDim;

  Range->NumDims = NumDims;
  for (iDim = 0; iDim < Range->NumDims; ++iDim) {
    Range->Begin[iDim] = Begin[iDim];
    Range->End[iDim] = End[iDim];
  }
  for (iDim = Range->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    Range->Begin[iDim] = 0;
    Range->End[iDim] = 1;
  }

}

static inline bool ovkRangeEquals(const ovk_range *LeftRange, const ovk_range *RightRange) {

  return LeftRange->NumDims == RightRange->NumDims &&
    LeftRange->Begin[0] == RightRange->Begin[0] &&
    LeftRange->Begin[1] == RightRange->Begin[1] &&
    LeftRange->Begin[2] == RightRange->Begin[2] &&
    LeftRange->End[0] == RightRange->End[0] &&
    LeftRange->End[1] == RightRange->End[1] &&
    LeftRange->End[2] == RightRange->End[2];

}

static inline void ovkRangeSize(const ovk_range *Range, int *Size) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    Size[iDim] = Range->End[iDim] - Range->Begin[iDim];
  }

}

static inline void ovkRangeCount(const ovk_range *Range, long long *Count) {

  int iDim;

  *Count = 1;

  for (iDim = 0; iDim < Range->NumDims; ++iDim) {
    *Count *= Range->End[iDim] - Range->Begin[iDim];
  }

}

static inline bool ovkRangeIsEmpty(const ovk_range *Range) {

  return
    Range->End[0] <= Range->Begin[0] ||
    Range->End[1] <= Range->Begin[1] ||
    Range->End[2] <= Range->Begin[2];

}

static inline void ovkRangeTupleToIndex(const ovk_range *Range, ovk_array_layout Layout,
  const int *Tuple, long long *Index) {

  int iDim;

  long long Offset[OVK_MAX_DIMS];
  for (iDim = 0; iDim < Range->NumDims; ++iDim) {
    Offset[iDim] = Tuple[iDim]-Range->Begin[iDim];
  }
  for (iDim = Range->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    Offset[iDim] = 0;
  }

  long long Size[OVK_MAX_DIMS];

  switch (Layout) {
  case OVK_ROW_MAJOR:
    Size[1] = Range->End[1]-Range->Begin[1];
    Size[2] = Range->End[2]-Range->Begin[2];
    *Index = Size[2]*(Size[1]*Offset[0] + Offset[1]) + Offset[2];
    break;
  case OVK_COLUMN_MAJOR:
    Size[0] = Range->End[0]-Range->Begin[0];
    Size[1] = Range->End[1]-Range->Begin[1];
    *Index = Offset[0] + Size[0]*(Offset[1] + Size[1]*Offset[2]);
    break;
  }

}

static inline void ovkRangeIndexToTuple(const ovk_range *Range, ovk_array_layout Layout,
  long long Index, int *Tuple) {

  int iDim;
  long long Stride[OVK_MAX_DIMS];
  long long ReducedIndex;

  switch (Layout) {
  case OVK_ROW_MAJOR:
    Stride[2] = 1;
    Stride[1] = Range->End[2]-Range->Begin[2];
    Stride[0] = Stride[1]*(Range->End[1]-Range->Begin[1]);
    ReducedIndex = Index;
    for (iDim = 0; iDim < Range->NumDims; ++iDim) {
      int Offset = (int)(ReducedIndex/Stride[iDim]);
      Tuple[iDim] = Range->Begin[iDim] + Offset;
      ReducedIndex -= Offset*Stride[iDim];
    }
    break;
  case OVK_COLUMN_MAJOR:
    Stride[0] = 1;
    Stride[1] = Range->End[0]-Range->Begin[0];
    Stride[2] = Stride[1]*(Range->End[1]-Range->Begin[1]);
    ReducedIndex = Index;
    for (iDim = Range->NumDims-1; iDim >= 0; --iDim) {
      int Offset = (int)(ReducedIndex/Stride[iDim]);
      Tuple[iDim] = Range->Begin[iDim] + Offset;
      ReducedIndex -= Offset*Stride[iDim];
    }
    break;
  }

  for (iDim = Range->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    Tuple[iDim] = 0;
  }

}

static inline bool ovkRangeContains(const ovk_range *Range, const int *Point) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Point[iDim] < Range->Begin[iDim] || Point[iDim] >= Range->End[iDim]) return false;
  }

  return true;

}

static inline bool ovkRangeIncludes(const ovk_range *Range, const ovk_range *OtherRange) {

  return ovkRangeIsEmpty(OtherRange) ||
    (Range->Begin[0] <= OtherRange->Begin[0] && Range->End[0] >= OtherRange->End[0] &&
    Range->Begin[1] <= OtherRange->Begin[1] && Range->End[1] >= OtherRange->End[1] &&
    Range->Begin[2] <= OtherRange->Begin[2] && Range->End[2] >= OtherRange->End[2]);

}

static inline bool ovkRangeOverlaps(const ovk_range *LeftRange, const ovk_range *RightRange) {

  return !ovkRangeIsEmpty(LeftRange) && !ovkRangeIsEmpty(RightRange) &&
    RightRange->End[0] > LeftRange->Begin[0] && LeftRange->End[0] > RightRange->Begin[0] &&
    RightRange->End[1] > LeftRange->Begin[1] && LeftRange->End[1] > RightRange->Begin[1] &&
    RightRange->End[2] > LeftRange->Begin[2] && LeftRange->End[2] > RightRange->Begin[2];

}

static inline void ovkRangeUnion(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *UnionRange) {

  int iDim;

  if (ovkRangeIsEmpty(LeftRange)) {
    *UnionRange = *RightRange;
  } else if (ovkRangeIsEmpty(RightRange)) {
    *UnionRange = *LeftRange;
  } else {
    UnionRange->NumDims = LeftRange->NumDims;
    for (iDim = 0; iDim < UnionRange->NumDims; ++iDim) {
      UnionRange->Begin[iDim] = ovk_min(LeftRange->Begin[iDim], RightRange->Begin[iDim]);
      UnionRange->End[iDim] = ovk_max(LeftRange->End[iDim], RightRange->End[iDim]);
    }
    for (iDim = UnionRange->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
      UnionRange->Begin[iDim] = 0;
      UnionRange->End[iDim] = 1;
    }
  }

}

static inline void ovkRangeIntersect(const ovk_range *LeftRange, const ovk_range *RightRange,
  ovk_range *IntersectRange) {

  int iDim;

  IntersectRange->NumDims = LeftRange->NumDims;
  for (iDim = 0; iDim < IntersectRange->NumDims; ++iDim) {
    IntersectRange->Begin[iDim] = ovk_max(LeftRange->Begin[iDim], RightRange->Begin[iDim]);
    IntersectRange->End[iDim] = ovk_min(LeftRange->End[iDim], RightRange->End[iDim]);
  }
  for (iDim = IntersectRange->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    IntersectRange->Begin[iDim] = 0;
    IntersectRange->End[iDim] = 1;
  }

}

static inline void ovkRangeClamp(const ovk_range *Range, int *Point) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Point[iDim] < Range->Begin[iDim]) {
      Point[iDim] = Range->Begin[iDim];
    } else if (Point[iDim] >= Range->End[iDim]) {
      Point[iDim] = Range->End[iDim]-1;
    }
  }

}

static inline void ovkRangeExtend(const ovk_range *Range, const int *Point, ovk_range *ExtendRange)
  {

  int iDim;

  ExtendRange->NumDims = Range->NumDims;
  if (ovkRangeIsEmpty(Range)) {
    for (iDim = 0; iDim < Range->NumDims; ++iDim) {
      ExtendRange->Begin[iDim] = Point[iDim];
      ExtendRange->End[iDim] = Point[iDim]+1;
    }
  } else {
    for (iDim = 0; iDim < Range->NumDims; ++iDim) {
      ExtendRange->Begin[iDim] = ovk_min(Range->Begin[iDim], Point[iDim]);
      ExtendRange->End[iDim] = ovk_max(Range->End[iDim], Point[iDim]+1);
    }
  }
  for (iDim = Range->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    ExtendRange->Begin[iDim] = 0;
    ExtendRange->End[iDim] = 1;
  }

}

#ifdef __cplusplus
}
#endif
