// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifdef __cplusplus
extern "C" {
#endif

static inline ovk_range ovkDefaultRange(int NumDims) {

  int iDim;

  ovk_range Range;

  Range.NumDims_ = NumDims;
  for (iDim = 0; iDim < Range.NumDims_; ++iDim) {
    Range.Begin_[iDim] = 0;
    Range.End_[iDim] = 0;
  }
  for (iDim = Range.NumDims_; iDim < OVK_MAX_DIMS; ++iDim) {
    Range.Begin_[iDim] = 0;
    Range.End_[iDim] = 1;
  }

  return Range;

}

static inline ovk_range ovkRange(int NumDims, const int *Begin, const int *End) {

  int iDim;

  ovk_range Range;

  Range.NumDims_ = NumDims;
  for (iDim = 0; iDim < Range.NumDims_; ++iDim) {
    Range.Begin_[iDim] = Begin[iDim];
    Range.End_[iDim] = End[iDim];
  }
  for (iDim = Range.NumDims_; iDim < OVK_MAX_DIMS; ++iDim) {
    Range.Begin_[iDim] = 0;
    Range.End_[iDim] = 1;
  }

  return Range;

}

static inline int ovkRangeDimension(const ovk_range *Range) {

  return Range->NumDims_;

}

static inline void ovkRangeBegin(const ovk_range *Range, int N, int *Begin) {

  int iDim;

  for (iDim = 0; iDim < ovk_min(OVK_MAX_DIMS, N); ++iDim) {
    Begin[iDim] = Range->Begin_[iDim];
  }

}

static inline int ovkRangeBeginDim(const ovk_range *Range, int iDim) {

  return Range->Begin_[iDim];

}

static inline const int *ovkRangeBeginData(const ovk_range *Range) {

  return Range->Begin_;

}

static inline void ovkRangeEnd(const ovk_range *Range, int N, int *End) {

  int iDim;

  for (iDim = 0; iDim < ovk_min(OVK_MAX_DIMS, N); ++iDim) {
    End[iDim] = Range->End_[iDim];
  }

}

static inline int ovkRangeEndDim(const ovk_range *Range, int iDim) {

  return Range->End_[iDim];

}

static inline const int *ovkRangeEndData(const ovk_range *Range) {

  return Range->End_;

}

static inline void ovkRangeSize(const ovk_range *Range, int N, int *Size) {

  int iDim;

  for (iDim = 0; iDim < ovk_min(OVK_MAX_DIMS, N); ++iDim) {
    Size[iDim] = Range->End_[iDim] - Range->Begin_[iDim];
  }

}

static inline int ovkRangeSizeDim(const ovk_range *Range, int iDim) {

  return Range->End_[iDim] - Range->Begin_[iDim];

}

static inline long long ovkRangeCount(const ovk_range *Range) {

  int iDim;

  long long Count = 1;

  for (iDim = 0; iDim < Range->NumDims_; ++iDim) {
    Count *= Range->End_[iDim] - Range->Begin_[iDim];
  }

  return Count;

}

static inline bool ovkRangeIsEmpty(const ovk_range *Range) {

  return
    Range->End_[0] <= Range->Begin_[0] ||
    Range->End_[1] <= Range->Begin_[1] ||
    Range->End_[2] <= Range->Begin_[2];

}

static inline bool ovkRangesAreEqual(const ovk_range *LeftRange, const ovk_range *RightRange) {

  return LeftRange->NumDims_ == RightRange->NumDims_ &&
    LeftRange->Begin_[0] == RightRange->Begin_[0] &&
    LeftRange->Begin_[1] == RightRange->Begin_[1] &&
    LeftRange->Begin_[2] == RightRange->Begin_[2] &&
    LeftRange->End_[0] == RightRange->End_[0] &&
    LeftRange->End_[1] == RightRange->End_[1] &&
    LeftRange->End_[2] == RightRange->End_[2];

}

static inline bool ovkRangeContains(const ovk_range *Range, const int *Point) {

  int iDim;

  for (iDim = 0; iDim < Range->NumDims_; ++iDim) {
    if (Point[iDim] < Range->Begin_[iDim] || Point[iDim] >= Range->End_[iDim]) return false;
  }

  return true;

}

static inline bool ovkRangeIncludes(const ovk_range *Range, const ovk_range *OtherRange) {

  return ovkRangeIsEmpty(OtherRange) ||
    (Range->Begin_[0] <= OtherRange->Begin_[0] && Range->End_[0] >= OtherRange->End_[0] &&
    Range->Begin_[1] <= OtherRange->Begin_[1] && Range->End_[1] >= OtherRange->End_[1] &&
    Range->Begin_[2] <= OtherRange->Begin_[2] && Range->End_[2] >= OtherRange->End_[2]);

}

static inline void ovkExtendRange(ovk_range *Range, const int *Point) {

  int iDim;

  if (ovkRangeIsEmpty(Range)) {
    for (iDim = 0; iDim < Range->NumDims_; ++iDim) {
      Range->Begin_[iDim] = Point[iDim];
      Range->End_[iDim] = Point[iDim]+1;
    }
  } else {
    for (iDim = 0; iDim < Range->NumDims_; ++iDim) {
      Range->Begin_[iDim] = ovk_min(Range->Begin_[iDim], Point[iDim]);
      Range->End_[iDim] = ovk_max(Range->End_[iDim], Point[iDim]+1);
    }
  }

}

static inline bool ovkRangesOverlap(const ovk_range *LeftRange, const ovk_range *RightRange) {

  return !ovkRangeIsEmpty(LeftRange) && !ovkRangeIsEmpty(RightRange) &&
    RightRange->End_[0] > LeftRange->Begin_[0] && LeftRange->End_[0] > RightRange->Begin_[0] &&
    RightRange->End_[1] > LeftRange->Begin_[1] && LeftRange->End_[1] > RightRange->Begin_[1] &&
    RightRange->End_[2] > LeftRange->Begin_[2] && LeftRange->End_[2] > RightRange->Begin_[2];

}

static inline ovk_range ovkUnionRanges(const ovk_range *LeftRange, const ovk_range *RightRange) {

  int iDim;

  ovk_range UnionRange;

  if (ovkRangeIsEmpty(LeftRange)) {
    UnionRange = *RightRange;
  } else if (ovkRangeIsEmpty(RightRange)) {
    UnionRange = *LeftRange;
  } else {
    UnionRange.NumDims_ = LeftRange->NumDims_;
    for (iDim = 0; iDim < UnionRange.NumDims_; ++iDim) {
      UnionRange.Begin_[iDim] = ovk_min(LeftRange->Begin_[iDim], RightRange->Begin_[iDim]);
      UnionRange.End_[iDim] = ovk_max(LeftRange->End_[iDim], RightRange->End_[iDim]);
    }
    for (iDim = UnionRange.NumDims_; iDim < OVK_MAX_DIMS; ++iDim) {
      UnionRange.Begin_[iDim] = 0;
      UnionRange.End_[iDim] = 1;
    }
  }

  return UnionRange;

}

static inline ovk_range ovkIntersectRanges(const ovk_range *LeftRange, const ovk_range *RightRange) {

  int iDim;

  ovk_range IntersectRange;

  IntersectRange.NumDims_ = LeftRange->NumDims_;
  for (iDim = 0; iDim < IntersectRange.NumDims_; ++iDim) {
    IntersectRange.Begin_[iDim] = ovk_max(LeftRange->Begin_[iDim], RightRange->Begin_[iDim]);
    IntersectRange.End_[iDim] = ovk_min(LeftRange->End_[iDim], RightRange->End_[iDim]);
  }
  for (iDim = IntersectRange.NumDims_; iDim < OVK_MAX_DIMS; ++iDim) {
    IntersectRange.Begin_[iDim] = 0;
    IntersectRange.End_[iDim] = 1;
  }

  return IntersectRange;

}

static inline void ovkClampToRange(int *Point, const ovk_range *Range) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Point[iDim] < Range->Begin_[iDim]) {
      Point[iDim] = Range->Begin_[iDim];
    } else if (Point[iDim] >= Range->End_[iDim]) {
      Point[iDim] = Range->End_[iDim]-1;
    }
  }

}

#ifdef __cplusplus
}
#endif
