// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline range MakeEmptyRange(int NumDims) {

  range Range;

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    Range.Begin(iDim) = 0;
    Range.End(iDim) = 0;
  }
  for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
    Range.Begin(iDim) = 0;
    Range.End(iDim) = 1;
  }

  return Range;

}

inline void ExtendRange(range &Range, const tuple<int> &Tuple) {

  if (!Range.Empty()) {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Range.Begin(iDim) = Min(Range.Begin(iDim), Tuple(iDim));
      Range.End(iDim) = Max(Range.End(iDim), Tuple(iDim)+1);
    }
  } else {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Range.Begin(iDim) = Tuple(iDim);
      Range.End(iDim) = Tuple(iDim)+1;
    }
  }

}

inline bool RangesOverlap(const range &LeftRange, const range &RightRange) {

  return !LeftRange.Empty() && !RightRange.Empty() &&
    RightRange.End(0) > LeftRange.Begin(0) && LeftRange.End(0) > RightRange.Begin(0) &&
    RightRange.End(1) > LeftRange.Begin(1) && LeftRange.End(1) > RightRange.Begin(1) &&
    RightRange.End(2) > LeftRange.Begin(2) && LeftRange.End(2) > RightRange.Begin(2);

}

inline range UnionRanges(const range &LeftRange, const range &RightRange) {

  range UnionRange;

  if (LeftRange.Empty()) {
    UnionRange = RightRange;
  } else if (RightRange.Empty()) {
    UnionRange = LeftRange;
  } else {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      UnionRange.Begin(iDim) = Min(LeftRange.Begin(iDim), RightRange.Begin(iDim));
      UnionRange.End(iDim) = Max(LeftRange.End(iDim), RightRange.End(iDim));
    }
  }

  return UnionRange;

}

inline range IntersectRanges(const range &LeftRange, const range &RightRange) {

  range IntersectRange;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    IntersectRange.Begin(iDim) = Max(LeftRange.Begin(iDim), RightRange.Begin(iDim));
    IntersectRange.End(iDim) = Min(LeftRange.End(iDim), RightRange.End(iDim));
  }

  return IntersectRange;

}

inline void ClampToRange(const range &Range, tuple<int> &Tuple) {

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    if (Tuple(iDim) < Range.Begin(iDim)) {
      Tuple(iDim) = Range.Begin(iDim);
    } else if (Tuple(iDim) >= Range.End(iDim)) {
      Tuple(iDim) = Range.End(iDim)-1;
    }
  }

}

}
