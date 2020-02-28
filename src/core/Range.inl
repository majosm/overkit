// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
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

inline range ExtendRange(const range &Range, const tuple<int> &Point) {

  range Result;

  if (!Range.Empty()) {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Result.Begin(iDim) = Min(Range.Begin(iDim), Point(iDim));
      Result.End(iDim) = Max(Range.End(iDim), Point(iDim)+1);
    }
  } else {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Result.Begin(iDim) = Point(iDim);
      Result.End(iDim) = Point(iDim)+1;
    }
  }

  return Result;

}

inline bool RangesOverlap(const range &LeftRange, const range &RightRange) {

  return !LeftRange.Empty() && !RightRange.Empty() &&
    RightRange.End(0) > LeftRange.Begin(0) && LeftRange.End(0) > RightRange.Begin(0) &&
    RightRange.End(1) > LeftRange.Begin(1) && LeftRange.End(1) > RightRange.Begin(1) &&
    RightRange.End(2) > LeftRange.Begin(2) && LeftRange.End(2) > RightRange.Begin(2);

}

inline range UnionRanges(const range &LeftRange, const range &RightRange) {

  range Result;

  if (LeftRange.Empty()) {
    Result = RightRange;
  } else if (RightRange.Empty()) {
    Result = LeftRange;
  } else {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Result.Begin(iDim) = Min(LeftRange.Begin(iDim), RightRange.Begin(iDim));
      Result.End(iDim) = Max(LeftRange.End(iDim), RightRange.End(iDim));
    }
  }

  return Result;

}

inline range IntersectRanges(const range &LeftRange, const range &RightRange) {

  range Result;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Result.Begin(iDim) = Max(LeftRange.Begin(iDim), RightRange.Begin(iDim));
    Result.End(iDim) = Min(LeftRange.End(iDim), RightRange.End(iDim));
  }

  return Result;

}

inline tuple<int> ClampToRange(const range &Range, const tuple<int> &Point) {

  tuple<int> Result;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    if (Point(iDim) < Range.Begin(iDim)) {
      Result(iDim) = Range.Begin(iDim);
    } else if (Point(iDim) >= Range.End(iDim)) {
      Result(iDim) = Range.End(iDim)-1;
    } else {
      Result(iDim) = Point(iDim);
    }
  }

  return Result;

}

}
