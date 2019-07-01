// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline box MakeEmptyBox(int NumDims) {

  box Box;

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    Box.Begin(iDim) = 0.;
    Box.End(iDim) = -1.;
  }
  for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
    Box.Begin(iDim) = 0.;
    Box.End(iDim) = 0.;
  }

  return Box;

}

inline box ExtendBox(const box &Box, const tuple<double> &Point) {

  box Result;

  if (!Box.Empty()) {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Result.Begin(iDim) = Min(Box.Begin(iDim), Point(iDim));
      Result.End(iDim) = Max(Box.End(iDim), Point(iDim));
    }
  } else {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Result.Begin(iDim) = Point(iDim);
      Result.End(iDim) = Point(iDim);
    }
  }

  return Result;

}

inline bool BoxesOverlap(const box &LeftBox, const box &RightBox) {

  return !LeftBox.Empty() && !RightBox.Empty() &&
    RightBox.End(0) >= LeftBox.Begin(0) && LeftBox.End(0) >= RightBox.Begin(0) &&
    RightBox.End(1) >= LeftBox.Begin(1) && LeftBox.End(1) >= RightBox.Begin(1) &&
    RightBox.End(2) >= LeftBox.Begin(2) && LeftBox.End(2) >= RightBox.Begin(2);

}

inline box UnionBoxes(const box &LeftBox, const box &RightBox) {

  box Result;

  if (LeftBox.Empty()) {
    Result = RightBox;
  } else if (RightBox.Empty()) {
    Result = LeftBox;
  } else {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Result.Begin(iDim) = Min(LeftBox.Begin(iDim), RightBox.Begin(iDim));
      Result.End(iDim) = Max(LeftBox.End(iDim), RightBox.End(iDim));
    }
  }

  return Result;

}

inline box IntersectBoxes(const box &LeftBox, const box &RightBox) {

  box Result;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Result.Begin(iDim) = Max(LeftBox.Begin(iDim), RightBox.Begin(iDim));
    Result.End(iDim) = Min(LeftBox.End(iDim), RightBox.End(iDim));
  }

  return Result;

}

inline tuple<double> ClampToBox(const box &Box, const tuple<double> &Point) {

  tuple<double> Result;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    if (Point(iDim) < Box.Begin(iDim)) {
      Result(iDim) = Box.Begin(iDim);
    } else if (Point(iDim) > Box.End(iDim)) {
      Result(iDim) = Box.End(iDim);
    } else {
      Result(iDim) = Point(iDim);
    }
  }

  return Result;

}

inline box MoveBox(const box &Box, const tuple<double> &Amount) {

  box Result;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Result.Begin(iDim) = Box.Begin(iDim) + Amount(iDim);
    Result.End(iDim) = Box.End(iDim) + Amount(iDim);
  }

  return Result;

}

inline box GrowBox(const box &Box, const tuple<double> &Amount) {

  box Result;

  if (Box.Empty()) {
    Result = Box;
  } else {
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Result.Begin(iDim) = Box.Begin(iDim) - Amount(iDim);
      Result.End(iDim) = Box.End(iDim) + Amount(iDim);
    }
  }

  return Result;

}

inline box ScaleBox(const box &Box, const tuple<double> &Factor) {

  box Result;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    double Center = 0.5 * (Box.Begin(iDim) + Box.End(iDim));
    double HalfSize = 0.5 * Factor(iDim) * (Box.End(iDim) - Box.Begin(iDim));
    Result.Begin(iDim) = Center - HalfSize;
    Result.End(iDim) = Center + HalfSize;
  }

  return Result;

}

}
