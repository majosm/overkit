// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline void DefaultBox(box &Box, int NumDims) {

  ovkDefaultBox(&Box, NumDims);

}

template <typename DoubleArrayType1, typename DoubleArrayType2> inline void ovkSetBox(box &Box,
  int NumDims, const DoubleArrayType1 &Begin, const DoubleArrayType2 &End) {

  ovkSetBox(&Box, NumDims, &Begin[0], &End[0]);

}

template <typename DoubleArrayType> inline void ovkBoxSize(const box &Box, DoubleArrayType &Size) {

  ovkBoxSize(&Box, &Size[0]);

}

inline void ovkBoxVolume(const box &Box, double &Volume) {

  ovkBoxVolume(&Box, &Volume);

}

inline bool ovkBoxIsEmpty(const box &Box) {

  return ovkBoxIsEmpty(&Box);

}

template <typename DoubleArrayType> inline bool ovkBoxContains(const box &Box, const DoubleArrayType
  &Point) {

  return ovkBoxContains(&Box, &Point[0]);

}

inline bool ovkBoxOverlaps(const box &LeftBox, const box &RightBox) {

  return ovkBoxOverlaps(&LeftBox, &RightBox);

}

inline void ovkBoxUnion(const box &LeftBox, const box &RightBox, box &UnionBox) {

  ovkBoxUnion(&LeftBox, &RightBox, &UnionBox);

}

inline void ovkBoxIntersect(const box &LeftBox, const box &RightBox, box &IntersectBox) {

  ovkBoxIntersect(&LeftBox, &RightBox, &IntersectBox);

}

template <typename DoubleArrayType> inline void ovkBoxClamp(const box &Box, DoubleArrayType &Point) {

  ovkBoxClamp(&Box, &Point[0]);

}

template <typename DoubleArrayType> inline void ovkBoxMove(const box &Box, const DoubleArrayType
  &Amount, box &MoveBox) {

  ovkBoxMove(&Box, &Amount[0], &MoveBox);

}

template <typename DoubleArrayType> inline void ovkBoxGrow(const box &Box, const DoubleArrayType
  &Amount, box &GrowBox) {

  ovkBoxGrow(&Box, &Amount[0], &GrowBox);

}

template <typename DoubleArrayType> inline void ovkBoxScale(const box &Box, const DoubleArrayType
  &Factor, box &ScaleBox) {

  ovkBoxScale(&Box, &Factor[0], &ScaleBox);

}

template <typename DoubleArrayType> inline void ovkBoxExtend(const box &Box, const DoubleArrayType
  &Point, box &ExtendBox) {

  ovkBoxExtend(&Box, &Point[0], &ExtendBox);

}

}

inline bool operator==(const ovk::box &LeftBox, const ovk::box &RightBox) {

  return ovkBoxEquals(&LeftBox, &RightBox);

}

inline bool operator!=(const ovk::box &LeftBox, const ovk::box &RightBox) {

  return !ovkBoxEquals(&LeftBox, &RightBox);

}
