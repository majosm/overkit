// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_BOX_HPP_INCLUDED
#define OVK_CORE_BOX_HPP_INCLUDED

#include <ovk/core/BoxBase.h>
#include <ovk/core/Global.hpp>

namespace ovk {

using box = ovk_box;

inline void DefaultBox(box &Box, int NumDims);
template <typename DoubleArrayType1, typename DoubleArrayType2> inline void ovkSetBox(box &Box,
  int NumDims, const DoubleArrayType1 &Begin, const DoubleArrayType2 &End);
template <typename DoubleArrayType> inline void ovkBoxSize(const box &Box, DoubleArrayType &Size);
inline void ovkBoxVolume(const box &Box, double &Volume);
inline bool ovkBoxIsEmpty(const box &Box);
template <typename DoubleArrayType> inline bool ovkBoxContains(const box &Box, const DoubleArrayType
  &Point);
inline bool ovkBoxOverlaps(const box &LeftBox, const box &RightBox);
inline void ovkBoxUnion(const box &LeftBox, const box &RightBox, box &UnionBox);
inline void ovkBoxIntersect(const box &LeftBox, const box &RightBox, box &IntersectBox);
template <typename DoubleArrayType> inline void ovkBoxClamp(const box &Box, DoubleArrayType &Point);
template <typename DoubleArrayType> inline void ovkBoxMove(const box &Box, const DoubleArrayType
  &Amount, box &MoveBox);
template <typename DoubleArrayType> inline void ovkBoxGrow(const box &Box, const DoubleArrayType
  &Amount, box &GrowBox);
template <typename DoubleArrayType> inline void ovkBoxScale(const box &Box, const DoubleArrayType
  &Factor, box &ScaleBox);
template <typename DoubleArrayType> inline void ovkBoxExtend(const box &Box, const DoubleArrayType
  &Point, box &ExtendBox);

}

// Can't put these in the ovk namespace (ADL won't work since ovk_box isn't defined there)
inline bool operator==(const ovk::box &LeftBox, const ovk::box &RightBox);
inline bool operator!=(const ovk::box &LeftBox, const ovk::box &RightBox);

#include <ovk/core/Box.inl>

#endif
