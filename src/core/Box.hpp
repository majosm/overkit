// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_BOX_HPP_INCLUDED
#define OVK_CORE_BOX_HPP_INCLUDED

#include <ovk/core/BoxBase.h>
#include <ovk/core/Global.hpp>

namespace ovk {

using box = ovk_box;

inline void DefaultBox(box &Box, int NumDims) { ovkDefaultBox(&Box, NumDims); }
inline void ovkSetBox(box &Box, int NumDims, const double *Begin, const double *End) { ovkSetBox(&Box, NumDims, Begin, End); }
inline void ovkBoxSize(const box &Box, double *Size) { ovkBoxSize(&Box, Size); }
inline void ovkBoxVolume(const box &Box, double *Volume) { ovkBoxVolume(&Box, Volume); }
inline bool ovkBoxIsEmpty(const box &Box) { return ovkBoxIsEmpty(&Box); }
inline bool ovkBoxEquals(const box &LeftBox, const box &RightBox) { return ovkBoxEquals(&LeftBox, &RightBox); }
inline bool ovkBoxContains(const box &Box, const double *Point) { return ovkBoxContains(&Box, Point); }
inline bool ovkBoxOverlaps(const box &LeftBox, const box &RightBox) { return ovkBoxOverlaps(&LeftBox, &RightBox); }
inline void ovkBoxUnion(const box &LeftBox, const box &RightBox, box &UnionBox) { ovkBoxUnion(&LeftBox, &RightBox, &UnionBox); }
inline void ovkBoxIntersect(const box &LeftBox, const box &RightBox, box &IntersectBox) { ovkBoxIntersect(&LeftBox, &RightBox, &IntersectBox); }
inline void ovkBoxClamp(const box &Box, double *Point) { ovkBoxClamp(&Box, Point); }
inline void ovkBoxMove(const box &Box, const double *Amount, box &MoveBox) { ovkBoxMove(&Box, Amount, &MoveBox); }
inline void ovkBoxGrow(const box &Box, const double *Amount, box &GrowBox) { ovkBoxGrow(&Box, Amount, &GrowBox); }
inline void ovkBoxScale(const box &Box, const double *Factor, box &ScaleBox) { ovkBoxScale(&Box, Factor, &ScaleBox); }
inline void ovkBoxExtend(const box &Box, const double *Point, box &ExtendBox) { ovkBoxExtend(&Box, Point, &ExtendBox); }

}

#endif
