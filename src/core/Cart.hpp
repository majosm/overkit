// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CART_HPP_INCLUDED
#define OVK_CORE_CART_HPP_INCLUDED

#include <ovk/core/CartBase.h>
#include <ovk/core/Global.hpp>

namespace ovk {

using cart = ovk_cart;

inline void DefaultCart(cart &Cart, int NumDims) { ovkDefaultCart(&Cart, NumDims); }
inline void SetCart(cart &Cart, int NumDims, const int *Size, const bool *Periodic) { ovkSetCart(&Cart, NumDims, Size, Periodic); }
inline void CartCount(const cart &Cart, long long &Count) { ovkCartCount(&Cart, &Count); }
inline bool CartEquals(const cart &Cart, const cart &OtherCart) { return ovkCartEquals(&Cart, &OtherCart); }
inline void CartFindPeriod(const cart &Cart, const int *Point, int *Period) { ovkCartFindPeriod(&Cart, Point, Period); }
inline void CartPeriodicAdjust(const cart &Cart, const int *Tuple, int *AdjustedTuple) { ovkCartPeriodicAdjust(&Cart, Tuple, AdjustedTuple); }
inline void CartPointToCell(const cart &PointCart, cart &CellCart) { ovkCartPointToCell(&PointCart, &CellCart); }

}

#endif
