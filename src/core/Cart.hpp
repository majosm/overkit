// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CART_HPP_INCLUDED
#define OVK_CORE_CART_HPP_INCLUDED

#include <ovk/core/CartBase.h>
#include <ovk/core/Global.hpp>

namespace ovk {

using cart = ovk_cart;

inline void DefaultCart(cart &Cart, int NumDims);
template <typename IntArrayType, typename BoolArrayType> inline void SetCart(cart &Cart,
  int NumDims, const IntArrayType &Size, const BoolArrayType &Periodic);
template <typename IntegerType> inline void CartCount(const cart &Cart, IntegerType &Count);
template <typename IntArrayType1, typename IntArrayType2> inline void CartFindPeriod(const cart
  &Cart, const IntArrayType1 &Point, IntArrayType2 &Period);
template <typename IntArrayType1, typename IntArrayType2> inline void CartPeriodicAdjust(const cart
  &Cart, const IntArrayType1 &Tuple, IntArrayType2 &AdjustedTuple);
inline void CartPointToCell(const cart &PointCart, cart &CellCart);

}

// Can't put these in the ovk namespace (ADL won't work since ovk_cart isn't defined there)
inline bool operator==(const ovk::cart &LeftCart, const ovk::cart &RightCart);
inline bool operator!=(const ovk::cart &LeftCart, const ovk::cart &RightCart);

#include <ovk/core/Cart.inl>

#endif
