// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CART_HPP_INCLUDED
#define OVK_CORE_CART_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/CartBase.h>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>

#include <type_traits>

namespace ovk {

using cart = ovk_cart;

namespace cart_internal {

template <typename TupleType> static constexpr bool is_compatible_int_tuple_type() {
  return (core::IsArray<TupleType>() && std::is_same<core::array_value_type<TupleType>, int>::value
    && core::ArrayRank<TupleType>() == 1) || std::is_same<TupleType, int *>::value;
}

template <typename TupleType> static constexpr bool is_compatible_bool_tuple_type() {
  return (core::IsArray<TupleType>() && std::is_same<core::array_value_type<TupleType>, bool>::value
    && core::ArrayRank<TupleType>() == 1) || std::is_same<TupleType, bool *>::value;
}

}

inline void DefaultCart(cart &Cart, int NumDims);

template <typename SizeTupleType, typename PeriodicTupleType, OVK_FUNCDECL_REQUIRES(
  cart_internal::is_compatible_int_tuple_type<SizeTupleType>() &&
  cart_internal::is_compatible_bool_tuple_type<PeriodicTupleType>())>
  inline void SetCart(cart &Cart, int NumDims, const SizeTupleType &Size, const PeriodicTupleType
  &Periodic);

template <typename IntegerType, OVK_FUNCDECL_REQUIRES(std::is_integral<IntegerType>::value)> inline
  void CartCount(const cart &Cart, IntegerType &Count);

// && in order to bind to temporary pointers (i.e., const pointer to non-const data); not sure if
// there is a nicer way to do this
template <typename PointTupleType, typename PeriodTupleType, OVK_FUNCDECL_REQUIRES(
  cart_internal::is_compatible_int_tuple_type<PointTupleType>() &&
  cart_internal::is_compatible_int_tuple_type<typename std::decay<PeriodTupleType>::type>() &&
  !std::is_const<typename std::decay<PeriodTupleType>::type>::value)> inline void
  CartFindPeriod(const cart &Cart, const PointTupleType &Point, PeriodTupleType &&Period);

// && in order to bind to temporary pointers (i.e., const pointer to non-const data); not sure if
// there is a nicer way to do this
template <typename TupleTupleType, typename AdjustedTupleTupleType, OVK_FUNCDECL_REQUIRES(
  cart_internal::is_compatible_int_tuple_type<TupleTupleType>() &&
  cart_internal::is_compatible_int_tuple_type<typename std::decay<AdjustedTupleTupleType>::type>()
  && !std::is_const<typename std::decay<AdjustedTupleTupleType>::type>::value)> inline void
  CartPeriodicAdjust(const cart &Cart, const TupleTupleType &Tuple, AdjustedTupleTupleType
  &&AdjustedTuple);

inline void CartPointToCell(const cart &PointCart, cart &CellCart);

}

// Can't put these in the ovk namespace (ADL won't work since ovk_cart isn't defined there)
inline bool operator==(const ovk::cart &LeftCart, const ovk::cart &RightCart);
inline bool operator!=(const ovk::cart &LeftCart, const ovk::cart &RightCart);

#include <ovk/core/Cart.inl>

#endif
