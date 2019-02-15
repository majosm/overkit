// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline void DefaultCart(cart &Cart, int NumDims) {

  ovkDefaultCart(&Cart, NumDims);

}

template <typename SizeArrayType, typename PeriodicArrayType, OVK_FUNCDEF_REQUIRES(
  cart_internal::is_compatible_int_tuple_type<SizeArrayType>() &&
  cart_internal::is_compatible_bool_tuple_type<PeriodicArrayType>())>
  inline void SetCart(cart &Cart, int NumDims, const SizeArrayType &Size, const PeriodicArrayType
  &Periodic) {

  ovkSetCart(&Cart, NumDims, &Size[0], &Periodic[0]);

}

template <typename IntegerType, OVK_FUNCDEF_REQUIRES(std::is_integral<IntegerType>::value)> inline
  void CartCount(const cart &Cart, IntegerType &Count) {

  long long CountLongLong;
  ovkCartCount(&Cart, &CountLongLong);

  Count = IntegerType(CountLongLong);

}

template <typename PointArrayType, typename PeriodArrayType, OVK_FUNCDEF_REQUIRES(
  cart_internal::is_compatible_int_tuple_type<PointArrayType>() &&
  cart_internal::is_compatible_int_tuple_type<typename std::decay<PeriodArrayType>::type>() &&
  !std::is_const<typename std::decay<PeriodArrayType>::type>::value)> inline void CartFindPeriod(
  const cart &Cart, const PointArrayType &Point, PeriodArrayType &&Period) {

  ovkCartFindPeriod(&Cart, &Point[0], &Period[0]);

}

template <typename TupleArrayType, typename AdjustedTupleArrayType, OVK_FUNCDEF_REQUIRES(
  cart_internal::is_compatible_int_tuple_type<TupleArrayType>() &&
  cart_internal::is_compatible_int_tuple_type<typename std::decay<AdjustedTupleArrayType>::type>()
  && !std::is_const<typename std::decay<AdjustedTupleArrayType>::type>::value)> inline void
  CartPeriodicAdjust(const cart &Cart, const TupleArrayType &Tuple, AdjustedTupleArrayType
  &&AdjustedTuple) {

  ovkCartPeriodicAdjust(&Cart, &Tuple[0], &AdjustedTuple[0]);

}

inline void CartPointToCell(const cart &PointCart, cart &CellCart) {

  ovkCartPointToCell(&PointCart, &CellCart);

}

}

inline bool operator==(const ovk::cart &LeftCart, const ovk::cart &RightCart) {

  return ovkCartEquals(&LeftCart, &RightCart);

}

inline bool operator!=(const ovk::cart &LeftCart, const ovk::cart &RightCart) {

  return ovkCartEquals(&LeftCart, &RightCart);

}
