// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline void DefaultCart(cart &Cart, int NumDims) {

  ovkDefaultCart(&Cart, NumDims);

}

template <typename IntArrayType, typename BoolArrayType> inline void SetCart(cart &Cart,
  int NumDims, const IntArrayType &Size, const BoolArrayType &Periodic) {

  ovkSetCart(&Cart, NumDims, &Size[0], &Periodic[0]);

}

template <typename IntegerType> inline void CartCount(const cart &Cart, IntegerType &Count) {

  long long CountLongLong;
  ovkCartCount(&Cart, &CountLongLong);

  Count = IntegerType(CountLongLong);

}

template <typename IntArrayType1, typename IntArrayType2> inline void CartFindPeriod(const cart
  &Cart, const IntArrayType1 &Point, IntArrayType2 &Period) {

  ovkCartFindPeriod(&Cart, &Point[0], &Period[0]);

}

template <typename IntArrayType1, typename IntArrayType2> inline void CartPeriodicAdjust(const cart
  &Cart, const IntArrayType1 &Tuple, IntArrayType2 &AdjustedTuple) {

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
