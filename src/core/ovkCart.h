// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CART_INCLUDED
#define OVK_CORE_PUBLIC_CART_INCLUDED

#include <ovkGlobal.h>

struct ovk_cart;
typedef struct ovk_cart ovk_cart;

static inline void ovkCartDefault(ovk_cart *Cart, int NumDims);
static inline void ovkCartCount(const ovk_cart *Cart, size_t *Count);
static inline bool ovkCartEquals(const ovk_cart *Cart, const ovk_cart *OtherCart);
static inline void ovkCartPeriodicAdjust(const ovk_cart *Cart, const int *Tuple,
  int *AdjustedTuple);
static inline void ovkCartPointToCell(const ovk_cart *PointCart, ovk_cart *CellCart);

#endif
