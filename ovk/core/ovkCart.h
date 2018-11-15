// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CART_INCLUDED
#define OVK_CORE_PUBLIC_CART_INCLUDED

#include <ovk/core/ovkGlobal.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int nd;
  int size[OVK_MAX_DIMS];
  bool periodic[OVK_MAX_DIMS];
} ovk_cart;

static inline void ovkDefaultCart(ovk_cart *Cart, int NumDims);
static inline void ovkSetCart(ovk_cart *Cart, int NumDims, const int *Size, const bool *Periodic);
static inline void ovkCartCount(const ovk_cart *Cart, long long *Count);
static inline bool ovkCartEquals(const ovk_cart *Cart, const ovk_cart *OtherCart);
static inline void ovkCartFindPeriod(const ovk_cart *Cart, const int *Point, int *Period);
static inline void ovkCartPeriodicAdjust(const ovk_cart *Cart, const int *Tuple,
  int *AdjustedTuple);
static inline void ovkCartPointToCell(const ovk_cart *PointCart, ovk_cart *CellCart);

#ifdef __cplusplus
}
#endif

#endif
