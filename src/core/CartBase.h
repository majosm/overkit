// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CART_BASE_H_INCLUDED
#define OVK_CORE_CART_BASE_H_INCLUDED

#include <ovk/core/ConstantsBase.h>
#include <ovk/core/GlobalBase.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int NumDims;
  int Size[OVK_MAX_DIMS];
  bool Periodic[OVK_MAX_DIMS];
} ovk_cart;

static inline void ovkDefaultCart(ovk_cart *Cart, int NumDims);
static inline void ovkSetCart(ovk_cart *Cart, int NumDims, const int *Size, const bool *Periodic);
static inline void ovkCartCount(const ovk_cart *Cart, long long *Count);
static inline bool ovkCartEquals(const ovk_cart *Cart, const ovk_cart *OtherCart);
static inline void ovkCartFindPeriod(const ovk_cart *Cart, const int *Point, int *Period);
static inline void ovkCartPeriodicAdjust(const ovk_cart *Cart, const int *Tuple, int *AdjustedTuple);
static inline void ovkCartPointToCell(const ovk_cart *PointCart, ovk_cart *CellCart);

#ifdef __cplusplus
}
#endif

#include <ovk/core/CartBase.inl>

#endif
