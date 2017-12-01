// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CART_INL_INCLUDED
#define OVK_CORE_PUBLIC_CART_INL_INCLUDED

#include <ovkCart.h>

#include <ovkGlobal.h>

struct ovk_cart {
  int nd;
  int size[OVK_MAX_DIMS];
  bool periodic[OVK_MAX_DIMS];
};

static inline void ovkCartDefault(ovk_cart *Cart, int NumDims) {

  Cart->nd = NumDims;
  for (int d = 0; d < Cart->nd; ++d) {
    Cart->size[d] = 0;
    Cart->periodic[d] = false;
  }
  for (int d = Cart->nd; d < OVK_MAX_DIMS; ++d) {
    Cart->size[d] = 1;
    Cart->periodic[d] = false;
  }

}

static inline void ovkCartCount(const ovk_cart *Cart, size_t *Count) {

  *Count = 1;

  for (int d = 0; d < Cart->nd; ++d) {
    *Count *= (size_t)Cart->size[d];
  }

}

static inline bool ovkCartEquals(const ovk_cart *Cart, const ovk_cart *OtherCart) {

  return Cart->nd == OtherCart->nd &&
    Cart->size[0] == OtherCart->size[0] &&
    Cart->size[1] == OtherCart->size[1] &&
    Cart->size[2] == OtherCart->size[2] &&
    Cart->periodic[0] == OtherCart->periodic[0] &&
    Cart->periodic[1] == OtherCart->periodic[1] &&
    Cart->periodic[2] == OtherCart->periodic[2];

}

static inline void ovkCartPeriodicAdjust(const ovk_cart *Cart, const int *Tuple,
  int *AdjustedTuple) {

  for (int d = 0; d < Cart->nd; ++d) {
    if (Cart->periodic[d]) {
      AdjustedTuple[d] = Tuple[d] % Cart->size[d];
    } else {
      AdjustedTuple[d] = Tuple[d];
    }
  }

}

static inline void ovkCartPointToCell(const ovk_cart *PointCart, ovk_cart *CellCart) {

  CellCart->nd = PointCart->nd;
  for (int d = 0; d < PointCart->nd; ++d) {
    CellCart->size[d] = PointCart->periodic[d] ? PointCart->size[d] : PointCart->size[d]-1;
    CellCart->periodic[d] = PointCart->periodic[d];
  }
  for (int d = PointCart->nd; d < OVK_MAX_DIMS; ++d) {
    CellCart->size[d] = 1;
    CellCart->periodic[d] = false;
  }

}

#endif
