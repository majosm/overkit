// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CART_INL_INCLUDED
#define OVK_CORE_PUBLIC_CART_INL_INCLUDED

#include <ovk/core/ovkCart.h>

#include <ovk/core/ovkGlobal.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline void ovkDefaultCart(ovk_cart *Cart, int NumDims) {

  int iDim;

  Cart->nd = NumDims;
  for (iDim = 0; iDim < Cart->nd; ++iDim) {
    Cart->size[iDim] = 0;
    Cart->periodic[iDim] = false;
  }
  for (iDim = Cart->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    Cart->size[iDim] = 1;
    Cart->periodic[iDim] = false;
  }

}

static inline void ovkSetCart(ovk_cart *Cart, int NumDims, const int *Size, const bool *Periodic) {

  int iDim;

  Cart->nd = NumDims;
  for (iDim = 0; iDim < Cart->nd; ++iDim) {
    Cart->size[iDim] = Size[iDim];
    Cart->periodic[iDim] = Periodic[iDim];
  }
  for (iDim = Cart->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    Cart->size[iDim] = 1;
    Cart->periodic[iDim] = false;
  }

}

static inline void ovkCartCount(const ovk_cart *Cart, size_t *Count) {

  int iDim;

  *Count = 1;

  for (iDim = 0; iDim < Cart->nd; ++iDim) {
    *Count *= (size_t)Cart->size[iDim];
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

static inline void ovkCartFindPeriod(const ovk_cart *Cart, const int *Point, int *Period) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Cart->periodic[iDim]) {
      Period[iDim] = Point[iDim]/Cart->size[iDim] - (int)(Point[iDim] % Cart->size[iDim] < 0);
    } else {
      Period[iDim] = 0;
    }
  }

}

static inline void ovkCartPeriodicAdjust(const ovk_cart *Cart, const int *Tuple,
  int *AdjustedTuple) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Cart->periodic[iDim]) {
      AdjustedTuple[iDim] = Tuple[iDim] % Cart->size[iDim];
    } else {
      AdjustedTuple[iDim] = Tuple[iDim];
    }
  }

}

static inline void ovkCartPointToCell(const ovk_cart *PointCart, ovk_cart *CellCart) {

  int iDim;

  CellCart->nd = PointCart->nd;
  for (iDim = 0; iDim < PointCart->nd; ++iDim) {
    CellCart->size[iDim] = PointCart->periodic[iDim] ? PointCart->size[iDim] :
      PointCart->size[iDim]-1;
    CellCart->periodic[iDim] = PointCart->periodic[iDim];
  }
  for (iDim = PointCart->nd; iDim < OVK_MAX_DIMS; ++iDim) {
    CellCart->size[iDim] = 1;
    CellCart->periodic[iDim] = false;
  }

}

#ifdef __cplusplus
}
#endif

#endif
