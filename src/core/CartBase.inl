// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifdef __cplusplus
extern "C" {
#endif

static inline void ovkDefaultCart(ovk_cart *Cart, int NumDims) {

  int iDim;

  Cart->NumDims = NumDims;
  for (iDim = 0; iDim < Cart->NumDims; ++iDim) {
    Cart->Size[iDim] = 0;
    Cart->Periodic[iDim] = false;
  }
  for (iDim = Cart->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    Cart->Size[iDim] = 1;
    Cart->Periodic[iDim] = false;
  }

}

static inline void ovkSetCart(ovk_cart *Cart, int NumDims, const int *Size, const bool *Periodic) {

  int iDim;

  Cart->NumDims = NumDims;
  for (iDim = 0; iDim < Cart->NumDims; ++iDim) {
    Cart->Size[iDim] = Size[iDim];
    Cart->Periodic[iDim] = Periodic[iDim];
  }
  for (iDim = Cart->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    Cart->Size[iDim] = 1;
    Cart->Periodic[iDim] = false;
  }

}

static inline bool ovkCartEquals(const ovk_cart *LeftCart, const ovk_cart *RightCart) {

  return LeftCart->NumDims == RightCart->NumDims &&
    LeftCart->Size[0] == RightCart->Size[0] &&
    LeftCart->Size[1] == RightCart->Size[1] &&
    LeftCart->Size[2] == RightCart->Size[2] &&
    LeftCart->Periodic[0] == RightCart->Periodic[0] &&
    LeftCart->Periodic[1] == RightCart->Periodic[1] &&
    LeftCart->Periodic[2] == RightCart->Periodic[2];

}

static inline void ovkCartCount(const ovk_cart *Cart, long long *Count) {

  int iDim;

  *Count = 1;

  for (iDim = 0; iDim < Cart->NumDims; ++iDim) {
    *Count *= Cart->Size[iDim];
  }

}

static inline void ovkCartFindPeriod(const ovk_cart *Cart, const int *Point, int *Period) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Cart->Periodic[iDim]) {
      Period[iDim] = Point[iDim]/Cart->Size[iDim] - (int)(Point[iDim] % Cart->Size[iDim] < 0);
    } else {
      Period[iDim] = 0;
    }
  }

}

static inline void ovkCartPeriodicAdjust(const ovk_cart *Cart, const int *Tuple, int *AdjustedTuple) {

  int iDim;

  for (iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
    if (Cart->Periodic[iDim]) {
      int Mod = Tuple[iDim] % Cart->Size[iDim];
      AdjustedTuple[iDim] = Mod + Cart->Size[iDim] * (Mod < 0);
    } else {
      AdjustedTuple[iDim] = Tuple[iDim];
    }
  }

}

static inline void ovkCartPointToCell(const ovk_cart *PointCart, ovk_cart *CellCart) {

  int iDim;

  CellCart->NumDims = PointCart->NumDims;
  for (iDim = 0; iDim < PointCart->NumDims; ++iDim) {
    CellCart->Size[iDim] = PointCart->Periodic[iDim] ? PointCart->Size[iDim] :
      PointCart->Size[iDim]-1;
    CellCart->Periodic[iDim] = PointCart->Periodic[iDim];
  }
  for (iDim = PointCart->NumDims; iDim < OVK_MAX_DIMS; ++iDim) {
    CellCart->Size[iDim] = 1;
    CellCart->Periodic[iDim] = false;
  }

}

#ifdef __cplusplus
}
#endif
