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

static inline void ovkCartCount(const ovk_cart *Cart, long long *Count) {

  int iDim;

  *Count = 1;

  for (iDim = 0; iDim < Cart->NumDims; ++iDim) {
    *Count *= Cart->Size[iDim];
  }

}

static inline bool ovkCartEquals(const ovk_cart *Cart, const ovk_cart *OtherCart) {

  return Cart->NumDims == OtherCart->NumDims &&
    Cart->Size[0] == OtherCart->Size[0] &&
    Cart->Size[1] == OtherCart->Size[1] &&
    Cart->Size[2] == OtherCart->Size[2] &&
    Cart->Periodic[0] == OtherCart->Periodic[0] &&
    Cart->Periodic[1] == OtherCart->Periodic[1] &&
    Cart->Periodic[2] == OtherCart->Periodic[2];

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
