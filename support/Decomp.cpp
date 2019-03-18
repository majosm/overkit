// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "support/Decomp.hpp"

#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

namespace support {

ovk::range CartesianDecomp(const ovk::cart &Cart, ovk::core::comm_view Comm, const ovk::tuple<int>
  &CartDims) {

  int NumDims = Cart.Dimension();

  ovk::indexer<int, int, OVK_MAX_DIMS> Indexer(CartDims);
  ovk::tuple<int> CartCoords = Indexer.ToTuple(Comm.Rank());

  ovk::range LocalRange = ovk::MakeEmptyRange(NumDims);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    int NumPerRank = Cart.Range().Size(iDim)/CartDims[iDim];
    int Remainder = Cart.Range().Size(iDim) - CartDims[iDim]*NumPerRank;
    int Begin = Cart.Range().Begin(iDim);
    int Coord = CartCoords[iDim];
    LocalRange.Begin(iDim) = Begin + NumPerRank*Coord + std::min(Remainder, Coord);
    LocalRange.End(iDim) = Begin + NumPerRank*(Coord+1) + std::min(Remainder, Coord+1);
  }

  return LocalRange;

}

ovk::range TriangularDecomp(const ovk::cart &Cart, ovk::core::comm_view Comm) {

  int NumDims = Cart.Dimension();

  auto CountProcs = [NumDims](int NumPlanes) -> int {
    if (NumDims == 2) {
      // Sum of j+1 from 0 to NumPlanes-1
      return NumPlanes*(NumPlanes+1)/2;
    } else {
      // Sum of (j+1)^2 from 0 to NumPlanes-1
      return NumPlanes*(NumPlanes+1)*(2*NumPlanes+1)/6;
    }
  };

  int NumPlanes = 0;
  while (CountProcs(NumPlanes+1) <= Comm.Size()) {
    ++NumPlanes;
  }

  int iPlane = 0;
  while (CountProcs(iPlane+1) <= Comm.Rank()) {
    ++iPlane;
  }

  int iProcInPlane = Comm.Rank() - CountProcs(iPlane);

  ovk::tuple<int> TriangleDims = ovk::MakeUniformTuple<int>(1);
  for (int iDim = 0; iDim < NumDims-1; ++iDim) {
    TriangleDims[iDim] = iPlane+1;
  }
  TriangleDims[NumDims-1] = NumPlanes;

  ovk::tuple<int> PlaneDims = ovk::MakeUniformTuple<int>(1);
  for (int iDim = 0; iDim < NumDims-1; ++iDim) {
    PlaneDims[iDim] = iPlane+1;
  }
  ovk::indexer<int, int, OVK_MAX_DIMS> Indexer(PlaneDims);
  ovk::tuple<int> TriangleCoords = Indexer.ToTuple(iProcInPlane);
  TriangleCoords[NumDims-1] = iPlane;

  ovk::range LocalRange = ovk::MakeEmptyRange(NumDims);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    int NumPerRank = Cart.Range().Size(iDim)/TriangleDims[iDim];
    int Remainder = Cart.Range().Size(iDim) - TriangleDims[iDim]*NumPerRank;
    int Begin = Cart.Range().Begin(iDim);
    int Coord = TriangleCoords[iDim];
    LocalRange.Begin(iDim) = Begin + NumPerRank*Coord + std::min(Remainder, Coord);
    LocalRange.End(iDim) = Begin + NumPerRank*(Coord+1) + std::min(Remainder, Coord+1);
  }

  return LocalRange;

}

}
