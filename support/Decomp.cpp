// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "support/Decomp.hpp"

#include <ovk/core/ArrayOps.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

namespace support {

void DecomposeDomain(ovk::array_view<const long long> NumPointsPerGrid, int NumProcs,
  ovk::array_view<int,2> GridProcRanges) {

  int NumGrids = NumPointsPerGrid.Count();

  long long TotalPoints = ovk::ArraySum(NumPointsPerGrid);

  long long AvgPointsPerProc = ovk::Max(TotalPoints/NumProcs, 1LL);

  ovk::array<int> NumProcsPerGrid({NumGrids});

  int NumRemainingProcs = NumProcs;
  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    NumProcsPerGrid(iGrid) = ovk::Max(int(NumPointsPerGrid(iGrid)/AvgPointsPerProc), 1);
    NumRemainingProcs -= NumProcsPerGrid(iGrid);
  }

  while (NumRemainingProcs > 0) {
    long long MaxPointsPerProc = 0;
    int iMaxPointsPerProc = 0;
    for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
      long long PointsPerProc = NumPointsPerGrid(iGrid)/NumProcsPerGrid(iGrid);
      if (PointsPerProc > MaxPointsPerProc) {
        MaxPointsPerProc = PointsPerProc;
        iMaxPointsPerProc = iGrid;
      }
    }
    NumProcsPerGrid(iMaxPointsPerProc) += 1;
    --NumRemainingProcs;
  }

  int ProcOffset = 0;
  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    GridProcRanges(iGrid,0) = ProcOffset;
    GridProcRanges(iGrid,1) = ProcOffset + NumProcsPerGrid(iGrid);
    ProcOffset = ovk::Min(ProcOffset + GridProcRanges(iGrid,1)-GridProcRanges(iGrid,0), NumProcs-1);
  }

}

}

#ifdef __cplusplus
extern "C" {
#endif

void support_DecomposeDomain(int NumGrids, const long long *NumPointsPerGrid, int NumProcs, int
  *GridProcRanges) {

  support::DecomposeDomain({NumPointsPerGrid, {NumGrids}}, NumProcs, {GridProcRanges,
    {{NumGrids,2}}});

}

#ifdef __cplusplus
}
#endif

namespace support {

ovk::range CartesianDecomp(int NumDims, const ovk::range &GlobalRange, ovk::comm_view CartComm) {

  OVK_DEBUG_ASSERT(ovk::IsCartComm(CartComm), "Communicator is not Cartesian.");

  ovk::tuple<int> CartDims = ovk::GetCartCommDims(CartComm);
  ovk::tuple<int> CartCoords = ovk::GetCartCommCoords(CartComm);

  ovk::range LocalRange = ovk::MakeEmptyRange(NumDims);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    int NumPerRank = GlobalRange.Size(iDim)/CartDims[iDim];
    int Remainder = GlobalRange.Size(iDim) - CartDims[iDim]*NumPerRank;
    int Begin = GlobalRange.Begin(iDim);
    int Coord = CartCoords[iDim];
    LocalRange.Begin(iDim) = Begin + NumPerRank*Coord + ovk::Min(Remainder, Coord);
    LocalRange.End(iDim) = Begin + NumPerRank*(Coord+1) + ovk::Min(Remainder, Coord+1);
  }

  return LocalRange;

}

ovk::range TriangularDecomp(int NumDims, const ovk::range &GlobalRange, ovk::comm_view Comm) {

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

  ovk::tuple<int> TriangleDims = ovk::MakeUniformTuple<int>(NumDims, 1, 1);
  for (int iDim = 0; iDim < NumDims-1; ++iDim) {
    TriangleDims[iDim] = iPlane+1;
  }
  TriangleDims[NumDims-1] = NumPlanes;

  ovk::tuple<int> PlaneDims = ovk::MakeUniformTuple<int>(NumDims, 1, 1);
  for (int iDim = 0; iDim < NumDims-1; ++iDim) {
    PlaneDims[iDim] = iPlane+1;
  }
  ovk::indexer<int, int, OVK_MAX_DIMS> Indexer(PlaneDims);
  ovk::tuple<int> TriangleCoords = Indexer.ToTuple(iProcInPlane);
  TriangleCoords[NumDims-1] = iPlane;

  ovk::range LocalRange = ovk::MakeEmptyRange(NumDims);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    int NumPerRank = GlobalRange.Size(iDim)/TriangleDims[iDim];
    int Remainder = GlobalRange.Size(iDim) - TriangleDims[iDim]*NumPerRank;
    int Begin = GlobalRange.Begin(iDim);
    int Coord = TriangleCoords[iDim];
    LocalRange.Begin(iDim) = Begin + NumPerRank*Coord + ovk::Min(Remainder, Coord);
    LocalRange.End(iDim) = Begin + NumPerRank*(Coord+1) + ovk::Min(Remainder, Coord+1);
  }

  return LocalRange;

}

}

#ifdef __cplusplus
extern "C" {
#endif

void support_CartesianDecomp(int NumDims, const int *GlobalBegin, const int *GlobalEnd, MPI_Comm
  CartComm, int *LocalBegin, int *LocalEnd) {

  ovk::range LocalRange = support::CartesianDecomp(NumDims, {GlobalBegin, GlobalEnd}, CartComm);

  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    LocalBegin[iDim] = LocalRange.Begin(iDim);
    LocalEnd[iDim] = LocalRange.End(iDim);
  }

}

void support_TriangularDecomp(int NumDims, const int *GlobalBegin, const int *GlobalEnd, MPI_Comm
  Comm, int *LocalBegin, int *LocalEnd) {

  ovk::range LocalRange = support::TriangularDecomp(NumDims, {GlobalBegin, GlobalEnd}, Comm);

  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    LocalBegin[iDim] = LocalRange.Begin(iDim);
    LocalEnd[iDim] = LocalRange.End(iDim);
  }

}

#ifdef __cplusplus
}
#endif
