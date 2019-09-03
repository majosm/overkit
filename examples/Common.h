// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXAMPLES_COMMON_H_LOADED
#define OVK_EXAMPLES_COMMON_H_LOADED

#include <support/Constants.h>
#include <support/Decomp.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EXAMPLES_PI SUPPORT_PI

void examples_DecomposeDomain(int NumGrids, const long long *NumPointsPerGrid, int NumProcs, int
  *GridProcRanges) {

  support_DecomposeDomain(NumGrids, NumPointsPerGrid, NumProcs, GridProcRanges);

}

void examples_CreateCartesianDecompDims(int Size, int NumDims, int *Dims) {

  support_CreateCartesianDecompDims(Size, NumDims, Dims);

}

void examples_CartesianDecomp(int NumDims, const int *Size, MPI_Comm CartComm, int *LocalRange) {

  int Zero[] = {0,0,0};

  support_CartesianDecomp(NumDims, Zero, Size, CartComm, LocalRange, LocalRange+3);

}

#ifdef __cplusplus
}
#endif

#endif
