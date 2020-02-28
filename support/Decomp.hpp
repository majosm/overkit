// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_DECOMP_HPP_LOADED
#define OVK_SUPPORT_DECOMP_HPP_LOADED

#include <support/Decomp.h>

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

namespace support {

using ovk::core::CreateCartesianDecompDims;

void DecomposeDomain(ovk::array_view<const long long> NumPointsPerGrid, int NumProcs,
  ovk::array_view<int,2> GridProcRanges);

ovk::range CartesianDecomp(int NumDims, const ovk::range &GlobalRange, ovk::comm_view CartComm);
ovk::range TriangularDecomp(int NumDims, const ovk::range &GlobalRange, ovk::comm_view Comm);

}

#endif
