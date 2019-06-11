// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_DECOMP_HPP_LOADED
#define OVK_SUPPORT_DECOMP_HPP_LOADED

#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

namespace support {

ovk::range CartesianDecomp(const ovk::cart &Cart, ovk::comm_view Comm, const ovk::tuple<int>
  &CartDims);
ovk::range TriangularDecomp(const ovk::cart &Cart, ovk::comm_view Comm);

}

#endif
