// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DECOMP_HPP_INCLUDED
#define OVK_CORE_DECOMP_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/DistributedRegionHash.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Map.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

namespace ovk {
namespace core {

using decomp_hash = distributed_region_hash<int>;
using decomp_hash_bin = distributed_region_hash_bin<int>;
using decomp_hash_region_data = distributed_region_data<int>;

decomp_hash CreateDecompHash(int NumDims, comm_view Comm, const range &LocalRange);

array<int> DetectNeighbors(const cart &Cart, comm_view Comm, const range &LocalRange, const
  decomp_hash &Hash);

range ExtendLocalRange(const cart &Cart, const range &LocalRange, int ExtendAmount);

cart CartPointToCell(const cart &Cart);
range RangePointToCell(const cart &Cart, const range &LocalRange);

struct decomp_info {
  range LocalRange;
  range ExtendedRange;
};

map<int,decomp_info> RetrieveDecompInfo(comm_view Comm, array_view<const int> Ranks, const
  range &LocalRange, const range &ExtendedRange);

}}

#endif
