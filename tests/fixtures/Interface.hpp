// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_TESTS_FIXTURES_INTERFACE_HPP_LOADED
#define OVK_TESTS_FIXTURES_INTERFACE_HPP_LOADED

#include "support/Decomp.hpp"

#include <overkit.hpp>

namespace tests {

ovk::domain Interface2D(ovk::comm_view Comm, const ovk::box &Bounds, const ovk::tuple<int> &Size,
  const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage);

ovk::domain Interface2DManualConnectivity(ovk::comm_view Comm, const ovk::box &Bounds, const
  ovk::tuple<int> &Size, const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage);

ovk::domain Interface3D(ovk::comm_view Comm, const ovk::box &Bounds, const ovk::tuple<int> &Size,
  const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage);

ovk::domain Interface3DManualConnectivity(ovk::comm_view Comm, const ovk::box &Bounds, const
  ovk::tuple<int> &Size, const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage);

}

#endif
