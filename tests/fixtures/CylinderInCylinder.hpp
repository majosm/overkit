// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_TESTS_FIXTURES_CYLINDER_IN_CYLINDER_HPP_LOADED
#define OVK_TESTS_FIXTURES_CYLINDER_IN_CYLINDER_HPP_LOADED

#include <overkit.hpp>

namespace tests {

ovk::domain CylinderInCylinder(ovk::comm_view Comm, int Size, int OverlapAmount, bool Stagger,
  const ovk::tuple<bool> &DecompDirs, ovk::periodic_storage PeriodicStorage);

}

#endif
