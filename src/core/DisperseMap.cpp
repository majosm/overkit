// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/DisperseMap.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <utility>

namespace ovk {
namespace core {

disperse_map::disperse_map(array<int,2> Points):
  Points_(std::move(Points))
{}

}}
