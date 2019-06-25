// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISPERSE_MAP_HPP_INCLUDED
#define OVK_CORE_DISPERSE_MAP_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>

#include <mpi.h>

namespace ovk {
namespace core {

class disperse_map {

public:

  disperse_map() = default;
  disperse_map(array<int,2> Points);

  floating_ref<const disperse_map> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<disperse_map> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  long long Count() const { return Points_.Size(2); }

  const array<int,2> &Points() const { return Points_; }

private:

  floating_ref_generator FloatingRefGenerator_;

  array<int,2> Points_;

};

}}

#endif
