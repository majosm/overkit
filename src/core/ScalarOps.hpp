// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_SCALAR_OPS_HPP_INCLUDED
#define OVK_CORE_SCALAR_OPS_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/ScalarTraits.hpp>

namespace ovk {

template <typename T, OVK_FUNCTION_REQUIRES(core::IsScalar<T>())> constexpr OVK_FORCE_INLINE const T
  &Min(const T &Left, const T &Right) {
  return (Right < Left) ? Right : Left;
}

template <typename T, OVK_FUNCTION_REQUIRES(core::IsScalar<T>())> constexpr OVK_FORCE_INLINE const T
  &Max(const T &Left, const T &Right) {
  return (Left < Right) ? Right : Left;
}

}

#endif
