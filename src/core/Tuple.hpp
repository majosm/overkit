// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_TUPLE_HPP_INCLUDED
#define OVK_CORE_TUPLE_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>

namespace ovk {

template <typename T> using tuple = elem<T,MAX_DIMS>;

template <typename T> constexpr OVK_FORCE_INLINE tuple<T> MakeUniformTuple(const T &Value) {
  return MakeUniformElem<T,MAX_DIMS>(Value);
}

}

#endif
