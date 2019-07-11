// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_FIELD_HPP_INCLUDED
#define OVK_CORE_FIELD_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Global.hpp>

namespace ovk {

template <typename T> using field = array<T,MAX_DIMS,array_layout::COLUMN_MAJOR>;
template <typename T> using field_view = array_view<T,MAX_DIMS,array_layout::COLUMN_MAJOR>;

namespace core {
template <typename T> constexpr bool IsField() {
  return IsArray<T>() && ArrayRank<T>() == MAX_DIMS && ArrayLayout<T>() ==
    array_layout::COLUMN_MAJOR;
}
}

}

#endif
