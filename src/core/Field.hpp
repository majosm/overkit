// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_FIELD_HPP_INCLUDED
#define OVK_CORE_FIELD_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Global.hpp>

namespace ovk {

template <typename T> using field = array_c<T,MAX_DIMS>;
template <typename T> using field_view = array_view_c<T,MAX_DIMS>;

using field_indexer = indexer_c<long long,int,MAX_DIMS>;

namespace core {
template <typename T> constexpr bool IsField() {
  return IsArray<T>() && ArrayRank<T>() == MAX_DIMS && ArrayLayout<T>() ==
    array_layout::COLUMN_MAJOR;
}
}

}

#endif
