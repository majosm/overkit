// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ELEM_SET_HPP_INCLUDED
#define OVK_CORE_ELEM_SET_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Set.hpp>

namespace ovk {

template <typename KeyElementType, int KeyRank, array_layout Layout=array_layout::ROW_MAJOR> using
  elem_set = set<elem<KeyElementType,KeyRank>, elem_less<KeyElementType, KeyRank, Layout>>;

template <typename KeyElementType, int KeyRank> using elem_set_r = elem_set<KeyElementType, KeyRank,
  array_layout::ROW_MAJOR>;
template <typename KeyElementType, int KeyRank> using elem_set_c = elem_set<KeyElementType, KeyRank,
  array_layout::COLUMN_MAJOR>;

}

#endif
