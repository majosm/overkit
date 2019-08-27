// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ELEM_MAP_HPP_INCLUDED
#define OVK_CORE_ELEM_MAP_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Map.hpp>

namespace ovk {

template <typename KeyElementType, int KeyRank, typename ValueType, array_layout
  Layout=array_layout::ROW_MAJOR, bool Contiguous=MapContiguousDefault<ValueType>()> using
  elem_map = map<elem<KeyElementType,KeyRank>, ValueType, elem_less<KeyElementType, KeyRank,
  Layout>, Contiguous>;

template <typename KeyElementType, int KeyRank, typename ValueType, array_layout
  Layout=array_layout::ROW_MAJOR> using elem_map_contig = elem_map<KeyElementType, KeyRank,
  ValueType, Layout, true>;
template <typename KeyElementType, int KeyRank, typename ValueType, array_layout
  Layout=array_layout::ROW_MAJOR> using elem_map_noncontig = elem_map<KeyElementType, KeyRank,
  ValueType, Layout, false>;

template <typename KeyElementType, int KeyRank, typename ValueType, bool Contiguous=
  MapContiguousDefault<ValueType>()> using elem_map_r = elem_map<KeyElementType, KeyRank,
  ValueType, array_layout::ROW_MAJOR, Contiguous>;

template <typename KeyElementType, int KeyRank, typename ValueType> using elem_map_r_contig =
  elem_map_r<KeyElementType, KeyRank, ValueType, true>;
template <typename KeyElementType, int KeyRank, typename ValueType> using elem_map_r_noncontig =
  elem_map_r<KeyElementType, KeyRank, ValueType, false>;

template <typename KeyElementType, int KeyRank, typename ValueType, bool Contiguous=
  MapContiguousDefault<ValueType>()> using elem_map_c = elem_map<KeyElementType, KeyRank,
  ValueType, array_layout::COLUMN_MAJOR, Contiguous>;

template <typename KeyElementType, int KeyRank, typename ValueType> using elem_map_c_contig =
  elem_map_c<KeyElementType, KeyRank, ValueType, true>;
template <typename KeyElementType, int KeyRank, typename ValueType> using elem_map_c_noncontig =
  elem_map_c<KeyElementType, KeyRank, ValueType, false>;

}

#endif
