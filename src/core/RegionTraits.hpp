// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_REGION_TRAITS_HPP_INCLUDED
#define OVK_CORE_REGION_TRAITS_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/Tuple.hpp>

#include <type_traits>

namespace ovk {
namespace core {

template <typename T, typename=void> struct region_traits {
/*
  using coord_type = ???;
  static T MakeEmptyRegion(int NumDims) { return ???; }
  static T UnionRegions(const T &Left, const T &Right) { return ???; }
  static T IntersectRegions(const T &Left, const T &Right) { return ???; }
  static tuple<coord_type> GetRegionLowerCorner(const T &Region) { return ???; }
  static tuple<coord_type> GetRegionUpperCorner(const T &Region) { return ???; }
*/
};

namespace is_region_internal {
template <typename T> using maybe_int = int;
template <typename T> constexpr std::true_type Test(maybe_int<typename
  region_traits<T>::coord_type>) {
  return {};
}
template <typename T> constexpr std::false_type Test(...) { return {}; }
}
template <typename T> constexpr bool IsRegion() {
  return decltype(is_region_internal::Test<T>(0))::value;
}

}}

#endif
