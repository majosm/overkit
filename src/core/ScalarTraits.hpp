// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_SCALAR_TRAITS_HPP_INCLUDED
#define OVK_CORE_SCALAR_TRAITS_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>

#include <type_traits>

namespace ovk {

template <typename T, typename=void> struct scalar_traits {
  // Nothing here of substance at the moment; using type as a placeholder for IsScalar to detect
  // using type = T;
};

template <typename T> struct scalar_traits<T, OVK_SPECIALIZATION_REQUIRES(
  std::is_arithmetic<T>::value)> {
  using type = T;
};

namespace core {

namespace scalar_traits_internal {
template <typename T> constexpr std::true_type IsScalarTest(typename scalar_traits<T>::type *) {
  return {};
}
template <typename T> constexpr std::false_type IsScalarTest(...) { return {}; }
}
template <typename T> constexpr bool IsScalar() {
  return decltype(scalar_traits_internal::IsScalarTest<T>(nullptr))::value;
}

}

}

#endif
