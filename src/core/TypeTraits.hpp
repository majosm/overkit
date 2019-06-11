// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_TYPE_TRAITS_HPP_INCLUDED
#define OVK_CORE_TYPE_TRAITS_HPP_INCLUDED

#include <ovk/core/Global.hpp>

#include <cstdint>
#include <type_traits>

namespace ovk {
namespace core {

// Remove cv and reference qualifiers
template <typename T> using remove_cvref = typename std::remove_cv<typename std::remove_reference<T
  >::type>::type;

// Apply cv and reference qualifiers of one type to another type
namespace mimic_cvref_internal {
template <typename T, typename U> struct helper {
  using U1 = typename std::conditional<std::is_const<typename std::remove_reference<T>::type>
    ::value, typename std::add_const<U>::type, U>::type;
  using U2 = typename std::conditional<std::is_volatile<typename std::remove_reference<T>::type>
    ::value, typename std::add_volatile<U1>::type, U1>::type;
  using U3 = typename std::conditional<std::is_lvalue_reference<T>::value, typename
    std::add_lvalue_reference<U2>::type, U2>::type;
  using U4 = typename std::conditional<std::is_rvalue_reference<T>::value, typename
    std::add_rvalue_reference<U3>::type, U3>::type;
  using type = U4;
};
}
template <typename T, typename U> using mimic_cvref = typename mimic_cvref_internal::helper<T,
U>::type;

// Check if a function can be called with a specified set of argument types
namespace is_callable_with_internal {
template <typename FRef, typename... Args> constexpr std::true_type Test(decltype(
  std::declval<FRef>()(std::declval<Args>()...)) *) { return {}; }
template <typename FRef, typename... Args> constexpr std::false_type Test(...) { return {}; }
}
template <typename FRef, typename... Args> constexpr bool IsCallableWith() {
  return decltype(is_callable_with_internal::Test<FRef, Args...>(nullptr))::value;
}

// Unique integer value for every type
using type_id_type = std::uintptr_t;
namespace type_id_internal {
template <typename T> struct helper {
  static int dummy;
};
}
template <typename T> int type_id_internal::helper<T>::dummy;
template <typename T> constexpr type_id_type GetTypeID() {
  return reinterpret_cast<std::uintptr_t>(&type_id_internal::helper<T>::dummy);
}

}}

#endif
