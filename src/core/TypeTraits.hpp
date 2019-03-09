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

// Apply const & reference qualifiers of one type to another type
namespace mimicked_ref_internal {
template <typename T, typename U> struct helper {
  using type = U;
};
template <typename T, typename U> struct helper<T &, U> {
  using type = typename std::add_lvalue_reference<U>::type;
};
template <typename T, typename U> struct helper<const T &, U> {
  using type = typename std::add_lvalue_reference<typename std::add_const<U>::type>::type;
};
template <typename T, typename U> struct helper<T &&, U> {
  using type = typename std::add_rvalue_reference<U>::type;
};
template <typename T, typename U> struct helper<const T &&, U> {
  using type = typename std::add_rvalue_reference<typename std::add_const<U>::type>::type;
};
}
template <typename T, typename U> using mimicked_ref = typename mimicked_ref_internal::helper<T,
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
