// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
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

// Check if a type can be called with a specified set of argument types
namespace is_callable_with_internal {
template <typename FRef, typename... Args> constexpr std::true_type Test(decltype(
  std::declval<FRef>()(std::declval<Args>()...)) *) { return {}; }
template <typename FRef, typename... Args> constexpr std::false_type Test(...) { return {}; }
}
template <typename FRef, typename... Args> constexpr bool IsCallableWith() {
  return decltype(is_callable_with_internal::Test<FRef, Args...>(nullptr))::value;
}

// Check if a type can be called as if it was a function with a given interface
namespace is_callable_as_internal {
template <typename Result> constexpr int ImplicitConvertToResult(const Result &) { return 0; }
template <typename Signature> struct helper;
template <typename Result, typename... Args> struct helper<Result(Args...)> {
  template <typename FRef> static constexpr std::true_type Test(decltype(ImplicitConvertToResult<
    Result>(std::declval<FRef>()(std::declval<Args>()...))) *) { return {}; }
  template <typename FRef> static constexpr std::false_type Test(...) { return {}; }
};
}
template <typename FRef, typename Signature> constexpr bool IsCallableAs() {
  return decltype(is_callable_as_internal::helper<Signature>::template Test<FRef>(nullptr))::value;
}

// Check if a list of arguments to be forwarded to a type's constructor is an lvalue reference to
// an object of the same type
namespace is_copy_argument_internal {
template <typename T, typename... Args> struct helper : std::false_type {};
template <typename T, typename Arg> struct helper<T, Arg> : std::integral_constant<bool,
  std::is_same<remove_cvref<Arg>, T>::value && std::is_lvalue_reference<Arg>::value> {};
}
template <typename T, typename... Args> constexpr bool IsCopyArgument() {
  return is_copy_argument_internal::helper<T, Args...>::value;
}

// Check if a list of arguments to be forwarded to a type's constructor is an rvalue reference to
// an object of the same type
namespace is_move_argument_internal {
template <typename T, typename... Args> struct helper : std::false_type {};
template <typename T, typename Arg> struct helper<T, Arg> : std::integral_constant<bool,
  std::is_same<remove_cvref<Arg>, T>::value && std::is_rvalue_reference<Arg>::value> {};
}
template <typename T, typename... Args> constexpr bool IsMoveArgument() {
  return is_move_argument_internal::helper<T, Args...>::value;
}

// Check if a list of arguments to be forwarded to a type's constructor is a reference to an object
// of the same type
template <typename T, typename... Args> constexpr bool IsCopyOrMoveArgument() {
  return IsCopyArgument<T, Args...>() || IsMoveArgument<T, Args...>();
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
