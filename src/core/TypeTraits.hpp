// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_TYPE_TRAITS_HPP_INCLUDED
#define OVK_CORE_TYPE_TRAITS_HPP_INCLUDED

#include <ovk/core/Global.hpp>

#include <type_traits>

namespace ovk {
namespace core {

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

}}

#endif
