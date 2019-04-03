// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_REQUIRES_HPP_INCLUDED
#define OVK_CORE_REQUIRES_HPP_INCLUDED

#include <ovk/core/Global.hpp>

#include <type_traits>

namespace ovk {
namespace core {

namespace requires_internal {

enum class function_requires_type {
  value
};
constexpr const function_requires_type function_requires_value = function_requires_type::value;

template <bool Condition> struct function_requires {};
template <> struct function_requires<true> {
  using type = function_requires_type;
};

}

#define OVK_FUNCTION_REQUIRES(...) typename ::ovk::core::requires_internal::function_requires<\
  (__VA_ARGS__)>::type = ::ovk::core::requires_internal::function_requires_value
#define OVK_FUNCDECL_REQUIRES(...) typename ::ovk::core::requires_internal::function_requires<\
  (__VA_ARGS__)>::type = ::ovk::core::requires_internal::function_requires_value
#define OVK_FUNCDEF_REQUIRES(...) typename ::ovk::core::requires_internal::function_requires<\
  (__VA_ARGS__)>::type

namespace requires_internal {

enum class alias_requires_type {
  value
};
constexpr const alias_requires_type alias_requires_value = alias_requires_type::value;

template <bool Condition> struct alias_requires {};
template <> struct alias_requires<true> {
  using type = alias_requires_type;
};

}

#define OVK_ALIAS_REQUIRES(...) typename ::ovk::core::requires_internal::alias_requires<\
  (__VA_ARGS__)>::type = ::ovk::core::requires_internal::alias_requires_value

namespace requires_internal {

template <bool Condition> struct specialization_requires {};
template <> struct specialization_requires<true> {
  using type = void;
};

}

#define OVK_SPECIALIZATION_REQUIRES(...) typename ::ovk::core::requires_internal::\
  specialization_requires<(__VA_ARGS__)>::type

}}

#endif
