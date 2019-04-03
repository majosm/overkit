// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_TYPE_SEQUENCE_HPP_INCLUDED
#define OVK_CORE_TYPE_SEQUENCE_HPP_INCLUDED

#include <cstddef>

namespace ovk {
namespace core {

template <typename... Sequence>  struct type_sequence {
  static constexpr std::size_t size() { return sizeof...(Sequence); }
};

namespace type_sequence_internal {

template <typename Left, typename Right> struct repeated_type_sequence_merge;
template <typename... LeftSequence, typename... RightSequence>
struct repeated_type_sequence_merge<type_sequence<LeftSequence...>,
  type_sequence<RightSequence...>> {
  using type = type_sequence<LeftSequence..., RightSequence...>;
};

template <typename T, std::size_t N> struct repeated_type_sequence_generator {
  using type = typename repeated_type_sequence_merge<
    typename repeated_type_sequence_generator<T,N/2>::type,
    typename repeated_type_sequence_generator<T,N-N/2>::type
  >::type;
};

template <typename T> struct repeated_type_sequence_generator<T,0> {
  using type = type_sequence<>;
};

template <typename T> struct repeated_type_sequence_generator<T,1> {
  using type = type_sequence<T>;
};

}

template <typename T, std::size_t N> using repeated_type_sequence_of_size
  = typename type_sequence_internal::repeated_type_sequence_generator<T,N>::type;

}}

#endif
