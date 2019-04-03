// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_INTEGER_SEQUENCE_HPP_INCLUDED
#define OVK_CORE_INTEGER_SEQUENCE_HPP_INCLUDED

#include <ovk/core/Global.hpp>

#include <cstddef>

namespace ovk {
namespace core {

template <typename T, T... Sequence> struct integer_sequence {
  using value_type = T;
  static constexpr std::size_t size() { return sizeof...(Sequence); }
};

template <std::size_t... Sequence> using index_sequence = integer_sequence<std::size_t,
  Sequence...>;

namespace integer_sequence_internal {

template <typename LeftSequence, typename RightSequence> struct index_sequence_merge;
template <std::size_t... LeftIndices, std::size_t... RightIndices> struct index_sequence_merge<
  index_sequence<LeftIndices...>, index_sequence<RightIndices...>
> {
  using type = index_sequence<LeftIndices..., RightIndices...>;
};

template <typename Sequence, std::size_t N> struct index_sequence_shift;
template <std::size_t... Indices, std::size_t N> struct index_sequence_shift<index_sequence<
  Indices...>, N> {
  using type = index_sequence<(N + Indices)...>;
};

template <std::size_t N> struct index_sequence_generator {
  using type = typename index_sequence_merge<
    typename index_sequence_generator<N/2>::type,
    typename index_sequence_shift<
      typename index_sequence_generator<N-N/2>::type, N/2>::type
  >::type;
};

template <> struct index_sequence_generator<0> {
  using type = index_sequence<>;
};

template <> struct index_sequence_generator<1> {
  using type = index_sequence<0>;
};

}

template <std::size_t N> using index_sequence_of_size
  = typename integer_sequence_internal::index_sequence_generator<N>::type;

}}

#endif
