// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_FOR_EACH_HPP_INCLUDED
#define OVK_CORE_FOR_EACH_HPP_INCLUDED

#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/TypeSequence.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

namespace for_each_internal {

namespace is_callable_with_tuple_elements_internal {
template <typename FRef, typename... TupleElementTypes> constexpr bool Helper(type_sequence<
  TupleElementTypes...>) {
  return IsCallableWith<FRef, TupleElementTypes...>();
}
}
template <typename FRef, typename T, int N> constexpr bool IsCallableWithTupleElements() {
  return is_callable_with_tuple_elements_internal::Helper<FRef>(
    repeated_type_sequence_of_size<T, N>());
}

template <int N, int I, array_layout Layout> struct helper;
template <int N, int I> struct helper<N, I, array_layout::ROW_MAJOR> {
  template <typename T, typename FRef, typename... TupleElementTypes> static OVK_FORCE_INLINE void
    ForEach(const interval<T,N> &Interval, FRef &&Func, TupleElementTypes... TupleElements) {
    for (T i = Interval.Begin(N-1-I); i < Interval.End(N-1-I); ++i) {
      helper<N, I-1, array_layout::ROW_MAJOR>::ForEach(Interval, std::forward<FRef>(Func),
        TupleElements..., i);
    }
  }
};
template <int N, int I> struct helper<N, I, array_layout::COLUMN_MAJOR> {
  template <typename T, typename FRef, typename... TupleElementTypes> static OVK_FORCE_INLINE void
    ForEach(const interval<T,N> &Interval, FRef &&Func, TupleElementTypes... TupleElements) {
    for (T i = Interval.Begin(I); i < Interval.End(I); ++i) {
      helper<N, I-1, array_layout::COLUMN_MAJOR>::ForEach(Interval, std::forward<FRef>(Func),
        i, TupleElements...);
    }
  }
};
template <int N> struct helper<N, 0, array_layout::ROW_MAJOR> {
  template <typename T, typename FRef, typename... TupleElementTypes, OVK_FUNCTION_REQUIRES(
    IsCallableWithTupleElements<FRef &&, T, N>())> static OVK_FORCE_INLINE void ForEach(const
    interval<T,N> &Interval, FRef &&Func, TupleElementTypes... TupleElements) {
    for (T i = Interval.Begin(N-1); i < Interval.End(N-1); ++i) {
      std::forward<FRef>(Func)(TupleElements...,i);
    }
  }
  template <typename T, typename FRef, typename... TupleElementTypes, OVK_FUNCTION_REQUIRES(
    !IsCallableWithTupleElements<FRef &&, T, N>() && IsCallableWith<FRef &&, const elem<T,N> &>())>
    static OVK_FORCE_INLINE void ForEach(const interval<T,N> &Interval, FRef &&Func,
    TupleElementTypes... TupleElements) {
    for (T i = Interval.Begin(N-1); i < Interval.End(N-1); ++i) {
      std::forward<FRef>(Func)(elem<T,N>(TupleElements...,i));
    }
  }
};
template <int N> struct helper<N, 0, array_layout::COLUMN_MAJOR> {
  template <typename T, typename FRef, typename... TupleElementTypes, OVK_FUNCTION_REQUIRES(
    IsCallableWithTupleElements<FRef &&, T, N>())> static OVK_FORCE_INLINE void ForEach(const
    interval<T,N> &Interval, FRef &&Func, TupleElementTypes... TupleElements) {
    for (T i = Interval.Begin(0); i < Interval.End(0); ++i) {
      std::forward<FRef>(Func)(i,TupleElements...);
    }
  }
  template <typename T, typename FRef, typename... TupleElementTypes, OVK_FUNCTION_REQUIRES(
    !IsCallableWithTupleElements<FRef &&, T, N>() && IsCallableWith<FRef &&, const elem<T,N> &>())>
    static OVK_FORCE_INLINE void ForEach(const interval<T,N> &Interval, FRef &&Func,
    TupleElementTypes... TupleElements) {
    for (T i = Interval.Begin(0); i < Interval.End(0); ++i) {
      std::forward<FRef>(Func)(elem<T,N>(i,TupleElements...));
    }
  }
};

}

template <array_layout Layout=array_layout::ROW_MAJOR, typename T, int N, typename FRef,
  OVK_FUNCTION_REQUIRES(IsCallableWith<FRef &&, const elem<T,N> &>() || for_each_internal::
  IsCallableWithTupleElements<FRef &&, T, N>())> OVK_FORCE_INLINE void ForEach(const interval<T,N>
  &Interval, FRef &&Func) {
  // Have to use N-1 and count down because we can specialize on (N, 0) but not (N, N-1)
  return for_each_internal::helper<N, N-1, Layout>::ForEach(Interval, std::forward<FRef>(Func));
}

}}

#endif
