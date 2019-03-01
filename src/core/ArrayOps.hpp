// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ARRAY_OPS_HPP_INCLUDED
#define OVK_CORE_ARRAY_OPS_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/Requires.hpp>

#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(!std::is_const<T>::value
  )> void ArrayFill(array_view<T, Rank, Layout> View, const T &Value) {
  View.Fill(Value);
}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>())> void ArrayFill(ArrayType
  &Array, const array_value_type<ArrayType> &Value) {
  MakeArrayView(Array).Fill(Value);
}

template <typename T, int Rank, array_layout Layout, typename IterType, OVK_FUNCTION_REQUIRES(
  !std::is_const<T>::value && IsInputIterator<IterType>() && std::is_convertible<
  iterator_deref_type<IterType>, T>::value)> void ArrayFill(array_view<T, Rank, Layout> View,
  IterType First) {
  View.Fill(First);
}

template <typename ArrayType, typename IterType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>() &&
  IsInputIterator<IterType>() && std::is_convertible<iterator_deref_type<IterType>,
  array_value_type<ArrayType>>::value)> void ArrayFill(ArrayType &Array, IterType First) {
  MakeArrayView(Array).Fill(First);
}

template <typename T, typename U, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(
  !std::is_const<T>::value && std::is_convertible<typename std::remove_const<U>::type, T>::value)>
  void ArrayFill(array_view<T, Rank, Layout> View, array_view<U, Rank, Layout> SourceView) {
  View.Fill(SourceView);
}

template <typename ArrayType, typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(
  IsArray<ArrayType>() && ArrayHasFootprint<ArrayType, Rank, Layout>() && std::is_convertible<
  typename std::remove_const<T>::type, array_value_type<ArrayType>>::value)> void ArrayFill(
  ArrayType &Array, array_view<T, Rank, Layout> SourceView) {
  MakeArrayView(Array).Fill(SourceView);
}

template <typename T, int Rank, array_layout Layout, typename SourceArrayRefType,
  OVK_FUNCTION_REQUIRES(!std::is_const<T>::value && IsArray<typename std::decay<SourceArrayRefType
  >::type>() && !IsIterator<typename std::decay<SourceArrayRefType>::type>() && ArrayHasFootprint<
  typename std::decay<SourceArrayRefType>::type, Rank, Layout>() && std::is_convertible<
  array_access_type<SourceArrayRefType &&>, T>::value)> void ArrayFill(array_view<T, Rank, Layout>
  View, SourceArrayRefType &&SourceArray) {
  View.Fill(std::forward<SourceArrayRefType>(SourceArray));
}

template <typename ArrayType, typename SourceArrayRefType, OVK_FUNCTION_REQUIRES(IsArray<
  ArrayType>() && IsArray<typename std::decay<SourceArrayRefType>::value>() && !IsIterator<typename
  std::decay<SourceArrayRefType>::type>() && ArraysAreSimilar<ArrayType, typename std::decay<
  SourceArrayRefType>::type>() && std::is_convertible<array_access_type<SourceArrayRefType &&>,
  array_value_type<ArrayType>>::value)> void ArrayFill(ArrayType &Array, SourceArrayRefType
  &&SourceArray) {
  MakeArrayView(Array).Fill(std::forward<SourceArrayRefType>(SourceArray));
}

}}

#endif
