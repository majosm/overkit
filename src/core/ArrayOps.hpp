// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ARRAY_OPS_HPP_INCLUDED
#define OVK_CORE_ARRAY_OPS_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/Requires.hpp>

#include <type_traits>
#include <utility>

namespace ovk {

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(!std::is_const<T>::value
  )> void ArrayFill(array_view<T, Rank, Layout> View, const T &Value) {

  View.Fill(Value);

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>())> void ArrayFill(
  ArrayType &Array, const core::array_value_type<ArrayType> &Value) {

  long long NumValues = core::ArrayCount(Array);

  for (long long i = 0; i < NumValues; ++i) {
    Array[i] = Value;
  }

}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(!std::is_const<T>::value
  )> void ArrayFill(array_view<T, Rank, Layout> View, std::initializer_list<T> ValuesList) {

  View.Fill(ValuesList);

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>())> void ArrayFill(
  ArrayType &Array, std::initializer_list<core::array_value_type<ArrayType>> ValuesList) {

  long long NumValues = core::ArrayCount(Array);

  auto Iter = ValuesList.begin();
  for (long long i = 0; i < NumValues; ++i) {
    Array[i] = *Iter++;
  }

}

template <typename T, int Rank, array_layout Layout, typename IterType, OVK_FUNCTION_REQUIRES(
  !std::is_const<T>::value && core::IsInputIterator<IterType>() && std::is_convertible<
  core::iterator_reference_type<IterType>, T>::value)> void ArrayFill(array_view<T, Rank, Layout>
  View, IterType First) {

  View.Fill(First);

}

template <typename ArrayType, typename IterType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  core::IsInputIterator<IterType>() && std::is_convertible<core::iterator_reference_type<IterType>,
  core::array_value_type<ArrayType>>::value)> void ArrayFill(ArrayType &Array, IterType First) {

  long long NumValues = core::ArrayCount(Array);

  IterType Iter = First;
  for (long long i = 0; i < NumValues; ++i) {
    Array[i] = *Iter++;
  }

}

template <typename T, typename U, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(
  !std::is_const<T>::value && std::is_convertible<typename std::remove_const<U>::type, T>::value)>
  void ArrayFill(array_view<T, Rank, Layout> View, array_view<U, Rank, Layout> SourceView) {

  View.Fill(SourceView);

}

template <typename ArrayType, typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(
  core::IsArray<ArrayType>() && core::ArrayHasFootprint<ArrayType, Rank, Layout>() &&
  std::is_convertible<typename std::remove_const<T>::type, core::array_value_type<ArrayType>>::
  value)> void ArrayFill(ArrayType &Array, array_view<T, Rank, Layout> SourceView) {

  long long NumValues = core::ArrayCount(Array);

  for (long long i = 0; i < NumValues; ++i) {
    Array[i] = SourceView[i];
  }

}

template <typename T, int Rank, array_layout Layout, typename SourceArrayRefType,
  OVK_FUNCTION_REQUIRES(!std::is_const<T>::value && core::IsArray<core::remove_cvref<
  SourceArrayRefType>>() && !core::IsIterator<typename std::decay<SourceArrayRefType>::type>()
  && core::ArrayHasFootprint<core::remove_cvref<SourceArrayRefType>, Rank, Layout>() &&
  std::is_convertible<core::array_access_type<SourceArrayRefType &&>, T>::value)> void ArrayFill(
  array_view<T, Rank, Layout> View, SourceArrayRefType &&SourceArray) {

  View.Fill(std::forward<SourceArrayRefType>(SourceArray));

}

template <typename ArrayType, typename SourceArrayRefType, OVK_FUNCTION_REQUIRES(core::IsArray<
  ArrayType>() && core::IsArray<core::remove_cvref<SourceArrayRefType>>() && !core::IsIterator<
  typename std::decay<SourceArrayRefType>::type>() && core::ArraysAreSimilar<ArrayType,
  core::remove_cvref<SourceArrayRefType>>() && std::is_convertible<core::array_access_type<
  SourceArrayRefType &&>, core::array_value_type<ArrayType>>::value)> void ArrayFill(ArrayType
  &Array, SourceArrayRefType &&SourceArray) {

  long long NumValues = core::ArrayCount(Array);

  for (long long i = 0; i < NumValues; ++i) {
    Array[i] = static_cast<core::array_access_type<SourceArrayRefType &&>>(SourceArray[i]);
  }

}

}

#endif
