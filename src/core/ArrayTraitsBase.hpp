// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ARRAY_TRAITS_BASE_HPP_INCLUDED
#define OVK_CORE_ARRAY_TRAITS_BASE_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/IntegerSequence.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <type_traits>

namespace ovk {

// This stuff needs to be defined in here to resolve circular dependency (array traits <-> elem)

template <typename T, typename=void> struct array_traits {
/*
  using value_type = ???;

  static constexpr int Rank = ???;

  static constexpr array_layout Layout = ???;

  // For statically-sized arrays
  template <int iDim> static constexpr long long ExtentBegin() { return ???; }
  template <int iDim> static constexpr long long ExtentEnd() { return ???; }

  // For runtime-sized arrays
  template <int iDim> static long long ExtentBegin(const T &Array) { return ???; }
  template <int iDim> static long long ExtentEnd(const T &Array) { return ???; }

  static const T *Data(const T &Array) { return ???; }
  static T *Data(T &Array) { return ???; }
*/
};

namespace core {

namespace is_array_internal {
template <typename T> using maybe_int = int;
template <typename T> constexpr std::true_type Test(maybe_int<typename array_traits<T>::value_type>)
  {
  return {};
}
template <typename T> constexpr std::false_type Test(...) { return {}; }
}
template <typename T> constexpr bool IsArray() {
  return decltype(is_array_internal::Test<T>(0))::value;
}

template <typename T, OVK_FUNCTION_REQUIRES(IsArray<T>())> constexpr int ArrayRank() {
  return array_traits<T>::Rank;
}
template <typename T, OVK_FUNCTION_REQUIRES(!IsArray<T>())> constexpr int ArrayRank() {
  return 1;
}

namespace array_access_type_internal {
template <typename TRef, typename=void> struct helper;
template <typename TRef> struct helper<TRef, OVK_SPECIALIZATION_REQUIRES(IsArray<remove_cvref<TRef>
  >())> {
  using type = mimic_cvref<TRef, typename array_traits<remove_cvref<TRef>>::value_type>;
};
template <typename TRef> struct helper<TRef, OVK_SPECIALIZATION_REQUIRES(!IsArray<
  remove_cvref<TRef>>())> {
  using type = std::false_type;
};
}
template <typename TRef> using array_access_type = typename array_access_type_internal::helper<TRef
  >::type;

namespace array_has_static_extents_internal {
template <typename T> using maybe_int = int;
template <typename T> constexpr std::true_type Test(maybe_int<decltype(array_traits<T>::template
  ExtentBegin<0>())>) {
  return {};
}
template <typename T> constexpr std::false_type Test(...) { return {}; }
}
template <typename T> constexpr bool ArrayHasStaticExtents() {
  return decltype(array_has_static_extents_internal::Test<T>(0))::value;
}

namespace array_has_runtime_extents_internal {
template <typename T> using maybe_int = int;
template <typename T> constexpr std::true_type Test(maybe_int<decltype(array_traits<T>::template
  ExtentBegin<0>(std::declval<T>()))>) {
  return {};
}
template <typename T> constexpr std::false_type Test(...) { return {}; }
}
template <typename T> constexpr bool ArrayHasRuntimeExtents() {
  return decltype(array_has_runtime_extents_internal::Test<T>(0))::value;
}

namespace static_array_has_extents_begin_internal {
template <typename T, long long BeginElement, std::size_t Dim, OVK_FUNCTION_REQUIRES(IsArray<T>()
  && ArrayHasStaticExtents<T>())> constexpr bool Helper(integer_sequence<long long, BeginElement>,
  index_sequence<Dim>) {
  return array_traits<T>::template ExtentBegin<Dim>() == BeginElement;
}
template <typename T, long long BeginElement1, long long BeginElement2, long long...
  RemainingElements, std::size_t Dim1, std::size_t Dim2, std::size_t... RemainingDims,
  OVK_FUNCTION_REQUIRES(IsArray<T>() && ArrayHasStaticExtents<T>())> constexpr bool Helper(
  integer_sequence<long long, BeginElement1, BeginElement2, RemainingElements...>, index_sequence<
  Dim1, Dim2, RemainingDims...>) {
  return array_traits<T>::template ExtentBegin<Dim1>() == BeginElement1 && Helper<T>(
    integer_sequence<long long, BeginElement2, RemainingElements...>(), index_sequence<Dim2,
    RemainingDims...>());
}
template <typename T, typename BeginElementSequence, typename DimSequence, OVK_FUNCTION_REQUIRES(
  !IsArray<T>() || ArrayHasRuntimeExtents<T>())> constexpr bool Helper(BeginElementSequence,
  DimSequence) {
  return false;
}
}
template <typename T, long long... BeginElements> constexpr bool StaticArrayHasExtentsBegin() {
  return static_array_has_extents_begin_internal::Helper<T>(integer_sequence<long long,
    BeginElements...>(), index_sequence_of_size<sizeof...(BeginElements)>());
}

namespace static_array_has_extents_end_internal {
template <typename T, long long EndElement, std::size_t Dim, OVK_FUNCTION_REQUIRES(IsArray<T>() &&
  ArrayHasStaticExtents<T>())> constexpr bool Helper(integer_sequence<long long, EndElement>,
  index_sequence<Dim>) {
  return array_traits<T>::template ExtentEnd<Dim>() == EndElement;
}
template <typename T, long long EndElement1, long long EndElement2, long long...  RemainingElements,
  std::size_t Dim1, std::size_t Dim2, std::size_t... RemainingDims, OVK_FUNCTION_REQUIRES(
  IsArray<T>() && ArrayHasStaticExtents<T>())> constexpr bool Helper(integer_sequence<long long,
  EndElement1, EndElement2, RemainingElements...>, index_sequence<Dim1, Dim2, RemainingDims...
  >) {
  return array_traits<T>::template ExtentEnd<Dim1>() == EndElement1 && Helper<T>(integer_sequence<
    long long, EndElement2, RemainingElements...>(), index_sequence<Dim2, RemainingDims...>());
}
template <typename T, typename EndElementSequence, typename DimSequence, OVK_FUNCTION_REQUIRES(
  !IsArray<T>() || ArrayHasRuntimeExtents<T>())> constexpr bool Helper(EndElementSequence,
    DimSequence) {
  return false;
}
}
template <typename T, long long... EndElements> constexpr bool StaticArrayHasExtentsEnd() {
  return static_array_has_extents_end_internal::Helper<T>(integer_sequence<long long,
    EndElements...>(), index_sequence_of_size<sizeof...(EndElements)>());
}

namespace static_array_has_size_internal {
template <typename T, long long SizeElement, std::size_t Dim, OVK_FUNCTION_REQUIRES(IsArray<T>()
  && ArrayHasStaticExtents<T>())> constexpr bool Helper(integer_sequence<long long, SizeElement>,
  index_sequence<Dim>) {
  return (array_traits<T>::template ExtentEnd<Dim>() - array_traits<T>::template ExtentBegin<
    Dim>()) == SizeElement;
}
template <typename T, long long SizeElement1, long long SizeElement2, long long...
  RemainingElements, std::size_t Dim1, std::size_t Dim2, std::size_t... RemainingDims,
  OVK_FUNCTION_REQUIRES(IsArray<T>() && ArrayHasStaticExtents<T>())> constexpr bool Helper(
  integer_sequence<long long, SizeElement1, SizeElement2, RemainingElements...>, index_sequence<
  Dim1, Dim2, RemainingDims...>) {
  return (array_traits<T>::template ExtentEnd<Dim1>() - array_traits<T>::template ExtentBegin<
    Dim1>()) == SizeElement1 && Helper<T>(integer_sequence<long long, SizeElement2,
    RemainingElements...>(), index_sequence<Dim2, RemainingDims...>());
}
template <typename T, typename SizeElementSequence, typename DimSequence, OVK_FUNCTION_REQUIRES(
  !IsArray<T>() || ArrayHasRuntimeExtents<T>())> constexpr bool Helper(SizeElementSequence,
  DimSequence) {
  return false;
}
}
template <typename T, long long... SizeElements> constexpr bool StaticArrayHasSize() {
  return static_array_has_size_internal::Helper<T>(integer_sequence<long long, SizeElements...>(),
    index_sequence_of_size<sizeof...(SizeElements)>());
}

}

}

#endif
