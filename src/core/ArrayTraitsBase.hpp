// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
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
  template <int iDim> static constexpr long long Begin() { return ???; }
  template <int iDim> static constexpr long long End() { return ???; }

  // For runtime-sized arrays
  template <int iDim> static long long Begin(const T &Array) { return ???; }
  template <int iDim> static long long End(const T &Array) { return ???; }

  static const T *Data(const T &Array) { return ???; }
  static T *Data(T &Array) { return ???; }
*/
};

namespace core {

namespace array_traits_internal {
template <typename T> constexpr std::true_type IsArrayTest(typename array_traits<T>::value_type *) {
  return {};
}
template <typename T> constexpr std::false_type IsArrayTest(...) { return {}; }
}
template <typename T> constexpr bool IsArray() {
  return decltype(array_traits_internal::IsArrayTest<T>(nullptr))::value;
}

template <typename T, OVK_FUNCTION_REQUIRES(IsArray<T>())> constexpr int ArrayRank() {
  return array_traits<T>::Rank;
}
template <typename T, OVK_FUNCTION_REQUIRES(!IsArray<T>())> constexpr int ArrayRank() {
  return 1;
}

namespace array_traits_internal {
template <typename TRef, typename=void> struct access_type_helper;
template <typename TRef> struct access_type_helper<TRef, OVK_SPECIALIZATION_REQUIRES(IsArray<
  remove_cvref<TRef>>())> {
  using type = mimic_cvref<TRef, typename array_traits<remove_cvref<TRef>>::value_type>;
};
template <typename TRef> struct access_type_helper<TRef, OVK_SPECIALIZATION_REQUIRES(!IsArray<
  remove_cvref<TRef>>())> {
  using type = std::false_type;
};
}
template <typename TRef> using array_access_type = typename array_traits_internal::
  access_type_helper<TRef>::type;

namespace array_traits_internal {
template <typename T> constexpr std::true_type ArrayHasStaticExtentsTest(decltype(array_traits<T>::
  template Begin<0>()) *) {
  return {};
}
template <typename T> constexpr std::false_type ArrayHasStaticExtentsTest(...) { return {}; }
}
template <typename T> constexpr bool ArrayHasStaticExtents() {
  return decltype(array_traits_internal::ArrayHasStaticExtentsTest<T>(nullptr))::value;
}

namespace array_traits_internal {
template <typename T> constexpr std::true_type ArrayHasRuntimeExtentsTest(decltype(array_traits<T>::
  template Begin<0>(std::declval<T>())) *) {
  return {};
}
template <typename T> constexpr std::false_type ArrayHasRuntimeExtentsTest(...) { return {}; }
}
template <typename T> constexpr bool ArrayHasRuntimeExtents() {
  return decltype(array_traits_internal::ArrayHasRuntimeExtentsTest<T>(nullptr))::value;
}

namespace array_traits_internal {
template <typename T, long long BeginElement, std::size_t Index, OVK_FUNCTION_REQUIRES(
  IsArray<T>() && ArrayHasStaticExtents<T>())> constexpr bool StaticArrayHasBeginHelper(
  integer_sequence<long long, BeginElement>, index_sequence<Index>) {
  return array_traits<T>::template Begin<Index>() == BeginElement;
}
template <typename T, long long BeginElement1, long long BeginElement2, long long...
  RemainingElements, std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices,
  OVK_FUNCTION_REQUIRES(IsArray<T>() && ArrayHasStaticExtents<T>())> constexpr bool
  StaticArrayHasBeginHelper(integer_sequence<long long, BeginElement1, BeginElement2,
  RemainingElements...>, index_sequence<Index1, Index2, RemainingIndices...>) {
  return array_traits<T>::template Begin<Index1>() == BeginElement1 &&
    StaticArrayHasBeginHelper<T>(integer_sequence<long long, BeginElement2,
    RemainingElements...>(), index_sequence<Index2, RemainingIndices...>());
}
template <typename T, typename BeginElementSequence, typename IndexSequence, OVK_FUNCTION_REQUIRES(
  !IsArray<T>() || ArrayHasRuntimeExtents<T>())> constexpr bool StaticArrayHasBeginHelper(
    BeginElementSequence, IndexSequence) {
  return false;
}
}
template <typename T, long long... BeginElements> constexpr bool StaticArrayHasBegin() {
  return array_traits_internal::StaticArrayHasBeginHelper<T>(integer_sequence<
    long long, BeginElements...>(), index_sequence_of_size<sizeof...(BeginElements)>());
}

namespace array_traits_internal {
template <typename T, long long EndElement, std::size_t Index, OVK_FUNCTION_REQUIRES(
  IsArray<T>() && ArrayHasStaticExtents<T>())> constexpr bool StaticArrayHasEndHelper(
  integer_sequence<long long, EndElement>, index_sequence<Index>) {
  return array_traits<T>::template End<Index>() == EndElement;
}
template <typename T, long long EndElement1, long long EndElement2, long long...
  RemainingElements, std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices,
  OVK_FUNCTION_REQUIRES(IsArray<T>() && ArrayHasStaticExtents<T>())> constexpr bool
  StaticArrayHasEndHelper(integer_sequence<long long, EndElement1, EndElement2,
  RemainingElements...>, index_sequence<Index1, Index2, RemainingIndices...>) {
  return array_traits<T>::template End<Index1>() == EndElement1 &&
    StaticArrayHasEndHelper<T>(integer_sequence<long long, EndElement2,
    RemainingElements...>(), index_sequence<Index2, RemainingIndices...>());
}
template <typename T, typename EndElementSequence, typename IndexSequence, OVK_FUNCTION_REQUIRES(
  !IsArray<T>() || ArrayHasRuntimeExtents<T>())> constexpr bool StaticArrayHasEndHelper(
    EndElementSequence, IndexSequence) {
  return false;
}
}
template <typename T, long long... EndElements> constexpr bool StaticArrayHasEnd() {
  return array_traits_internal::StaticArrayHasEndHelper<T>(integer_sequence<
    long long, EndElements...>(), index_sequence_of_size<sizeof...(EndElements)>());
}

namespace array_traits_internal {
template <typename T, long long SizeElement, std::size_t Index, OVK_FUNCTION_REQUIRES(
  IsArray<T>() && ArrayHasStaticExtents<T>())> constexpr bool StaticArrayHasSizeHelper(
  integer_sequence<long long, SizeElement>, index_sequence<Index>) {
  return (array_traits<T>::template End<Index>() - array_traits<T>::template Begin<Index>()) ==
    SizeElement;
}
template <typename T, long long SizeElement1, long long SizeElement2, long long...
  RemainingElements, std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices,
  OVK_FUNCTION_REQUIRES(IsArray<T>() && ArrayHasStaticExtents<T>())> constexpr bool
  StaticArrayHasSizeHelper(integer_sequence<long long, SizeElement1, SizeElement2,
  RemainingElements...>, index_sequence<Index1, Index2, RemainingIndices...>) {
  return (array_traits<T>::template End<Index1>() - array_traits<T>::template Begin<Index1>()) ==
    SizeElement1 && StaticArrayHasSizeHelper<T>(integer_sequence<long long, SizeElement2,
    RemainingElements...>(), index_sequence<Index2, RemainingIndices...>());
}
template <typename T, typename SizeElementSequence, typename IndexSequence, OVK_FUNCTION_REQUIRES(
  !IsArray<T>() || ArrayHasRuntimeExtents<T>())> constexpr bool StaticArrayHasSizeHelper(
    SizeElementSequence, IndexSequence) {
  return false;
}
}
template <typename T, long long... SizeElements> constexpr bool StaticArrayHasSize() {
  return array_traits_internal::StaticArrayHasSizeHelper<T>(integer_sequence<
    long long, SizeElements...>(), index_sequence_of_size<sizeof...(SizeElements)>());
}

}

}

#endif
