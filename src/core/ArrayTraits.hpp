// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ARRAY_TRAITS_HPP_INCLUDED
#define OVK_CORE_ARRAY_TRAITS_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IntegerSequence.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <vector>
#include <type_traits>

namespace ovk {

template <typename T, typename=void> struct array_traits {
/*
  using value_type = ???;

  static constexpr int Rank = ???;

  static constexpr const array_layout Layout = ???;

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

// C-style array
template <typename T> struct array_traits<T, OVK_SPECIALIZATION_REQUIRES(std::is_array<T>::value)> {
  using value_type = typename std::remove_all_extents<T>::type;
  static constexpr int Rank = std::rank<T>::value;
  static constexpr const array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static constexpr long long Begin() { return 0; }
  template <int iDim> static constexpr long long End() { return std::extent<T,iDim>::value; }
  // Not sure if there's a better way to do this that works for general multidimensional arrays
  static const value_type *Data(const T &Array) {
    return reinterpret_cast<const value_type *>(&Array[0]);
  }
  static value_type *Data(T &Array) {
    return reinterpret_cast<value_type *>(&Array[0]);
  }
};

// std::array
template <typename T, std::size_t N> struct array_traits<std::array<T,N>> {
  using value_type = T;
  static constexpr int Rank = 1;
  static constexpr const array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static constexpr long long Begin() { return 0; }
  template <int> static constexpr long long End() { return N; }
  static const T *Data(const std::array<T,N> &Array) { return Array.data(); }
  static T *Data(std::array<T,N> &Array) { return Array.data(); }
};

// std::vector
template <typename T, typename Allocator> struct array_traits<std::vector<T, Allocator>> {
  using value_type = T;
  static constexpr int Rank = 1;
  static constexpr const array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static long long Begin(const std::vector<T, Allocator> &) { return 0; }
  template <int> static long long End(const std::vector<T, Allocator> &Vec) { return Vec.size(); }
  static const T *Data(const std::vector<T, Allocator> &Vec) { return Vec.data(); }
  static T *Data(std::vector<T, Allocator> &Vec) { return Vec.data(); }
};

// std::basic_string
template <typename CharT, typename Traits, typename Allocator> struct array_traits<
  std::basic_string<CharT, Traits, Allocator>> {
  using value_type = CharT;
  static constexpr int Rank = 1;
  static constexpr const array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static long long Begin(const std::basic_string<CharT, Traits, Allocator> &) {
    return 0;
  }
  template <int> static long long End(const std::basic_string<CharT, Traits, Allocator> &String) {
    return String.length();
  }
  static const CharT *Data(const std::basic_string<CharT, Traits, Allocator> &String) {
    return String.c_str();
  }
  // No non-const Data access
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

namespace array_traits_internal {
template <typename T, typename=void> struct value_type_helper;
template <typename T> struct value_type_helper<T, OVK_SPECIALIZATION_REQUIRES(IsArray<T>())> {
  using type = typename array_traits<T>::value_type;
};
template <typename T> struct value_type_helper<T, OVK_SPECIALIZATION_REQUIRES(!IsArray<T>())> {
  using type = std::false_type;
};
}
template <typename T> using array_value_type = typename array_traits_internal::value_type_helper<T>
  ::type;

template <typename T, OVK_FUNCTION_REQUIRES(IsArray<T>())> constexpr int ArrayRank() {
  return array_traits<T>::Rank;
}
template <typename T, OVK_FUNCTION_REQUIRES(!IsArray<T>())> constexpr int ArrayRank() {
  return 1;
}

template <typename T, OVK_FUNCTION_REQUIRES(IsArray<T>())> constexpr array_layout ArrayLayout() {
  return array_traits<T>::Layout;
}
template <typename T, OVK_FUNCTION_REQUIRES(!IsArray<T>())> constexpr array_layout ArrayLayout() {
  return array_layout::ROW_MAJOR;
}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(IsArray<T>())>
  constexpr bool ArrayHasFootprint() {
  return ArrayRank<T>() == Rank && (Rank == 1 || ArrayLayout<T>() == Layout);
}
template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(!IsArray<T>())>
  constexpr bool ArrayHasFootprint() {
  return false;
}

template <typename T, typename U, OVK_FUNCTION_REQUIRES(IsArray<T>() && IsArray<U>())>
  constexpr bool ArraysAreSimilar() {
  return ArrayRank<T>() == ArrayRank<U>() && (ArrayLayout<T>() == 1 || ArrayLayout<T>() ==
    ArrayLayout<U>());
}
template <typename T, typename U, OVK_FUNCTION_REQUIRES(!IsArray<T>() || !IsArray<U>())>
  constexpr bool ArraysAreSimilar() {
  return false;
}

namespace array_traits_internal {
template <typename TRef, typename=void> struct access_type_helper;
template <typename TRef> struct access_type_helper<TRef, OVK_SPECIALIZATION_REQUIRES(IsArray<
  remove_cvref<TRef>>())> {
  using type = mimicked_ref<TRef, typename array_traits<remove_cvref<TRef>>::value_type>;
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
