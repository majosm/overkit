// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ARRAY_TRAITS_HPP_INCLUDED
#define OVK_CORE_ARRAY_TRAITS_HPP_INCLUDED

#include <ovk/core/ArrayTraitsBase.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IntegerSequence.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ovk {

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
template <typename ArrayType, std::size_t... Indices> constexpr elem<long long,ArrayRank<
  ArrayType>()> StaticArrayBeginHelper(core::index_sequence<Indices...>) {
  return {array_traits<ArrayType>::template Begin<Indices>()...};
}
template <typename ArrayType, std::size_t... Indices> elem<long long,ArrayRank<ArrayType>()>
  RuntimeArrayBeginHelper(core::index_sequence<Indices...>, const ArrayType &Array) {
  return {array_traits<ArrayType>::template Begin<Indices>(Array)...};
}
}
template <typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>() &&
  ArrayHasStaticExtents<ArrayType>())> elem<long long,ArrayRank<ArrayType>()> constexpr ArrayBegin(
  const ArrayType &) {
  return array_traits_internal::StaticArrayBeginHelper<ArrayType>(core::index_sequence_of_size<
    ArrayRank<ArrayType>()>());
}
template <typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>() &&
  ArrayHasRuntimeExtents<ArrayType>())> elem<long long,ArrayRank<ArrayType>()> ArrayBegin(const
  ArrayType &Array) {
  return array_traits_internal::RuntimeArrayBeginHelper<ArrayType>(core::index_sequence_of_size<
    ArrayRank<ArrayType>()>(), Array);
}

namespace array_traits_internal {
template <typename ArrayType, std::size_t... Indices> constexpr elem<long long,ArrayRank<
  ArrayType>()> StaticArrayEndHelper(core::index_sequence<Indices...>) {
  return {array_traits<ArrayType>::template End<Indices>()...};
}
template <typename ArrayType, std::size_t... Indices> elem<long long,ArrayRank<ArrayType>()>
  RuntimeArrayEndHelper(core::index_sequence<Indices...>, const ArrayType &Array) {
  return {array_traits<ArrayType>::template End<Indices>(Array)...};
}
}
template <typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>() &&
  ArrayHasStaticExtents<ArrayType>())> elem<long long,ArrayRank<ArrayType>()> constexpr ArrayEnd(
  const ArrayType &) {
  return array_traits_internal::StaticArrayEndHelper<ArrayType>(core::index_sequence_of_size<
    ArrayRank<ArrayType>()>());
}
template <typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>() &&
  ArrayHasRuntimeExtents<ArrayType>())> elem<long long,ArrayRank<ArrayType>()> ArrayEnd(const
  ArrayType &Array) {
  return array_traits_internal::RuntimeArrayEndHelper<ArrayType>(core::index_sequence_of_size<
    ArrayRank<ArrayType>()>(), Array);
}

namespace array_traits_internal {
template <typename ArrayType, std::size_t... Indices> constexpr elem<long long,ArrayRank<
  ArrayType>()> StaticArraySizeHelper(core::index_sequence<Indices...>) {
  return {array_traits<ArrayType>::template End<Indices>() - array_traits<ArrayType>::template
  Begin<Indices>()...};
}
template <typename ArrayType, std::size_t... Indices> elem<long long,ArrayRank<ArrayType>()>
  RuntimeArraySizeHelper(core::index_sequence<Indices...>, const ArrayType &Array) {
  return {array_traits<ArrayType>::template End<Indices>(Array) - array_traits<ArrayType>::template
    Begin<Indices>(Array)...};
}
}
template <typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>() &&
  ArrayHasStaticExtents<ArrayType>())> elem<long long,ArrayRank<ArrayType>()> constexpr ArraySize(
  const ArrayType &) {
  return array_traits_internal::StaticArraySizeHelper<ArrayType>(core::index_sequence_of_size<
    ArrayRank<ArrayType>()>());
}
template <typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>() &&
  ArrayHasRuntimeExtents<ArrayType>())> elem<long long,ArrayRank<ArrayType>()> ArraySize(const
  ArrayType &Array) {
  return array_traits_internal::RuntimeArraySizeHelper<ArrayType>(core::index_sequence_of_size<
    ArrayRank<ArrayType>()>(), Array);
}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>())> interval<long long,
  ArrayRank<ArrayType>()> ArrayExtents(const ArrayType &Array) {
  return {ArrayBegin(Array), ArrayEnd(Array)};
}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<ArrayType>())> long long ArrayCount(
  const ArrayType &Array) {
  return interval<long long,ArrayRank<ArrayType>>(ArraySize(Array)).Count();
}

template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(IsArray<remove_cvref<ArrayRefType>>())> auto
  ArrayData(ArrayRefType &&Array) -> decltype(array_traits<remove_cvref<ArrayRefType>>::Data(
  std::forward<ArrayRefType>(Array))) {
  return array_traits<remove_cvref<ArrayRefType>>::Data(std::forward<ArrayRefType>(Array));
}

template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(IsArray<remove_cvref<ArrayRefType>>())> auto
  ArrayLinearBegin(ArrayRefType &&Array) -> decltype(MakeForwardingIterator<array_access_type<
  ArrayRefType &&>>(ArrayData(Array))) {
  return MakeForwardingIterator<array_access_type<ArrayRefType &&>>(ArrayData(Array));
}

template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(IsArray<remove_cvref<ArrayRefType>>())> auto
  ArrayLinearEnd(ArrayRefType &&Array) -> decltype(MakeForwardingIterator<array_access_type<
  ArrayRefType &&>>(core::ArrayData(Array)+ArrayCount(Array))) {
  return MakeForwardingIterator<array_access_type<ArrayRefType &&>>(ArrayData(Array)+
    ArrayCount(Array));
}

}

}

#endif
