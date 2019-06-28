// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ARRAY_TRAITS_HPP_INCLUDED
#define OVK_CORE_ARRAY_TRAITS_HPP_INCLUDED

#include <ovk/core/ArrayTraitsBase.hpp>
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
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static constexpr long long ExtentBegin() { return 0; }
  template <int iDim> static constexpr long long ExtentEnd() { return std::extent<T,iDim>::value; }
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
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static constexpr long long ExtentBegin() { return 0; }
  template <int> static constexpr long long ExtentEnd() { return N; }
  static const T *Data(const std::array<T,N> &Array) { return Array.data(); }
  static T *Data(std::array<T,N> &Array) { return Array.data(); }
};

// std::vector
template <typename T, typename Allocator> struct array_traits<std::vector<T, Allocator>> {
  using value_type = T;
  static constexpr int Rank = 1;
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static long long ExtentBegin(const std::vector<T, Allocator> &) { return 0; }
  template <int> static long long ExtentEnd(const std::vector<T, Allocator> &Vec) {
    return Vec.size();
  }
  static const T *Data(const std::vector<T, Allocator> &Vec) { return Vec.data(); }
  static T *Data(std::vector<T, Allocator> &Vec) { return Vec.data(); }
};

// std::basic_string
template <typename CharT, typename Traits, typename Allocator> struct array_traits<
  std::basic_string<CharT, Traits, Allocator>> {
  using value_type = CharT;
  static constexpr int Rank = 1;
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static long long ExtentBegin(const std::basic_string<CharT, Traits, Allocator> &) {
    return 0;
  }
  template <int> static long long ExtentEnd(const std::basic_string<CharT, Traits, Allocator>
    &String) {
    return String.length();
  }
  static const CharT *Data(const std::basic_string<CharT, Traits, Allocator> &String) {
    return String.c_str();
  }
  // No non-const Data access
};

namespace core {

namespace array_value_type_internal {
template <typename T, typename=void> struct helper;
template <typename T> struct helper<T, OVK_SPECIALIZATION_REQUIRES(IsArray<T>())> {
  using type = typename array_traits<T>::value_type;
};
template <typename T> struct helper<T, OVK_SPECIALIZATION_REQUIRES(!IsArray<T>())> {
  using type = std::false_type;
};
}
template <typename T> using array_value_type = typename array_value_type_internal::helper<T>::type;

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
  constexpr bool ArraysAreCongruent() {
  return ArrayRank<T>() == ArrayRank<U>() && (ArrayRank<T>() == 1 || ArrayLayout<T>() ==
    ArrayLayout<U>());
}
template <typename T, typename U, OVK_FUNCTION_REQUIRES(!IsArray<T>() || !IsArray<U>())>
  constexpr bool ArraysAreCongruent() {
  return false;
}

namespace array_extents_internal {
template <typename TupleElementType, typename ArrayType, std::size_t... Indices> constexpr
  interval<TupleElementType,ArrayRank<ArrayType>()> StaticHelper(core::index_sequence<Indices...>) {
  return {
    {TupleElementType(array_traits<ArrayType>::template ExtentBegin<Indices>())...},
    {TupleElementType(array_traits<ArrayType>::template ExtentEnd<Indices>())...}
  };
}
template <typename TupleElementType, typename ArrayType, std::size_t... Indices> interval<
  TupleElementType,ArrayRank<ArrayType>()> RuntimeHelper(core::index_sequence<Indices...>, const
  ArrayType &Array) {
  return {
    {TupleElementType(array_traits<ArrayType>::template ExtentBegin<Indices>(Array))...},
    {TupleElementType(array_traits<ArrayType>::template ExtentEnd<Indices>(Array))...}
  };
}
}
template <typename TupleElementType=long long, typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<
  ArrayType>() && ArrayHasStaticExtents<ArrayType>())> constexpr interval<TupleElementType,
  ArrayRank<ArrayType>()> ArrayExtents(const ArrayType &) {
  return array_extents_internal::StaticHelper<TupleElementType, ArrayType>(
    core::index_sequence_of_size<ArrayRank<ArrayType>()>());
}
template <typename TupleElementType=long long, typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<
  ArrayType>() && ArrayHasRuntimeExtents<ArrayType>())> interval<TupleElementType,ArrayRank<
  ArrayType>()> ArrayExtents(const ArrayType &Array) {
  return array_extents_internal::RuntimeHelper<TupleElementType, ArrayType>(
    core::index_sequence_of_size<ArrayRank<ArrayType>()>(), Array);
}

namespace array_size_internal {
template <typename TupleElementType, typename ArrayType, std::size_t... Indices> constexpr
  elem<TupleElementType,ArrayRank<ArrayType>()> StaticHelper(core::index_sequence<Indices...>) {
  return {TupleElementType(array_traits<ArrayType>::template ExtentEnd<Indices>() -
    array_traits<ArrayType>::template ExtentBegin<Indices>())...};
}
template <typename TupleElementType, typename ArrayType, std::size_t... Indices> elem<
  TupleElementType,ArrayRank<ArrayType>()> RuntimeHelper(core::index_sequence<Indices...>, const
  ArrayType &Array) {
  return {TupleElementType(array_traits<ArrayType>::template ExtentEnd<Indices>(Array) -
    array_traits<ArrayType>::template ExtentBegin<Indices>(Array))...};
}
}
template <typename TupleElementType=long long, typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<
  ArrayType>() && ArrayHasStaticExtents<ArrayType>())> constexpr elem<TupleElementType,ArrayRank<
  ArrayType>()> ArraySize(const ArrayType &) {
  return array_size_internal::StaticHelper<TupleElementType, ArrayType>(
    core::index_sequence_of_size<ArrayRank<ArrayType>()>());
}
template <typename TupleElementType=long long, typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<
  ArrayType>() && ArrayHasRuntimeExtents<ArrayType>())> elem<TupleElementType,ArrayRank<ArrayType>()
  > ArraySize(const ArrayType &Array) {
  return array_size_internal::RuntimeHelper<TupleElementType, ArrayType>(
    core::index_sequence_of_size<ArrayRank<ArrayType>()>(), Array);
}

namespace array_count_internal {
template <typename IndexType, typename ArrayType, int Dim, OVK_FUNCTION_REQUIRES(Dim ==
  ArrayRank<ArrayType>()-1)> constexpr IndexType StaticHelper() {
  return IndexType(array_traits<ArrayType>::template ExtentEnd<Dim>() -
    array_traits<ArrayType>::template ExtentBegin<Dim>());
}
template <typename IndexType, typename ArrayType, int Dim, OVK_FUNCTION_REQUIRES(Dim <
  ArrayRank<ArrayType>()-1)> constexpr IndexType StaticHelper() {
  return IndexType(array_traits<ArrayType>::template ExtentEnd<Dim>() - array_traits<ArrayType>::
  template ExtentBegin<Dim>()) * StaticHelper<IndexType, ArrayType, Dim+1>();
}
template <typename IndexType, typename ArrayType, int Dim, OVK_FUNCTION_REQUIRES(Dim ==
  ArrayRank<ArrayType>()-1)> IndexType RuntimeHelper(const ArrayType &Array) {
  return IndexType(array_traits<ArrayType>::template ExtentEnd<Dim>(Array) -
    array_traits<ArrayType>::template ExtentBegin<Dim>(Array));
}
template <typename IndexType, typename ArrayType, int Dim, OVK_FUNCTION_REQUIRES(Dim <
  ArrayRank<ArrayType>()-1)> IndexType RuntimeHelper(const ArrayType &Array) {
  return IndexType(array_traits<ArrayType>::template ExtentEnd<Dim>(Array) - array_traits<
    ArrayType>::template ExtentBegin<Dim>(Array)) * RuntimeHelper<IndexType, ArrayType, Dim+1>(
    Array);
}
}
template <typename IndexType=long long, typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<
  ArrayType>() && ArrayHasStaticExtents<ArrayType>())> constexpr IndexType ArrayCount(const
  ArrayType &) {
  return array_count_internal::StaticHelper<IndexType, ArrayType, 0>();
}
template <typename IndexType=long long, typename ArrayType, OVK_FUNCTION_REQUIRES(IsArray<
  ArrayType>() && ArrayHasRuntimeExtents<ArrayType>())> IndexType ArrayCount(const ArrayType
  &Array) {
  return array_count_internal::RuntimeHelper<IndexType, ArrayType, 0>(Array);
}

template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(IsArray<remove_cvref<ArrayRefType>>())> auto
  ArrayData(ArrayRefType &&Array) -> decltype(array_traits<remove_cvref<ArrayRefType>>::Data(
  std::forward<ArrayRefType>(Array))) {
  return array_traits<remove_cvref<ArrayRefType>>::Data(std::forward<ArrayRefType>(Array));
}

template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(IsArray<remove_cvref<ArrayRefType>>())> auto
  ArrayBegin(ArrayRefType &&Array) -> decltype(MakeForwardingIterator<array_access_type<
  ArrayRefType &&>>(ArrayData(Array))) {
  return MakeForwardingIterator<array_access_type<ArrayRefType &&>>(ArrayData(Array));
}

template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(IsArray<remove_cvref<ArrayRefType>>())> auto
  ArrayEnd(ArrayRefType &&Array) -> decltype(MakeForwardingIterator<array_access_type<
  ArrayRefType &&>>(core::ArrayData(Array)+ArrayCount(Array))) {
  return MakeForwardingIterator<array_access_type<ArrayRefType &&>>(ArrayData(Array)+
    ArrayCount(Array));
}

}

}

#endif
