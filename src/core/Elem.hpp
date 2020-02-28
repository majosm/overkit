// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ELEM_HPP_INCLUDED
#define OVK_CORE_ELEM_HPP_INCLUDED

#include <ovk/core/ArrayTraitsBase.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IntegerSequence.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/ScalarTraits.hpp>
#include <ovk/core/TypeSequence.hpp>

#include <type_traits>

namespace ovk {

namespace elem_internal {

// Want elem<T,1> to be castable to T
namespace cast_internal {
struct null_type {};
template <typename T, int N> struct scalar_cast_helper {
  using const_ref = null_type;
  using ref = null_type;
};
template <typename T> struct scalar_cast_helper<T,1> {
  using const_ref = const T &;
  using ref = T &;
};
template <typename T, int N> using scalar_cast_const_ref_type = typename scalar_cast_helper<T,N>::
  const_ref;
template <typename T, int N> using scalar_cast_ref_type = typename scalar_cast_helper<T,N>::ref;
template <typename T, int N> constexpr scalar_cast_const_ref_type<T,N> ScalarCast(const T
  (&Values)[N]) {
  return {};
}
template <typename T, int N> constexpr scalar_cast_ref_type<T,N> ScalarCast(T (&Values)[N]) {
  return {};
}
template <typename T> constexpr scalar_cast_const_ref_type<T,1> ScalarCast(const T (&Values)[1]) {
  return Values[0];
}
template <typename T> constexpr scalar_cast_ref_type<T,1> ScalarCast(T (&Values)[1]) {
  return Values[0];
}
}

template <typename U, typename T, int N> constexpr bool ImplicitlyConvertibleToElem() {
  return (core::IsRandomAccessIterator<U>() && std::is_convertible<core::iterator_reference_type<U>,
    T>::value) || (core::IsArray<U>() && std::is_convertible<core::array_access_type<const U &>,
    T>::value && core::ArrayRank<U>() == 1 && (core::ArrayHasRuntimeExtents<U>() ||
    (core::StaticArrayHasExtentsBegin<U, 0>() && core::StaticArrayHasExtentsEnd<U, N>())));
}

template <typename U, typename T, int N> constexpr bool ExplicitlyConvertibleToElem() {
  return (core::IsRandomAccessIterator<U>() && !std::is_convertible<core::iterator_reference_type<
    U>, T>::value && std::is_constructible<T, core::iterator_reference_type<U>>::value) ||
    (core::IsArray<U>() && !std::is_convertible<core::array_access_type<const U &>, T>::value &&
    std::is_constructible<T, core::array_access_type<const U &>>::value && core::ArrayRank<U>() == 1
    && (core::ArrayHasRuntimeExtents<U>() || (core::StaticArrayHasExtentsBegin<U, 0>() &&
    core::StaticArrayHasExtentsEnd<U, N>())));
}

template <typename T, int N, typename TSequence> class elem_base_1;

template <typename T, int N, typename... Ts> class elem_base_1<T, N, core::type_sequence<Ts...>> {

private:

  using scalar_cast_const_ref_type = cast_internal::scalar_cast_const_ref_type<T,N>;
  using scalar_cast_ref_type = cast_internal::scalar_cast_ref_type<T,N>;

public:

  static_assert(N > 0 && N <= 127, "Invalid elem rank.");

  using value_type = T;
  static constexpr int Rank = N;
  using iterator = value_type *;
  using const_iterator = const value_type *;

protected:

  value_type Values_[Rank];

public:

  constexpr OVK_FORCE_INLINE elem_base_1() = default;

  constexpr OVK_FORCE_INLINE elem_base_1(Ts... Values):
    Values_{Values...}
  {}

  template <typename U, OVK_FUNCTION_REQUIRES(!std::is_base_of<elem_base_1, U>::value &&
    ImplicitlyConvertibleToElem<U, value_type, Rank>())> constexpr OVK_FORCE_INLINE elem_base_1(
    const U &Other):
    elem_base_1(core::index_sequence_of_size<Rank>(), Other)
  {}

  template <typename U, OVK_FUNCTION_REQUIRES(!std::is_base_of<elem_base_1, U>::value &&
    !ImplicitlyConvertibleToElem<U, value_type, Rank>() && ExplicitlyConvertibleToElem<U,
    value_type, Rank>())> OVK_FORCE_INLINE explicit elem_base_1(const U &Other):
    elem_base_1(core::index_sequence_of_size<Rank>(), Other)
  {}

  constexpr OVK_FORCE_INLINE elem_base_1(const elem_base_1 &Other) = default;

  template <typename U, OVK_FUNCTION_REQUIRES(!std::is_base_of<elem_base_1, U>::value &&
    ImplicitlyConvertibleToElem<U, value_type, Rank>())> OVK_FORCE_INLINE elem_base_1
    &operator=(const U &Other) {
    Assign_(core::index_sequence_of_size<Rank>(), Other);
    return *this;
  }

  OVK_FORCE_INLINE elem_base_1 &operator=(const elem_base_1 &Other) = default;

  constexpr OVK_FORCE_INLINE const value_type &operator()(int iElement) const {
    return Values_[iElement];
  }
  OVK_FORCE_INLINE value_type &operator()(int iElement) { return Values_[iElement]; }

  constexpr OVK_FORCE_INLINE const value_type &operator[](int iElement) const {
    return Values_[iElement];
  }
  OVK_FORCE_INLINE value_type &operator[](int iElement) { return Values_[iElement]; }

  constexpr OVK_FORCE_INLINE const value_type *Data() const { return Values_; }
  OVK_FORCE_INLINE value_type *Data() { return Values_; }

  constexpr OVK_FORCE_INLINE const_iterator Begin() const { return Values_; }
  iterator Begin() { return Values_; }
  constexpr OVK_FORCE_INLINE const_iterator End() const { return Values_+Rank; }
  iterator End() { return Values_+Rank; }

  // Google Test doesn't use free begin/end functions and instead expects container to have
  // lowercase begin/end methods
  constexpr OVK_FORCE_INLINE const_iterator begin() const { return Begin(); }
  iterator begin() { return Begin(); }
  constexpr OVK_FORCE_INLINE const_iterator end() const { return End(); }
  iterator end() { return End(); }

  constexpr operator scalar_cast_const_ref_type() const {
    return cast_internal::ScalarCast(Values_);
  }
  operator scalar_cast_ref_type() {
    return cast_internal::ScalarCast(Values_);
  }

private:

  template <std::size_t... Indices, typename U> constexpr OVK_FORCE_INLINE elem_base_1(
    core::index_sequence<Indices...>, const U &Other):
    Values_{value_type(Other[Indices])...}
  {}

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices, typename U>
    OVK_FORCE_INLINE void Assign_(core::index_sequence<Index1, Index2, RemainingIndices...>,
    const U &Other) {
    Values_[Index1] = value_type(Other[Index1]);
    Assign_(core::index_sequence<Index2, RemainingIndices...>(), Other);
  }

  template <std::size_t Index, typename U> OVK_FORCE_INLINE void Assign_(core::index_sequence<
    Index>, const U &Other) {
    Values_[Index] = value_type(Other[Index]);
  }

};

template <typename T, int N, typename TSequence, typename=void> class elem_base_2;

template <typename T, int N, typename... Ts> class elem_base_2<T, N, core::type_sequence<Ts...>,
  OVK_SPECIALIZATION_REQUIRES(core::IsScalar<T>())> : public elem_base_1<T, N,
  core::type_sequence<Ts...>> {

protected:

  using parent_type = elem_base_1<T, N, core::type_sequence<Ts...>>;
  using parent_type::Values_;

public:

  using typename parent_type::value_type;
  using parent_type::Rank;
  using typename parent_type::iterator;
  using typename parent_type::const_iterator;
  using parent_type::parent_type;

  template <typename U, OVK_FUNCTION_REQUIRES(!std::is_base_of<elem_base_2, U>::value &&
    ImplicitlyConvertibleToElem<U, value_type, Rank>())> elem_base_2 &operator+=(const U &Other) {
    PlusEquals_(core::index_sequence_of_size<Rank>(), Other);
    return *this;
  }

  elem_base_2 &operator+=(const elem_base_2 &Other) {
    PlusEquals_(core::index_sequence_of_size<Rank>(), Other);
    return *this;
  }

  template <typename U, OVK_FUNCTION_REQUIRES(!std::is_base_of<elem_base_2, U>::value &&
    ImplicitlyConvertibleToElem<U, value_type, Rank>())> elem_base_2 &operator-=(const U &Other) {
    MinusEquals_(core::index_sequence_of_size<Rank>(), Other);
    return *this;
  }

  elem_base_2 &operator-=(const elem_base_2 &Other) {
    MinusEquals_(core::index_sequence_of_size<Rank>(), Other);
    return *this;
  }

private:

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices, typename U>
    OVK_FORCE_INLINE void PlusEquals_(core::index_sequence<Index1, Index2, RemainingIndices...>,
    const U &Other) {
    Values_[Index1] += value_type(Other[Index1]);
    PlusEquals_(core::index_sequence<Index2, RemainingIndices...>(), Other);
  }

  template <std::size_t Index, typename U> OVK_FORCE_INLINE void PlusEquals_(core::index_sequence<
    Index>, const U &Other) {
    Values_[Index] += value_type(Other[Index]);
  }

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices, typename U>
    OVK_FORCE_INLINE void MinusEquals_(core::index_sequence<Index1, Index2, RemainingIndices...>,
    const U &Other) {
    Values_[Index1] -= value_type(Other[Index1]);
    MinusEquals_(core::index_sequence<Index2, RemainingIndices...>(), Other);
  }

  template <std::size_t Index, typename U> OVK_FORCE_INLINE void MinusEquals_(core::index_sequence<
    Index>, const U &Other) {
    Values_[Index] -= value_type(Other[Index]);
  }

};

template <typename T, int N, typename... Ts> class elem_base_2<T, N, core::type_sequence<Ts...>,
  OVK_SPECIALIZATION_REQUIRES(!core::IsScalar<T>())> : public elem_base_1<T, N,
  core::type_sequence<Ts...>> {

protected:

  using parent_type = elem_base_1<T, N, core::type_sequence<Ts...>>;
  using parent_type::Values_;

public:

  using typename parent_type::value_type;
  using parent_type::Rank;
  using typename parent_type::iterator;
  using typename parent_type::const_iterator;
  using parent_type::parent_type;

  template <typename... Args> elem_base_2 &operator+=(Args &&...) {
    // Assert on Args instead of just 'false' so it doesn't trigger unless method is instantiated
    static_assert(int(sizeof...(Args)) < 0, "Cannot use elem::operator+= for non-scalar types.");
    return *this;
  }

  template <typename... Args> elem_base_2 &operator-=(Args &&...) {
    // Assert on Args instead of just 'false' so it doesn't trigger unless method is instantiated
    static_assert(int(sizeof...(Args)) < 0, "Cannot use elem::operator-= for non-scalar types.");
    return *this;
  }

};

}

template <typename T, int N=1> class elem : public elem_internal::elem_base_2<T, N,
  core::repeated_type_sequence_of_size<T,N>> {

private:

  using parent_type = elem_internal::elem_base_2<T, N, core::repeated_type_sequence_of_size<T,N>>;
  using parent_type::Values_;

public:

  using typename parent_type::value_type;
  using parent_type::Rank;
  using typename parent_type::iterator;
  using typename parent_type::const_iterator;
  using parent_type::parent_type;
  using parent_type::operator+=;
  using parent_type::operator-=;

private:

  friend class core::test_helper<elem>;

};

template <typename T, int N> struct array_traits<elem<T,N>> {
  using value_type = T;
  static constexpr int Rank = 1;
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static constexpr long long ExtentBegin() { return 0; }
  template <int> static constexpr long long ExtentEnd() { return N; }
  static const T *Data(const elem<T,N> &Elem) { return Elem.Data(); }
  static T *Data(elem<T,N> &Elem) { return Elem.Data(); }
};

template <typename T, int N> constexpr OVK_FORCE_INLINE typename elem<T,N>::iterator begin(
  elem<T,N> &Elem);
template <typename T, int N> constexpr OVK_FORCE_INLINE typename elem<T,N>::const_iterator begin(
  const elem<T,N> &Elem);

template <typename T, int N> constexpr OVK_FORCE_INLINE typename elem<T,N>::iterator end(
  elem<T,N> &Elem);
template <typename T, int N> constexpr OVK_FORCE_INLINE typename elem<T,N>::const_iterator end(
  const elem<T,N> &Elem);

// Can't do this with a constructor due to issue with Intel 16 (nested brace initialization
// ignores explicit keyword)
template <typename T, int N> constexpr OVK_FORCE_INLINE elem<T,N> MakeUniformElem(T Value);

// May appear dangerous since all rank 1 array types are convertible to elem, but templating
// prevents implicit conversions
template <typename T, int N> constexpr OVK_FORCE_INLINE bool operator==(const elem<T,N> &Left,
  const elem<T,N> &Right);
template <typename T, int N> constexpr OVK_FORCE_INLINE bool operator!=(const elem<T,N> &Left,
  const elem<T,N> &Right);

template <typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>())> constexpr OVK_FORCE_INLINE
  elem<T,N> operator+(const elem<T,N> &Left, const elem<T,N> &Right);
template <typename T, int N, typename U, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>() &&
  !std::is_same<U, elem<T,N>>::value && elem_internal::ImplicitlyConvertibleToElem<U, T, N>())>
  constexpr OVK_FORCE_INLINE elem<T,N> operator+(const elem<T,N> &Left, const U &Right);
template <typename U, typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>() &&
  !std::is_same<U, elem<T,N>>::value && elem_internal::ImplicitlyConvertibleToElem<U, T, N>())>
  constexpr OVK_FORCE_INLINE elem<T,N> operator+(const U &Left, const elem<T,N> &Right);

template <typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>())> constexpr OVK_FORCE_INLINE
  elem<T,N> operator-(const elem<T,N> &Left, const elem<T,N> &Right);
template <typename T, int N, typename U, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>() &&
  !std::is_same<U, elem<T,N>>::value && elem_internal::ImplicitlyConvertibleToElem<U, T, N>())>
  constexpr OVK_FORCE_INLINE elem<T,N> operator-(const elem<T,N> &Left, const U &Right);
template <typename U, typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>() &&
  !std::is_same<U, elem<T,N>>::value && elem_internal::ImplicitlyConvertibleToElem<U, T, N>())>
  constexpr OVK_FORCE_INLINE elem<T,N> operator-(const U &Left, const elem<T,N> &Right);

template <typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>())> constexpr OVK_FORCE_INLINE
  elem<T,N> Min(const elem<T,N> &Left, const elem<T,N> &Right);
template <typename T, int N, typename U, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>() &&
  !std::is_same<U, elem<T,N>>::value && elem_internal::ImplicitlyConvertibleToElem<U, T, N>())>
  constexpr OVK_FORCE_INLINE elem<T,N> Min(const elem<T,N> &Left, const U &Right);
template <typename U, typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>() &&
  !std::is_same<U, elem<T,N>>::value && elem_internal::ImplicitlyConvertibleToElem<U, T, N>())>
  constexpr OVK_FORCE_INLINE elem<T,N> Min(const U &Left, const elem<T,N> &Right);

template <typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>())> constexpr OVK_FORCE_INLINE
  elem<T,N> Max(const elem<T,N> &Left, const elem<T,N> &Right);
template <typename T, int N, typename U, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>() &&
  !std::is_same<U, elem<T,N>>::value && elem_internal::ImplicitlyConvertibleToElem<U, T, N>())>
  constexpr OVK_FORCE_INLINE elem<T,N> Max(const elem<T,N> &Left, const U &Right);
template <typename U, typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>() &&
  !std::is_same<U, elem<T,N>>::value && elem_internal::ImplicitlyConvertibleToElem<U, T, N>())>
  constexpr OVK_FORCE_INLINE elem<T,N> Max(const U &Left, const elem<T,N> &Right);

template <typename T, int N1, int N2> constexpr OVK_FORCE_INLINE elem<T,N1+N2> ConcatElems(const
  elem<T,N1> &Left, const elem<T,N2> &Right);

template <typename T, int N, array_layout Layout=array_layout::ROW_MAJOR> class elem_less;

template <typename T, int N> class elem_less<T, N, array_layout::ROW_MAJOR> {

public:

  using value_type = T;
  static constexpr int Rank = N;
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  using elem_type = elem<value_type,Rank>;

  bool operator()(const elem_type &Left, const elem_type &Right) const {
    return Compare_(Left, Right, index_tag<0>());
  }

private:

  template <int iElement> struct index_tag {};

  template <int iElement> OVK_FORCE_INLINE bool Compare_(const elem_type &Left, const elem_type
    &Right, index_tag<iElement>) const {
    return Left(iElement) < Right(iElement) || (Left(iElement) == Right(iElement) && Compare_(
      Left, Right, index_tag<iElement+1>()));
  }

  OVK_FORCE_INLINE bool Compare_(const elem_type &Left, const elem_type &Right, index_tag<Rank-1>)
    const {
    return Left(Rank-1) < Right(Rank-1);
  }

};

template <typename T, int N> class elem_less<T, N, array_layout::COLUMN_MAJOR> {

public:

  using value_type = T;
  static constexpr int Rank = N;
  static constexpr array_layout Layout = array_layout::COLUMN_MAJOR;
  using elem_type = elem<value_type,Rank>;

  bool operator()(const elem_type &Left, const elem_type &Right) const {
    return Compare_(Left, Right, index_tag<Rank-1>());
  }

private:

  template <int iElement> struct index_tag {};

  template <int iElement> OVK_FORCE_INLINE bool Compare_(const elem_type &Left, const elem_type
    &Right, index_tag<iElement>) const {
    return Left(iElement) < Right(iElement) || (Left(iElement) == Right(iElement) && Compare_(
      Left, Right, index_tag<iElement-1>()));
  }

  OVK_FORCE_INLINE bool Compare_(const elem_type &Left, const elem_type &Right, index_tag<0>) const
    {
    return Left(0) < Right(0);
  }

};

template <typename T, int N> using elem_less_r = elem_less<T, N, array_layout::ROW_MAJOR>;
template <typename T, int N> using elem_less_c = elem_less<T, N, array_layout::COLUMN_MAJOR>;

}

#include <ovk/core/Elem.inl>

#endif
