// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
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

template <typename T, int N, typename TSequence> class elem_base;

template <typename T, int N, typename... Ts> class elem_base<T, N, core::type_sequence<Ts...>> {

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

  constexpr OVK_FORCE_INLINE elem_base() = default;

  constexpr OVK_FORCE_INLINE elem_base(Ts... Values):
    Values_{Values...}
  {}

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsRandomAccessIterator<IterType>() &&
    std::is_convertible<core::iterator_deref_type<IterType>, value_type>::value)> constexpr
    OVK_FORCE_INLINE elem_base(IterType First):
    elem_base(core::index_sequence_of_size<Rank>(), First)
  {}

  template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
    !std::is_base_of<elem_base, ArrayType>::value && !core::IsIterator<typename
    std::decay<ArrayType>::type>() && std::is_convertible<core::array_access_type<const ArrayType
    &>, value_type>::value && core::ArrayRank<ArrayType>() == 1 && (core::ArrayHasRuntimeExtents<
    ArrayType>() || (core::StaticArrayHasBegin<ArrayType, 0>() && core::StaticArrayHasEnd<ArrayType,
    N>())))> constexpr OVK_FORCE_INLINE elem_base(const ArrayType &Array):
    elem_base(core::index_sequence_of_size<Rank>(), Array)
  {}

  template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
    !std::is_base_of<elem_base, ArrayType>::value && !core::IsIterator<typename std::decay<
    ArrayType>::type>() && !std::is_convertible<core::array_access_type<const ArrayType &>,
    value_type>::value && std::is_constructible<value_type, core::array_access_type<const ArrayType
    &>>::value && core::ArrayRank<ArrayType>() == 1 && (core::ArrayHasRuntimeExtents<ArrayType>() ||
    (core::StaticArrayHasBegin<ArrayType, 0>() && core::StaticArrayHasEnd<ArrayType, N>())))>
    constexpr OVK_FORCE_INLINE explicit elem_base(const ArrayType &Array):
    elem_base(core::index_sequence_of_size<Rank>(), Array)
  {}

  constexpr OVK_FORCE_INLINE elem_base(const elem_base &Other) = default;

  OVK_FORCE_INLINE elem_base &operator=(const elem_base &Other) {
    Assign_(core::index_sequence_of_size<Rank>(), Other);
    return *this;
  }

  OVK_FORCE_INLINE elem_base &operator=(elem_base &&Other) noexcept {
    Assign_(core::index_sequence_of_size<Rank>(), Other);
    return *this;
  }

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

  template <std::size_t... Indices, typename ArrayOrIterType> constexpr OVK_FORCE_INLINE
    elem_base(core::index_sequence<Indices...>, const ArrayOrIterType &ArrayOrIter):
    Values_{value_type(ArrayOrIter[Indices])...}
  {}

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices,
    typename ArrayOrIterType> OVK_FORCE_INLINE void Assign_(core::index_sequence<Index1, Index2,
    RemainingIndices...>, const ArrayOrIterType &ArrayOrIter) {
    Values_[Index1] = value_type(ArrayOrIter[Index1]);
    Assign_(core::index_sequence<Index2, RemainingIndices...>(), ArrayOrIter);
  }

  template <std::size_t Index, typename ArrayOrIterType> OVK_FORCE_INLINE void Assign_(
    core::index_sequence<Index>, const ArrayOrIterType &ArrayOrIter) {
    Values_[Index] = value_type(ArrayOrIter[Index]);
  }

};

}

template <typename T, int N=1> class elem : public elem_internal::elem_base<T, N,
  core::repeated_type_sequence_of_size<T,N>> {

private:

  using parent_type = elem_internal::elem_base<T, N, core::repeated_type_sequence_of_size<T,N>>;
  using parent_type::Values_;

public:

  using typename parent_type::value_type;
  using parent_type::Rank;
  using typename parent_type::iterator;
  using typename parent_type::const_iterator;
  using parent_type::parent_type;

private:

  friend class core::test_helper<elem>;

};

template <typename T, int N> struct array_traits<elem<T,N>> {
  using value_type = T;
  static constexpr int Rank = 1;
  static constexpr const array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static constexpr long long Begin() { return 0; }
  template <int> static constexpr long long End() { return N; }
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

template <typename T, int N> constexpr OVK_FORCE_INLINE bool operator==(const elem<T,N> &Left,
  const elem<T,N> &Right);
template <typename T, int N> constexpr OVK_FORCE_INLINE bool operator!=(const elem<T,N> &Left,
  const elem<T,N> &Right);

template <typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>())> OVK_FORCE_INLINE elem<T,N>
  operator+(const elem<T,N> &Left, const elem<T,N> &Right);
template <typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>())> OVK_FORCE_INLINE elem<T,N>
  operator-(const elem<T,N> &Left, const elem<T,N> &Right);

template <typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>())> OVK_FORCE_INLINE elem<T,N>
  Min(const elem<T,N> &Left, const elem<T,N> &Right);
template <typename T, int N, OVK_FUNCDECL_REQUIRES(core::IsScalar<T>())> OVK_FORCE_INLINE elem<T,N>
  Max(const elem<T,N> &Left, const elem<T,N> &Right);

}

#include <ovk/core/Elem.inl>

#endif
