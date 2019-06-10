// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ARRAY_VIEW_HPP_INCLUDED
#define OVK_CORE_ARRAY_VIEW_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/IntegerSequence.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeSequence.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <utility>

namespace ovk {

template <typename T, int Rank, array_layout Layout> class array_view;

namespace array_view_internal {
template <typename ArrayRefType> using array_view_type = array_view<typename std::remove_reference<
  core::mimic_cvref<ArrayRefType, core::array_value_type<core::remove_cvref<ArrayRefType>>>
  >::type, core::ArrayRank<core::remove_cvref<ArrayRefType>>(), core::ArrayLayout<
  core::remove_cvref<ArrayRefType>>()>;
}
template <typename ArrayRefType, OVK_FUNCDECL_REQUIRES(core::IsArray<core::remove_cvref<
  ArrayRefType>>())> array_view_internal::array_view_type<ArrayRefType &&> MakeArrayView(
  ArrayRefType &&Array);
template <typename ArrayRefType, OVK_FUNCDECL_REQUIRES(core::IsArray<core::remove_cvref<
  ArrayRefType>>())> array_view_internal::array_view_type<ArrayRefType &&> MakeArrayView(
  ArrayRefType &&Array, const typename array_view_internal::array_view_type<ArrayRefType &&>::
  interval_type &Extents);

namespace core {

template <typename PtrType, typename T> constexpr bool PtrIsViewableAs() {
  return std::is_convertible<PtrType, T *>::value;
}

template <typename TRef, typename U, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(
  IsArray<remove_cvref<TRef>>())> constexpr bool ArrayIsViewableAs() {
  using array_type = remove_cvref<TRef>;
  return ArrayHasFootprint<array_type, Rank, Layout>() && std::is_convertible<typename
    std::remove_reference<mimic_cvref<TRef, array_value_type<array_type>>>::type *, U *>::value;
}
template <typename TRef, typename U, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(
  !core::IsArray<remove_cvref<TRef>>())> constexpr bool ArrayIsViewableAs() {
  return false;
}

}

namespace array_view_internal {

template <typename T, int Rank, array_layout Layout, typename TupleElementTypeSequence> class
  array_view_base_1;

template <typename T, int Rank, array_layout Layout, typename... TupleElementTypes>
  class array_view_base_1<T, Rank, Layout, core::type_sequence<TupleElementTypes...>> {

public:

  using value_type = T;
  using index_type = long long;
  using tuple_element_type = long long;
  using interval_type = interval<tuple_element_type,Rank>;
  using indexer_type = indexer<index_type, tuple_element_type, Rank, Layout>;

protected:

  value_type *Ptr_;
  interval_type Extents_;
  indexer_type Indexer_;
  index_type NumValues_;

public:

  constexpr array_view_base_1():
    Ptr_(nullptr),
    Extents_(MakeEmptyInterval<tuple_element_type,Rank>()),
    Indexer_(Extents_),
    NumValues_(0)
  {}

  constexpr array_view_base_1(value_type *Ptr, const interval_type &Extents):
    Ptr_(Ptr),
    Extents_(Extents),
    Indexer_(Extents),
    NumValues_(Extents.Count())
  {}

  constexpr OVK_FORCE_INLINE value_type &operator()(TupleElementTypes... TupleElements) const {
    return Ptr_[Indexer_.ToIndex(TupleElements...)];
  }

  constexpr OVK_FORCE_INLINE value_type *Data(TupleElementTypes... TupleElements) const {
    return Ptr_+Indexer_.ToIndex(TupleElements...);
  }

};

template <typename T, int Rank, array_layout Layout, typename TupleElementTypeSequence,
  typename=void> class array_view_base_2;

template <typename T, int Rank, array_layout Layout, typename TupleElementTypeSequence> class
  array_view_base_2<T, Rank, Layout, TupleElementTypeSequence, OVK_SPECIALIZATION_REQUIRES(
  !std::is_const<T>::value)> : public array_view_base_1<T, Rank, Layout, TupleElementTypeSequence>
  {

  using parent_type = array_view_base_1<T, Rank, Layout, TupleElementTypeSequence>; 

protected:

  using parent_type::Ptr_;
  using parent_type::Extents_;
  using parent_type::Indexer_;
  using parent_type::NumValues_;

public:

  using value_type = T;
  using index_type = long long;

  using parent_type::parent_type;
  using parent_type::operator();
  using parent_type::Data;

  const array_view_base_2 &Fill(const value_type &Value) const {
    for (index_type i = 0; i < NumValues_; ++i) {
      Ptr_[i] = Value;
    }
    return *this;
  }

  const array_view_base_2 &Fill(std::initializer_list<value_type> ValuesList) const {
    auto Iter = ValuesList.begin();
    for (index_type i = 0; i < NumValues_; ++i) {
      Ptr_[i] = *Iter++;
    }
    return *this;
  }

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsInputIterator<IterType>() &&
    std::is_convertible<core::iterator_reference_type<IterType>, value_type>::value)> const
    array_view_base_2 &Fill(IterType First) const {
    IterType Iter = First;
    for (index_type i = 0; i < NumValues_; ++i) {
      Ptr_[i] = *Iter++;
    }
    return *this;
  }

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U
    >::type, value_type>::value)> const array_view_base_2 &Fill(array_view<U, Rank, Layout>
    SourceView) const {
    for (index_type i = 0; i < NumValues_; ++i) {
      Ptr_[i] = SourceView[i];
    }
    return *this;
  }

  template <typename SourceArrayRefType, OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<
    SourceArrayRefType>>() && !core::IsIterator<core::remove_cvref<SourceArrayRefType>>() &&
    core::ArrayHasFootprint<core::remove_cvref<SourceArrayRefType>, Rank, Layout>() &&
    std::is_convertible<core::array_access_type<SourceArrayRefType &&>, value_type>::value)> const
    array_view_base_2 &Fill(SourceArrayRefType &&SourceArray) const {
    for (index_type i = 0; i < NumValues_; ++i) {
      Ptr_[i] = static_cast<core::array_access_type<SourceArrayRefType &&>>(SourceArray[i]);
    }
    return *this;
  }

};

template <typename T, int Rank, array_layout Layout, typename TupleElementTypeSequence> class
  array_view_base_2<T, Rank, Layout, TupleElementTypeSequence, OVK_SPECIALIZATION_REQUIRES(
  std::is_const<T>::value)> : public array_view_base_1<T, Rank, Layout, TupleElementTypeSequence>
  {

  using parent_type = array_view_base_1<T, Rank, Layout, TupleElementTypeSequence>;

protected:

  using parent_type::Ptr_;
  using parent_type::Extents_;
  using parent_type::Indexer_;
  using parent_type::NumValues_;

public:

  using parent_type::parent_type;
  using parent_type::operator();
  using parent_type::Data;

  template <typename... Args> const array_view_base_2 &Fill(Args &&...) const {
    // Assert on Args instead of just 'false' so it doesn't trigger unless method is instantiated
    static_assert(int(sizeof...(Args)) < 0, "Cannot use array_view::Fill when value type "
      "is constant.");
    return *this;
  }

};

}

template <typename T, int Rank_=1, array_layout Layout_=array_layout::ROW_MAJOR> class array_view :
  private array_view_internal::array_view_base_2<T, Rank_, Layout_,
  core::repeated_type_sequence_of_size<long long,Rank_>> {

  using parent_type = array_view_internal::array_view_base_2<T, Rank_, Layout_,
    core::repeated_type_sequence_of_size<long long,Rank_>>;

  using parent_type::Ptr_;
  using parent_type::Extents_;
  using parent_type::Indexer_;
  using parent_type::NumValues_;

public:

  using value_type = T;
  static constexpr const int Rank = Rank_;
  static constexpr const array_layout Layout = Layout_;
  using index_type = long long;
  using tuple_element_type = long long;
  using tuple_type = elem<tuple_element_type,Rank>;
  using interval_type = interval<tuple_element_type,Rank>;
  using indexer_type = indexer<index_type, tuple_element_type, Rank, Layout>;
  using iterator = value_type *;

  using parent_type::operator();
  using parent_type::Data;
  using parent_type::Fill;

  constexpr array_view() = default;

  constexpr array_view(value_type *Ptr, const interval_type &Extents):
    parent_type(Ptr, Extents)
  {}

  // Intel 17 didn't like using Rank and Layout instead of Rank_ and Layout_
  template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(core::ArrayIsViewableAs<ArrayRefType &&,
    value_type, Rank_, Layout_>())> constexpr array_view(ArrayRefType &&Array):
    parent_type(core::ArrayData(Array), core::ArrayExtents(Array))
  {}

  // Intel 17 didn't like using Rank and Layout instead of Rank_ and Layout_
  template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(core::ArrayIsViewableAs<ArrayRefType &&,
    value_type, Rank_, Layout_>())> constexpr array_view(ArrayRefType &&Array, const
    interval_type &Extents):
    parent_type(core::ArrayData(Array), Extents)
  {}

  template <typename U, OVK_FUNCTION_REQUIRES(!std::is_same<U, value_type>::value &&
    core::PtrIsViewableAs<U *, value_type>())> constexpr array_view(const array_view<U, Rank,
    Layout> &Other):
    parent_type(Other.Ptr_, Other.Extents_)
  {}

  template <typename U, OVK_FUNCTION_REQUIRES(!std::is_same<U, value_type>::value &&
    core::PtrIsViewableAs<U *, value_type>())> constexpr array_view(const array_view<U, Rank,
    Layout> &Other, const interval_type &Extents):
    parent_type(Other.Ptr_, Extents)
  {}

  constexpr array_view(const array_view &Other) = default;

  constexpr array_view(const array_view &Other, const interval_type &Extents):
    parent_type(Other.Ptr_, Extents)
  {}

  constexpr array_view(array_view &&Other) noexcept = default;

  array_view &operator=(const array_view &Other) = default;
  array_view &operator=(array_view &&Other) noexcept = default;

  void Reset() { *this = array_view(); }

  explicit constexpr OVK_FORCE_INLINE operator bool() const { return Ptr_; }

  constexpr OVK_FORCE_INLINE const interval_type &Extents() const { return Extents_; }

  constexpr OVK_FORCE_INLINE const tuple_type &Begin() const { return Extents_.Begin(); }
  constexpr OVK_FORCE_INLINE tuple_element_type Begin(int iDim) const {
    return Extents_.Begin(iDim);
  }

  constexpr OVK_FORCE_INLINE const tuple_type &End() const { return Extents_.End(); }
  constexpr OVK_FORCE_INLINE tuple_element_type End(int iDim) const {
    return Extents_.End(iDim);
  }

  constexpr OVK_FORCE_INLINE tuple_type Size() const { return Extents_.Size(); }
  constexpr OVK_FORCE_INLINE tuple_element_type Size(int iDim) const {
    return Extents_.Size(iDim);
  }

  constexpr OVK_FORCE_INLINE index_type Count() const { return NumValues_; }

  constexpr OVK_FORCE_INLINE bool Empty() const { return Extents_.Empty(); }

  constexpr OVK_FORCE_INLINE const indexer_type &Indexer() const { return Indexer_; }

  // Want to use iterator directly here instead of constructing intermediate elem type
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsRandomAccessIterator<IterType>() &&
    std::is_convertible<core::iterator_reference_type<IterType>, tuple_element_type>::value)>
    constexpr OVK_FORCE_INLINE value_type &operator()(IterType First) const {
    return Ptr_[Indexer_.ToIndex(First)];
  }

  // Want to use array directly here instead of constructing intermediate elem type
  // Intel 17 didn't like using Rank instead of Rank_
  template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
    !core::IsIterator<typename std::decay<ArrayType>::type>() && std::is_convertible<
    core::array_access_type<const ArrayType &>, tuple_element_type>::value && core::ArrayRank<
    ArrayType>() == 1 && (core::ArrayHasRuntimeExtents<ArrayType>() || (core::StaticArrayHasBegin<
    ArrayType,0>() && core::StaticArrayHasEnd<ArrayType,Rank_>())))> constexpr OVK_FORCE_INLINE
    value_type &operator()(const ArrayType &Array) const {
    return Ptr_[Indexer_.ToIndex(Array)];
  }

  constexpr OVK_FORCE_INLINE value_type &operator[](int Index) const { return Ptr_[Index]; }

  constexpr OVK_FORCE_INLINE value_type *Data() const { return Ptr_; }

  // Want to use iterator directly here instead of constructing intermediate elem type
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsRandomAccessIterator<IterType>() &&
    std::is_convertible<core::iterator_reference_type<IterType>, tuple_element_type>::value)>
    constexpr OVK_FORCE_INLINE value_type *Data(IterType First) const {
    return Ptr_+Indexer_.ToIndex(First);
  }

  // Want to use array directly here instead of constructing intermediate elem type
  // Intel 17 didn't like using Rank instead of Rank_
  template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
    !core::IsIterator<typename std::decay<ArrayType>::type>() && std::is_convertible<
    core::array_access_type<const ArrayType &>, tuple_element_type>::value && core::ArrayRank<
    ArrayType>() == 1 && (core::ArrayHasRuntimeExtents<ArrayType>() || (core::StaticArrayHasBegin<
    ArrayType,0>() && core::StaticArrayHasEnd<ArrayType,Rank_>())))> constexpr OVK_FORCE_INLINE
    value_type *Data(const ArrayType &Array) const {
    return Ptr_+Indexer_.ToIndex(Array);
  }

  constexpr OVK_FORCE_INLINE iterator LinearBegin() const { return Ptr_; }
  constexpr OVK_FORCE_INLINE iterator LinearEnd() const { return Ptr_ + NumValues_; }

  // Google Test doesn't use free begin/end functions and instead expects container to have
  // lowercase begin/end methods
  constexpr OVK_FORCE_INLINE iterator begin() const { return LinearBegin(); }
  constexpr OVK_FORCE_INLINE iterator end() const { return LinearEnd(); }

private:

  template <typename U, int OtherRank, array_layout OtherLayout> friend class array_view;

  friend class core::test_helper<array_view>;

};

template <typename T, int Rank, array_layout Layout> constexpr OVK_FORCE_INLINE typename
  array_view<T, Rank, Layout>::iterator begin(const array_view<T, Rank, Layout> &View) {
  return View.LinearBegin();
}

template <typename T, int Rank, array_layout Layout> constexpr OVK_FORCE_INLINE typename
  array_view<T, Rank, Layout>::iterator end(const array_view<T, Rank, Layout> &View) {
  return View.LinearEnd();
}

template <typename ArrayRefType, OVK_FUNCDEF_REQUIRES(core::IsArray<core::remove_cvref<
  ArrayRefType>>())> array_view_internal::array_view_type<ArrayRefType &&> MakeArrayView(
  ArrayRefType &&Array) {
  return {Array};
}

template <typename ArrayRefType, OVK_FUNCDEF_REQUIRES(core::IsArray<core::remove_cvref<
  ArrayRefType>>())> array_view_internal::array_view_type<ArrayRefType &&> MakeArrayView(
  ArrayRefType &&Array, const typename array_view_internal::array_view_type<ArrayRefType &&>::
  interval_type &Extents) {
  return {Array, Extents};
}

template <typename T, typename U, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(std::is_same<
  typename std::remove_const<T>::type, typename std::remove_const<U>::type>::value)> constexpr
  bool operator==(const array_view<T, Rank, Layout> &Left, const array_view<U, Rank, Layout>
  &Right) {
  return Left.Data() == Right.Data() && Left.Extents() == Right.Extents();
}
template <typename T, typename U, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(std::is_same<
  typename std::remove_const<T>::type, typename std::remove_const<U>::type>::value)> constexpr
  bool operator!=(const array_view<T, Rank, Layout> &Left, const array_view<U, Rank, Layout>
  &Right) {
  return !(Left == Right);
}

template <typename T, int Rank, array_layout Layout> constexpr bool operator==(const array_view<T,
  Rank, Layout> &Left, std::nullptr_t) {
  return Left.Data() == nullptr;
}
template <typename T, int Rank, array_layout Layout> constexpr bool operator==(std::nullptr_t,
  const array_view<T, Rank, Layout> &Right) {
  return nullptr == Right.Data();
}
template <typename T, int Rank, array_layout Layout> constexpr bool operator!=(const array_view<T,
  Rank, Layout> &Left, std::nullptr_t) {
  return !(Left == nullptr);
}
template <typename T, int Rank, array_layout Layout> constexpr bool operator!=(std::nullptr_t,
  const array_view<T, Rank, Layout> &Right) {
  return !(nullptr == Right);
}

}

#endif
