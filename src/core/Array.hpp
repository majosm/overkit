// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ARRAY_HPP_INCLUDED
#define OVK_CORE_ARRAY_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/DebugBase.h>
#include <ovk/core/Elem.hpp>
#include <ovk/core/ForEach.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/TypeSequence.hpp>
#include <ovk/core/Vector.hpp>

#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {

namespace array_internal {

template <typename T, int Rank, array_layout Layout, typename TupleElementTypeSequence>
  class array_base_1;

template <typename T, int Rank, array_layout Layout, typename... TupleElementTypes>
  class array_base_1<T, Rank, Layout, core::type_sequence<TupleElementTypes...>> {

public:

  using value_type = T;

protected:

  using view_type = array_view<value_type,Rank,Layout>;

  core::vector<value_type> Values_;
  view_type View_;

public:

  using index_type = long long;
  using tuple_element_type = long long;
  using tuple_type = elem<tuple_element_type,Rank>;
  using interval_type = interval<index_type,Rank>;
  using indexer_type = indexer<index_type, tuple_element_type, Rank, Layout>;
  using iterator = value_type *;
  using const_iterator = const value_type *;

  array_base_1():
    View_(Values_.Data(), MakeEmptyInterval<index_type,Rank>())
  {}

  array_base_1(const interval<index_type,Rank> &Extents):
    Values_(Extents.Count()),
    View_(Values_.Data(), Extents)
  {}

  array_base_1(const interval<index_type,Rank> &Extents, const value_type &Value):
    Values_(Extents.Count(), Value),
    View_(Values_.Data(), Extents)
  {}

  array_base_1(const interval<index_type,Rank> &Extents, std::initializer_list<value_type>
    ValuesList):
    Values_(ValuesList),
    View_(Values_.Data(), Extents)
  {
    OVK_DEBUG_ASSERT((long long)(Values_.Count()) == View_.Count(), "Incorrect number of values in "
      "initializer list");
  }

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsInputIterator<IterType>())>
    array_base_1(const interval<index_type,Rank> &Extents, IterType First) {
    index_type NumValues = Extents.Count();
    Values_.Reserve(NumValues);
    IterType Iter = First;
    for (index_type i = 0; i < NumValues; ++i) {
      Values_.Append(*Iter++);
    }
    View_ = view_type(Values_.Data(), Extents);
  }

  array_base_1(const array_base_1 &Other):
    Values_(Other.Values_),
    View_(Values_.Data(), Other.View_.Extents())
  {}

  array_base_1(const interval_type &Extents, const array_base_1 &Other):
    Values_(Other.Values_),
    View_(Values_.Data(), Extents)
  {}

  array_base_1(array_base_1 &&Other) noexcept:
    Values_(std::move(Other.Values_)),
    View_(Other.View_)
  {
    Other.View_ = view_type();
  }

  array_base_1(const interval_type &Extents, array_base_1 &&Other) noexcept:
    Values_(std::move(Other.Values_)),
    View_(Values_.Data(), Extents)
  {
    Other.View_ = view_type();
  }

  array_base_1 &operator=(const array_base_1 &Other) {
    return Assign(Other);
  }

  array_base_1 &operator=(array_base_1 &&Other) noexcept {
    return Assign(std::move(Other));
  }

  array_base_1 &Assign(const interval_type &Extents, const value_type &Value) {
    Values_.Assign(Extents.Count(), Value);
    View_ = view_type(Values_.Data(), Extents);
    return *this;
  }

  array_base_1 &Assign(const interval_type &Extents, std::initializer_list<value_type> ValuesList) {
    Values_.Assign(ValuesList);
    View_ = view_type(Values_.Data(), Extents);
    OVK_DEBUG_ASSERT((long long)(Values_.Count()) == View_.Count(), "Incorrect number of values in "
      "initializer list");
    return *this;
  }

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsInputIterator<IterType>())>
    array_base_1 &Assign(const interval_type &Extents, IterType First) {
    index_type NumValuesOld = View_.Count();
    index_type NumValues = Extents.Count();
    Values_.Reserve(NumValues);
    IterType Iter = First;
    for (auto &Value : Values_) {
      Value = *Iter++;
    }
    for (index_type i = 0; i < NumValues-NumValuesOld; ++i) {
      Values_.Append(*Iter++);
    }
    View_ = view_type(Values_.Data(), Extents);
    return *this;
  }

  array_base_1 &Assign(const array_base_1 &Other) {
    Values_ = Other.Values_;
    View_ = view_type(Values_.Data(), Other.View_.Extents());
    return *this;
  }

  array_base_1 &Assign(array_base_1 &&Other) noexcept {
    Values_ = std::move(Other.Values_);
    View_ = Other.View_;
    Other.View_ = view_type();
    return *this;
  }

  array_base_1 &Reserve(index_type NumValues) {
    Values_.Reserve(NumValues);
    View_ = view_type(Values_.Data(), View_.Extents());
    return *this;
  }

  array_base_1 &Resize(const interval_type &Extents) {
    using std::swap;
    indexer_type NewIndexer(Extents);
    index_type NumValues = Extents.Count();
    if (NewIndexer != View_.Indexer()) {
      index_type NumValuesOld = View_.Count();
      core::vector<value_type> OldValues(NumValuesOld);
      for (index_type i = 0; i < NumValuesOld; ++i) {
        swap(Values_[i], OldValues[i]);
      }
      view_type OldView(OldValues.Data(), View_.Extents());
      Values_.Resize(NumValues);
      View_ = view_type(Values_.Data(), Extents);
      core::ForEach<Layout>(OldView.Extents(), [&](const tuple_type &Tuple) {
        if (Extents.Contains(Tuple)) {
          View_(Tuple) = std::move(OldView(Tuple));
        }
      });
    } else {
      Values_.Resize(NumValues);
      View_ = view_type(Values_.Data(), Extents);
    }
    return *this;
  }

  array_base_1 &Resize(const interval_type &Extents, const value_type &Value) {
    using std::swap;
    indexer_type NewIndexer(Extents);
    index_type NumValues = Extents.Count();
    if (NewIndexer != View_.Indexer()) {
      index_type NumValuesOld = View_.Count();
      core::vector<value_type> OldValues(NumValuesOld, Value);
      for (index_type i = 0; i < NumValuesOld; ++i) {
        swap(Values_[i], OldValues[i]);
      }
      view_type OldView(OldValues.Data(), View_.Extents());
      Values_.Resize(NumValues, Value);
      View_ = view_type(Values_.Data(), Extents);
      core::ForEach<Layout>(OldView.Extents(), [&](const tuple_type &Tuple) {
        if (Extents.Contains(Tuple)) {
          View_(Tuple) = std::move(OldView(Tuple));
        }
      });
    } else {
      Values_.Resize(NumValues, Value);
      View_ = view_type(Values_.Data(), Extents);
    }
    return *this;
  }

  array_base_1 &Clear() {
    Values_.Clear();
    View_ = view_type();
    return *this;
  }

  OVK_FORCE_INLINE const value_type &operator()(TupleElementTypes... TupleElements) const {
    return View_(TupleElements...);
  }
  OVK_FORCE_INLINE value_type &operator()(TupleElementTypes... TupleElements) {
    return View_(TupleElements...);
  }

  OVK_FORCE_INLINE const value_type *Data(TupleElementTypes... TupleElements) const {
    return View_.Data(TupleElements...);
  }
  OVK_FORCE_INLINE value_type *Data(TupleElementTypes... TupleElements) {
    return View_.Data(TupleElements...);
  }

};

template <typename T, int Rank, array_layout Layout, typename=void> class array_base_2;

template <typename T, int Rank, array_layout Layout> class array_base_2<T, Rank, Layout,
 OVK_SPECIALIZATION_REQUIRES(Rank==1)> : public array_base_1<T, Rank, Layout,
  core::repeated_type_sequence_of_size<long long, Rank>> {

private:

  using parent_type = array_base_1<T, Rank, Layout, core::repeated_type_sequence_of_size<
    long long, Rank>>;

protected:

  using typename parent_type::value_type;
  using typename parent_type::index_type;
  using typename parent_type::view_type;
  using typename parent_type::iterator;
  using typename parent_type::const_iterator;
  using parent_type::Values_;
  using parent_type::View_;

public:

  using parent_type::parent_type;
  using parent_type::Assign;
  using parent_type::Reserve;
  using parent_type::Resize;
  using parent_type::Clear;
  using parent_type::operator();
  using parent_type::Data;

  value_type &Append(const value_type &Value) {
    Values_.Append(Value);
    View_ = view_type(Values_.Data(), {View_.Begin(0), View_.End(0)+1});
    return Values_.Back();
  }

  value_type &Append(value_type &&Value) {
    Values_.Append(std::move(Value));
    View_ = view_type(Values_.Data(), {View_.Begin(0), View_.End(0)+1});
    return Values_.Back();
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> value_type &Append(Args
    &&... Arguments) {
    Values_.Append(std::forward<Args>(Arguments)...);
    View_ = view_type(Values_.Data(), {View_.Begin(0), View_.End(0)+1});
    return Values_.Back();
  }

  iterator Insert(const_iterator Pos, const value_type &Value) {
    auto ValuesPos = Values_.Begin() + (Pos - View_.LinearBegin());
    auto ValuesIter = Values_.Insert(ValuesPos, Value);
    View_ = view_type(Values_.Data(), {View_.Begin(0), View_.End(0)+1});
    return View_.LinearBegin() + (ValuesIter - Values_.Begin());
  }

  iterator Insert(const_iterator Pos, value_type &&Value) {
    auto ValuesPos = Values_.Begin() + (Pos - View_.LinearBegin());
    auto ValuesIter = Values_.Insert(ValuesPos, std::move(Value));
    View_ = view_type(Values_.Data(), {View_.Begin(0), View_.End(0)+1});
    return View_.LinearBegin() + (ValuesIter - Values_.Begin());
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> iterator Insert(
    const_iterator Pos, Args &&... Arguments) {
    auto ValuesPos = Values_.Begin() + (Pos - View_.LinearBegin());
    auto ValuesIter = Values_.Insert(ValuesPos, std::forward<Args>(Arguments)...);
    View_ = view_type(Values_.Data(), {View_.Begin(0), View_.End(0)+1});
    return View_.LinearBegin() + (ValuesIter - Values_.Begin());
  }

  iterator Erase(const_iterator Pos) {
    auto ValuesPos = Values_.Begin() + (Pos - View_.LinearBegin());
    auto ValuesIter = Values_.Erase(ValuesPos);
    View_ = view_type(Values_.Data(), {View_.Begin(0), View_.End(0)-1});
    return View_.LinearBegin() + (ValuesIter - Values_.Begin());
  }

};

template <typename T, int Rank, array_layout Layout> class array_base_2<T, Rank, Layout,
 OVK_SPECIALIZATION_REQUIRES(Rank>1)> : public array_base_1<T, Rank, Layout,
  core::repeated_type_sequence_of_size<long long, Rank>> {

private:

  using parent_type = array_base_1<T, Rank, Layout, core::repeated_type_sequence_of_size<
    long long, Rank>>;

protected:

  using typename parent_type::value_type;
  using typename parent_type::index_type;
  using typename parent_type::view_type;
  using typename parent_type::iterator;
  using typename parent_type::const_iterator;
  using parent_type::Values_;
  using parent_type::View_;

public:

  using parent_type::parent_type;
  using parent_type::Assign;
  using parent_type::Reserve;
  using parent_type::Resize;
  using parent_type::Clear;
  using parent_type::operator();
  using parent_type::Data;

  // Don't know how to define this sensibly for ranks > 1
  template <typename... Args> value_type &Append(Args &&...) {
    // Assert on Args instead of just 'false' so it doesn't trigger unless method is instantiated
    static_assert(int(sizeof...(Args)) < 0, "Cannot use array::Append for ranks > 1.");
    return Values_.Back();
  }

  // Don't know how to define this sensibly for ranks > 1
  template <typename... Args> iterator Insert(Args &&...) {
    // Assert on Args instead of just 'false' so it doesn't trigger unless method is instantiated
    static_assert(int(sizeof...(Args)) < 0, "Cannot use array::Insert for ranks > 1.");
    return {};
  }

  // Don't know how to define this sensibly for ranks > 1
  template <typename... Args> iterator Erase(Args &&...) {
    // Assert on Args instead of just 'false' so it doesn't trigger unless method is instantiated
    static_assert(int(sizeof...(Args)) < 0, "Cannot use array::Erase for ranks > 1.");
    return {};
  }

};

}

template <typename T, int Rank_=1, array_layout Layout_=array_layout::ROW_MAJOR> class array :
  protected array_internal::array_base_2<T, Rank_, Layout_> {

private:

  using parent_type = array_internal::array_base_2<T, Rank_, Layout_>;

  using parent_type::Values_;
  using parent_type::View_;

public:

  using value_type = T;
  static constexpr const int Rank = Rank_;
  static constexpr const array_layout Layout = Layout_;
  using index_type = long long;
  using tuple_element_type = long long;
  using tuple_type = elem<tuple_element_type,Rank>;
  using interval_type = interval<index_type,Rank>;
  using indexer_type = indexer<index_type, tuple_element_type, Rank, Layout>;
  using view_type = array_view<value_type, Rank, Layout>;
  using const_view_type = array_view<const value_type, Rank, Layout>;
  using iterator = value_type *;
  using const_iterator = const value_type *;

  using parent_type::operator();
  using parent_type::Data;
  using parent_type::Append;
  using parent_type::Insert;
  using parent_type::Erase;

  array() = default;

  explicit array(const interval_type &Extents):
    parent_type(Extents)
  {}

  array(const interval_type &Extents, const value_type &Value):
    parent_type(Extents, Value)
  {}

  array(const interval_type &Extents, std::initializer_list<value_type> ValuesList):
    parent_type(Extents, ValuesList)
  {}

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsInputIterator<IterType>() &&
    std::is_convertible<core::iterator_deref_type<IterType>, value_type>::value)> array(const
    interval_type &Extents, IterType First):
    parent_type(Extents, First)
  {}

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U>::
    type, value_type>::value)> array(array_view<U, Rank, Layout> View):
    parent_type(View.Extents(), View.LinearBegin())
  {}

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U>::
    type, value_type>::value)> array(const interval_type &Extents, array_view<U, Rank, Layout>
    View):
    parent_type(Extents, View.LinearBegin())
  {}

  // Intel 17 didn't like using Rank and Layout instead of Rank_ and Layout_
  template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<
    ArrayRefType>>() && !std::is_same<core::remove_cvref<ArrayRefType>, array>::value &&
    core::ArrayHasFootprint<core::remove_cvref<ArrayRefType>, Rank_, Layout_>() &&
    std::is_convertible<core::array_access_type<ArrayRefType &&>, value_type>::value)>
    array(ArrayRefType &&Array):
    parent_type(core::ArrayExtents(Array), core::ArrayLinearBegin(std::forward<ArrayRefType>(
      Array)))
  {}

  // Intel 17 didn't like using Rank and Layout instead of Rank_ and Layout_
  template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<
    ArrayRefType>>() && !std::is_same<core::remove_cvref<ArrayRefType>, array>::value &&
    !core::IsIterator<typename std::decay<ArrayRefType>::type>() && core::ArrayHasFootprint<
    core::remove_cvref<ArrayRefType>, Rank_, Layout_>() && std::is_convertible<
    core::array_access_type<ArrayRefType &&>, value_type>::value)> array(const interval_type
    &Extents, ArrayRefType &&Array):
    parent_type(Extents, core::ArrayLinearBegin(std::forward<ArrayRefType>(Array)))
  {}

  array(const array &Other) = default;

  array(const interval_type &Extents, const array &Other):
    parent_type(Extents, Other)
  {}

  array(array &&Other) noexcept = default;

  array(const interval_type &Extents, array &&Other) noexcept:
    parent_type(Extents, std::move(Other))
  {}

  array &operator=(const array &Other) = default;

  array &operator=(array &&Other) noexcept = default;

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U>::
    type, value_type>::value)> array &operator=(array_view<U, Rank, Layout> View) {
    parent_type::Assign(View.Extents(), View.LinearBegin());
    return *this;
  }

  // Intel 17 didn't like using Rank and Layout instead of Rank_ and Layout_
  template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<
    ArrayRefType>>() && !std::is_same<core::remove_cvref<ArrayRefType>, array>::value &&
    core::ArrayHasFootprint<core::remove_cvref<ArrayRefType>, Rank_, Layout_>() &&
    std::is_convertible<core::array_access_type<ArrayRefType &&>, value_type>::value)> array
    &operator=(ArrayRefType &&Array) {
    parent_type::Assign(core::ArrayExtents(Array), core::ArrayLinearBegin(std::forward<
      ArrayRefType>(Array)));
    return *this;
  }

  array &Assign(const interval_type &Extents, const value_type &Value) {
    parent_type::Assign(Extents, Value);
    return *this;
  }

  array &Assign(const interval_type &Extents, std::initializer_list<value_type> ValuesList) {
    parent_type::Assign(Extents, ValuesList);
    return *this;
  }

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsInputIterator<IterType>() &&
    std::is_convertible<core::iterator_deref_type<IterType>, value_type>::value)> array
    &Assign(const interval_type &Extents, IterType First) {
    parent_type::Assign(Extents, First);
    return *this;
  }

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U>::
    type, value_type>::value)> array &Assign(array_view<U, Rank, Layout> View) {
    parent_type::Assign(View.Extents(), View.LinearBegin());
    return *this;
  }

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U>::
    type, value_type>::value)> array &Assign(const interval_type &Extents, array_view<U, Rank,
    Layout> View) {
    parent_type::Assign(Extents, View.LinearBegin());
    return *this;
  }

  // Intel 17 didn't like using Rank and Layout instead of Rank_ and Layout_
  template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<
    ArrayRefType>>() && !std::is_same<core::remove_cvref<ArrayRefType>, array>::value &&
    core::ArrayHasFootprint<core::remove_cvref<ArrayRefType>, Rank_, Layout_>() &&
    std::is_convertible<core::array_access_type<ArrayRefType &&>, value_type>::value)>
    array &Assign(ArrayRefType &&Array) {
    parent_type::Assign(core::ArrayExtents(Array), core::ArrayLinearBegin(std::forward<
      ArrayRefType>(Array)));
    return *this;
  }

  // Intel 17 didn't like using Rank and Layout instead of Rank_ and Layout_
  template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<
    ArrayRefType>>() && !std::is_same<core::remove_cvref<ArrayRefType>, array>::value &&
    !core::IsIterator<typename std::decay<ArrayRefType>::type>() && core::ArrayHasFootprint<
    core::remove_cvref<ArrayRefType>, Rank_, Layout_>() && std::is_convertible<
    core::array_access_type<ArrayRefType &&>, value_type>::value)> array &Assign(const
    interval_type &Extents, ArrayRefType &&Array) {
    parent_type::Assign(Extents, core::ArrayLinearBegin(std::forward<ArrayRefType>(Array)));
    return *this;
  }

  array &Assign(const array &Other) {
    parent_type::Assign(Other);
    return *this;
  }

  array &Assign(const interval_type &Extents, const array &Other) {
    parent_type::Assign(Extents, Other);
    return *this;
  }

  array &Assign(array &&Other) noexcept {
    parent_type::Assign(std::move(Other));
    return *this;
  }

  array &Assign(const interval_type &Extents, array &&Other) noexcept {
    parent_type::Assign(Extents, std::move(Other));
    return *this;
  }

  array &Reserve(index_type NumValues) {
    parent_type::Reserve(NumValues);
    return *this;
  }

  array &Resize(const interval_type &Extents) {
    parent_type::Resize(Extents);
    return *this;
  }

  array &Resize(const interval_type &Extents, const value_type &Value) {
    parent_type::Resize(Extents, Value);
    return *this;
  }

  array &Clear() {
    parent_type::Clear();
    return *this;
  }

  array &Fill(const value_type &Value) {
    View_.Fill(Value);
    return *this;
  }

  array &Fill(std::initializer_list<value_type> Values) {
    View_.Fill(Values);
    return *this;
  }

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsInputIterator<IterType>() &&
    std::is_convertible<core::iterator_deref_type<IterType>, value_type>::value)> array
    &Fill(IterType First) {
    View_.Fill(First);
    return *this;
  }

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U>::
    type, value_type>::value)> array &Fill(array_view<U, Rank, Layout> View) {
    View_.Fill(View);
    return *this;
  }

  // Intel 17 didn't like using Rank and Layout instead of Rank_ and Layout_
  template <typename ArrayRefType, OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<
    ArrayRefType>>() && !core::IsIterator<typename std::decay<ArrayRefType>::type>() &&
    core::ArrayHasFootprint<core::remove_cvref<ArrayRefType>, Rank_, Layout_>() &&
    std::is_convertible<core::array_access_type<ArrayRefType &&>, value_type>::value)> array
    &Fill(ArrayRefType &&Array) {
    View_.Fill(std::forward<ArrayRefType>(Array));
    return *this;
  }

  OVK_FORCE_INLINE const interval_type &Extents() const { return View_.Extents(); }

  OVK_FORCE_INLINE const tuple_type &Begin() const { return View_.Begin(); }
  OVK_FORCE_INLINE tuple_element_type Begin(int iDim) const { return View_.Begin(iDim); }

  OVK_FORCE_INLINE const tuple_type &End() const { return View_.End(); }
  OVK_FORCE_INLINE tuple_element_type End(int iDim) const { return View_.End(iDim); }

  tuple_type Size() const { return View_.Size(); }
  tuple_element_type Size(int iDim) const { return View_.Size(iDim); }

  index_type Count() const { return View_.Count(); }

  bool Empty() const { return View_.Empty(); }

  index_type Capacity() const { return index_type(Values_.Capacity()); }

  const indexer_type &Indexer() { return View_.Indexer(); }

  // Want to use iterator directly here instead of constructing intermediate elem type
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsRandomAccessIterator<IterType>() &&
    std::is_convertible<core::iterator_deref_type<IterType>, tuple_element_type>::value)>
    OVK_FORCE_INLINE const value_type &operator()(IterType First) const {
    return View_(First);
  }
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsRandomAccessIterator<IterType>() &&
    std::is_convertible<core::iterator_deref_type<IterType>, tuple_element_type>::value)>
    OVK_FORCE_INLINE value_type &operator()(IterType First) {
    return View_(First);
  }

  // Want to use array directly here instead of constructing intermediate elem type
  // Intel 17 didn't like using Rank instead of Rank_
  template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
    !core::IsIterator<typename std::decay<ArrayType>::type>() && std::is_convertible<
    core::array_access_type<const ArrayType &>, tuple_element_type>::value && core::ArrayRank<
    ArrayType>() == 1 && (core::ArrayHasRuntimeExtents<ArrayType>() || (core::StaticArrayHasBegin<
    ArrayType,0>() && core::StaticArrayHasEnd<ArrayType,Rank_>())))> OVK_FORCE_INLINE const
    value_type &operator()(const ArrayType &Array) const {
    return View_(Array);
  }
  template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
    !core::IsIterator<typename std::decay<ArrayType>::type>() && std::is_convertible<
    core::array_access_type<const ArrayType &>, tuple_element_type>::value && core::ArrayRank<
    ArrayType>() == 1 && (core::ArrayHasRuntimeExtents<ArrayType>() || (core::StaticArrayHasBegin<
    ArrayType,0>() && core::StaticArrayHasEnd<ArrayType,Rank_>())))> OVK_FORCE_INLINE value_type
    &operator()(const ArrayType &Array) {
    return View_(Array);
  }

  const value_type &operator[](int Index) const { return View_[Index]; }
  value_type &operator[](int Index) { return View_[Index]; }

  OVK_FORCE_INLINE const value_type *Data() const { return View_.Data(); }
  OVK_FORCE_INLINE value_type *Data() { return View_.Data(); }
  // TODO: Change this to be like operator() above
  OVK_FORCE_INLINE const value_type *Data(const tuple_type &Tuple) const {
    return View_.Data(Tuple);
  }
  OVK_FORCE_INLINE value_type *Data(const tuple_type &Tuple) {
    return View_.Data(Tuple);
  }

  const_iterator LinearBegin() const { return View_.LinearBegin(); }
  iterator LinearBegin() { return View_.LinearBegin(); }
  const_iterator LinearEnd() const { return View_.LinearEnd(); }
  iterator LinearEnd() { return View_.LinearEnd(); }

  // Google Test doesn't use free begin/end functions and instead expects container to have
  // lowercase begin/end methods
  const_iterator begin() const { return LinearBegin(); }
  iterator begin() { return LinearBegin(); }
  const_iterator end() const { return LinearEnd(); }
  iterator end() { return LinearEnd(); }

private:

  friend class core::test_helper<array>;

};

template <typename T, int Rank_, array_layout Layout_> struct array_traits<array<T, Rank_, Layout_>>
  {
  using value_type = T;
  static constexpr const int Rank = Rank_;
  static constexpr const array_layout Layout = Layout_;
  template <int iDim> static long long Begin(const array<T, Rank, Layout> &Array) {
    return Array.Begin(iDim);
  }
  template <int iDim> static long long End(const array<T, Rank, Layout> &Array) {
    return Array.End(iDim);
  }
  static const T *Data(const array<T, Rank, Layout> &Array) { return Array.Data(); }
  static T *Data(array<T, Rank, Layout> &Array) { return Array.Data(); }
};

template <typename T, int Rank, array_layout Layout> typename array<T, Rank, Layout>::iterator
  begin(array<T, Rank, Layout> &Array) {
  return Array.LinearBegin();
}

template <typename T, int Rank, array_layout Layout> typename array<T, Rank, Layout>::const_iterator
  begin(const array<T, Rank, Layout> &Array) {
  return Array.LinearBegin();
}

template <typename T, int Rank, array_layout Layout> typename array<T, Rank, Layout>::iterator
  end(array<T, Rank, Layout> &Array) {
  return Array.LinearEnd();
}

template <typename T, int Rank, array_layout Layout> typename array<T, Rank, Layout>::const_iterator
  end(const array<T, Rank, Layout> &Array) {
  return Array.LinearEnd();
}

}

#endif
