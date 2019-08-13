// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_SET_HPP_INCLUDED
#define OVK_CORE_SET_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/PointerIterator.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <algorithm>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <utility>

namespace ovk {

template <typename ValueType, typename CompareType=std::less<ValueType>> class set {

public:

  using value_type = ValueType;
  using compare_type = CompareType;
  using index_type = long long;
  // Values aren't mutable, so iterator is const
  using iterator = core::pointer_iterator<set, const value_type *>;
  using const_iterator = iterator;

  set() = default;

  explicit set(compare_type Compare):
    Compare_(std::move(Compare))
  {}

  set(std::initializer_list<value_type> ValuesList) {
    Values_.Reserve(ValuesList.size());
    for (auto &Value : ValuesList) {
      Insert(Value);
    }
  }

  // Don't care about input iterators enough to implement a second overload
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>() &&
    std::is_convertible<core::iterator_reference_type<IterType>, value_type>::value)>
    set(IterType Begin, IterType End) {
    Values_.Reserve(std::distance(Begin, End));
    IterType Iter = Begin;
    while (Iter != End) {
      Insert(*Iter);
      ++Iter;
    }
  }

  // Don't care about input iterators enough to implement a second overload
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>() &&
    std::is_convertible<core::iterator_reference_type<IterType>, value_type>::value)>
    set(IterType Begin, IterType End, compare_type Compare):
    Compare_(std::move(Compare))
  {
    Values_.Reserve(std::distance(Begin, End));
    IterType Iter = Begin;
    while (Iter != End) {
      Insert(*Iter);
      ++Iter;
    }
  }

  set &operator=(std::initializer_list<value_type> ValuesList) {
    return Assign(ValuesList);
  }

  set &Assign(std::initializer_list<value_type> ValuesList) {
    Values_.Clear();
    Values_.Reserve(ValuesList.size());
    for (auto &Value : ValuesList) {
      Insert(Value);
    }
    return *this;
  }

  // Don't care about input iterators enough to implement a second overload
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>()
    && std::is_convertible<core::iterator_reference_type<IterType>, value_type>::value)>
    set &Assign(IterType Begin, IterType End) {
    Values_.Clear();
    Values_.Reserve(std::distance(Begin, End));
    IterType Iter = Begin;
    while (Iter != End) {
      Insert(*Iter);
      ++Iter;
    }
    return *this;
  }

  void Reserve(index_type Count) {
    Values_.Reserve(Count);
  }

  void Insert(const value_type &Value) {
    auto ValuesIter = LowerBound_(Value);
    if (ValuesIter == Values_.End() || Compare_(Value, *ValuesIter)) {
      ValuesIter = Values_.Insert(ValuesIter, Value);
    }
  }

  iterator Insert(iterator LowerBoundIter, const value_type &Value) {
    auto ValuesIter = Values_.Begin() + (LowerBoundIter - Begin());
    if (ValuesIter == Values_.End() || Compare_(Value, *ValuesIter)) {
      ValuesIter = Values_.Insert(ValuesIter, Value);
    }
    return iterator(Values_.Data() + (ValuesIter - Values_.Begin()));
  }

  void Erase(const value_type &Value) {
    auto ValuesIter = LowerBound_(Value);
    if (ValuesIter != Values_.End() && !Compare_(Value, *ValuesIter)) {
      Values_.Erase(ValuesIter);
    }
  }

  iterator Erase(iterator Pos) {
    auto ValuesIter = Values_.Begin() + (Pos - Begin());
    ValuesIter = Values_.Erase(ValuesIter);
    return iterator(Values_.Data() + (ValuesIter - Values_.Begin()));
  }

  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableAs<F &&, bool(const value_type &)>())>
    void EraseIf(F &&Predicate) {
    auto ValuesIter = Values_.Begin();
    while (ValuesIter != Values_.End()) {
      if (std::forward<F>(Predicate)(*ValuesIter)) {
        ValuesIter = Values_.Erase(ValuesIter);
      } else {
        ++ValuesIter;
      }
    }
  }

  void Clear() {
    Values_.Clear();
  }

  bool Contains(const value_type &Value) const {
    auto ValuesIter = LowerBound_(Value);
    return ValuesIter != Values_.End() && !Compare_(Value, *ValuesIter);
  }

  iterator Find(const value_type &Value) const {
    auto ValuesIter = LowerBound_(Value);
    if (ValuesIter != Values_.End() && Compare_(Value, *ValuesIter)) {
      ValuesIter = Values_.End();
    }
    return iterator(Values_.Data() + (ValuesIter - Values_.Begin()));
  }

  iterator LowerBound(const value_type &Value) const {
    auto ValuesIter = LowerBound_(Value);
    return iterator(Values_.Data() + (ValuesIter - Values_.Begin()));
  }

  iterator UpperBound(const value_type &Value) const {
    auto ValuesIter = UpperBound_(Value);
    return iterator(Values_.Data() + (ValuesIter - Values_.Begin()));
  }

  index_type Count() const { return Values_.Count(); }

  bool Empty() const { return Values_.Empty(); }

  index_type Capacity() const { return Values_.Capacity(); }

  const value_type &operator[](index_type Index) const { return Values_(Index); }

  const value_type *Data() const { return Values_.Data(); }

  iterator Begin() const { return iterator(Values_.Data()); }

  iterator End() const { return iterator(Values_.Data() + Values_.Count()); }

  const compare_type &Compare() const { return Compare_; }

  bool Compare(const value_type &Left, const value_type &Right) const {
    return Compare_(Left, Right);
  }

private:

  array<value_type> Values_;
  compare_type Compare_;

  typename array<value_type>::const_iterator LowerBound_(const value_type &Value) const {
    return std::lower_bound(Values_.Begin(), Values_.End(), Value, Compare_);
  }

  typename array<value_type>::iterator LowerBound_(const value_type &Value) {
    return std::lower_bound(Values_.Begin(), Values_.End(), Value, Compare_);
  }

  typename array<value_type>::const_iterator UpperBound_(const value_type &Value) const {
    return std::upper_bound(Values_.Begin(), Values_.End(), Value, Compare_);
  }

  typename array<value_type>::iterator UpperBound_(const value_type &Value) {
    return std::upper_bound(Values_.Begin(), Values_.End(), Value, Compare_);
  }

  friend class core::test_helper<set>;

};

template <typename ValueType, typename CompareType> typename set<ValueType, CompareType>::iterator
  begin(const set<ValueType, CompareType> &Set) {
  return Set.Begin();
}

template <typename ValueType, typename CompareType> typename set<ValueType, CompareType>::iterator
  end(const set<ValueType, CompareType> &Set) {
  return Set.End();
}

template <typename ValueType, typename CompareType> struct array_traits<set<ValueType, CompareType>>
  {
  using set_type = set<ValueType, CompareType>;
  using value_type = typename set_type::value_type;
  static constexpr int Rank = 1;
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static long long ExtentBegin(const set_type &) { return 0; }
  template <int> static long long ExtentEnd(const set_type &Set) { return Set.Count(); }
  static const value_type *Data(const set_type &Set) { return Set.Data(); }
  // No non-const Data access
};

}

#endif
