// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_VECTOR_HPP_INCLUDED
#define OVK_CORE_VECTOR_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/PointerIterator.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <initializer_list>
#include <iterator>
#include <memory>
#include <utility>
#include <type_traits>
#include <vector>

namespace ovk {
namespace core {

namespace vector_internal {

enum class not_bool : bool { FALSE=false, TRUE=true };
template <typename T> struct no_bool_helper { using type = T; };
template <> struct no_bool_helper<bool> { using type = not_bool; };
template <typename T> using no_bool = typename no_bool_helper<T>::type;

template <typename Allocator, typename T, OVK_FUNCTION_REQUIRES(!std::is_same<T, bool>::value)>
  std::vector<T, Allocator> ConstructFromInitializerList(std::initializer_list<T> ValuesList) {
  return {ValuesList};
}
template <typename Allocator> std::vector<not_bool, Allocator> ConstructFromInitializerList(
  std::initializer_list<bool> ValuesList) {
  std::vector<not_bool, Allocator> Values;
  Values.reserve(ValuesList.size());
  auto Iter = ValuesList.begin();
  while (Iter != ValuesList.end()) {
    Values.push_back(static_cast<not_bool>(*Iter++));
  }
  return Values;
}

template <typename Allocator, typename T, OVK_FUNCTION_REQUIRES(!std::is_same<T, bool>::value)>
  void AssignInitializerList(std::vector<T, Allocator> &Values, std::initializer_list<T> ValuesList) {
  return Values.assign(ValuesList);
}
template <typename Allocator> void AssignInitializerList(std::vector<not_bool, Allocator> &Values,
  std::initializer_list<bool> ValuesList) {
  Values.reserve(ValuesList.size());
  auto Iter = ValuesList.begin();
  for (auto &Value : Values) {
    Value = static_cast<not_bool>(*Iter++);
  }
  while (Iter != ValuesList.end()) {
    Values.push_back(static_cast<not_bool>(*Iter++));
  }
}

template <typename Allocator, typename T, typename IterType, OVK_FUNCTION_REQUIRES(!std::is_same<T,
  bool>::value)> void AssignIterators(std::vector<T, Allocator> &Values, IterType Begin, IterType
  End) {
  return Values.assign(Begin, End);
}
template <typename Allocator, typename IterType, OVK_FUNCTION_REQUIRES(IsForwardIterator<IterType
  >())> void AssignIterators(std::vector<not_bool, Allocator> &Values, IterType Begin, IterType End) {
  Values.reserve(std::distance(Begin, End));
  auto Iter = Begin;
  for (auto &Value : Values) {
    Value = static_cast<not_bool>(*Iter++);
  }
  while (Iter != End) {
    Values.push_back(static_cast<not_bool>(*Iter++));
  }
}
template <typename Allocator, typename IterType, OVK_FUNCTION_REQUIRES(!IsForwardIterator<IterType
  >() && IsInputIterator<IterType>())> void AssignIterators(std::vector<not_bool, Allocator>
  &Values, IterType Begin, IterType End) {
  auto Iter = Begin;
  for (auto &Value : Values) {
    Value = static_cast<not_bool>(*Iter++);
  }
  while (Iter != End) {
    Values.push_back(static_cast<not_bool>(*Iter++));
  }
}

}

// Wrapper around std::vector to avoid the abomination that is std::vector<bool>
template <typename T, typename Allocator=std::allocator<T>> class vector {

public:

  using value_type = T;
  using allocator_type = Allocator;

private:

  using storage_value_type = vector_internal::no_bool<value_type>;
  using storage_allocator_type = typename std::allocator_traits<Allocator>::template rebind_alloc<
    storage_value_type>;
  using storage_type = std::vector<storage_value_type, storage_allocator_type>;
  using storage_iterator = typename storage_type::iterator;
  using const_storage_iterator = typename storage_type::const_iterator;

public:

  using index_type = long long;
  using iterator = core::pointer_iterator<vector, value_type *>;
  using const_iterator = core::pointer_iterator<vector, const value_type *>;

  vector() = default;

  explicit vector(index_type NumValues):
    Values_(NumValues)
  {}

  vector(index_type NumValues, const value_type &Value):
    Values_(NumValues, reinterpret_cast<const storage_value_type &>(Value))
  {}

  vector(std::initializer_list<value_type> ValuesList):
    Values_(vector_internal::ConstructFromInitializerList<storage_allocator_type>(ValuesList))
  {}

  vector(const vector &Other) = default;
  // Can't use "= default" because std::vector is not noexcept movable until C++17
  vector(vector &&Other) noexcept:
    Values_(std::move(Other.Values_))
  {}

  vector &operator=(const vector &Other) = default;
  // Can't use "= default" because std::vector is not noexcept movable until C++17
  vector &operator=(vector &&Other) noexcept {
    Values_ = std::move(Other.Values_);
    return *this;
  }

  vector &operator=(std::initializer_list<value_type> ValuesList) {
    vector_internal::AssignInitializerList(Values_, ValuesList);
    return *this;
  }

  vector &Assign(index_type NumValues, const value_type &Value) {
    Values_.assign(NumValues, reinterpret_cast<const storage_value_type &>(Value));
    return *this;
  }

  vector &Assign(std::initializer_list<value_type> ValuesList) {
    vector_internal::AssignInitializerList(Values_, ValuesList);
    return *this;
  }

  template <typename IterType, OVK_FUNCTION_REQUIRES(IsInputIterator<IterType>() &&
    std::is_same<iterator_value_type<IterType>, value_type>::value)> vector &Assign(IterType Begin,
    IterType End) {
    vector_internal::AssignIterators(Values_, Begin, End);
    return *this;
  }

  void Reserve(index_type NumValues) {
    Values_.reserve(NumValues);
  }

  void Resize(index_type NumValues) {
    Values_.resize(NumValues);
  }

  void Resize(index_type NumValues, const value_type &Value) {
    Values_.resize(NumValues, reinterpret_cast<const storage_value_type &>(Value));
  }

  void Clear() {
    Values_.clear();
  }

  value_type &Append(const value_type &Value) {
    Values_.emplace_back(reinterpret_cast<const storage_value_type &>(Value));
    return reinterpret_cast<value_type &>(Values_.back());
  }

  value_type &Append(value_type &&Value) {
    Values_.emplace_back(std::move(reinterpret_cast<storage_value_type &>(Value)));
    return reinterpret_cast<value_type &>(Values_.back());
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> value_type &Append(Args
    &&... Arguments) {
    value_type Value(std::forward<Args>(Arguments)...);
    Values_.emplace_back(std::move(reinterpret_cast<storage_value_type &>(Value)));
    return reinterpret_cast<value_type &>(Values_.back());
  }

  value_type &Insert(index_type iValue, const value_type &Value) {
    Values_.insert(Values_.begin()+iValue, reinterpret_cast<const storage_value_type &>(Value));
    return *(Begin() + iValue);
  }

  value_type &Insert(index_type iValue, value_type &&Value) {
    Values_.insert(Values_.begin()+iValue, std::move(reinterpret_cast<storage_value_type
      &>(Value)));
    return *(Begin() + iValue);
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> value_type &Insert(
    index_type iValue, Args &&... Arguments) {
    value_type Value(std::forward<Args>(Arguments)...);
    Values_.insert(Values_.begin()+iValue, std::move(reinterpret_cast<storage_value_type
      &>(Value)));
    return *(Begin() + iValue);
  }

  iterator Insert(const_iterator Pos, const value_type &Value) {
    index_type iValue = index_type(Pos - Begin());
    Values_.insert(Values_.begin()+iValue, reinterpret_cast<const storage_value_type &>(Value));
    return Begin() + iValue;
  }

  iterator Insert(const_iterator Pos, value_type &&Value) {
    index_type iValue = index_type(Pos - Begin());
    Values_.insert(Values_.begin()+iValue, std::move(reinterpret_cast<storage_value_type
      &>(Value)));
    return Begin() + iValue;
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> iterator Insert(
    const_iterator Pos, Args &&... Arguments) {
    value_type Value(std::forward<Args>(Arguments)...);
    index_type iValue = index_type(Pos - Begin());
    Values_.insert(Values_.begin()+iValue, std::move(reinterpret_cast<storage_value_type &>(
      Value)));
    return Begin() + iValue;
  }

  void Erase(index_type iValue) {
    Values_.erase(Values_.begin()+iValue);
  }

  iterator Erase(const_iterator Pos) {
    index_type iValue = index_type(Pos - Begin());
    Values_.erase(Values_.begin()+iValue);
    return Begin() + iValue;
  }

  index_type Count() const { return index_type(Values_.size()); }
  index_type Capacity() const { return index_type(Values_.capacity()); }

  bool Empty() const { return Values_.empty(); }

  OVK_FORCE_INLINE const value_type &operator[](index_type iValue) const {
    return reinterpret_cast<const value_type &>(Values_[iValue]);
  }
  OVK_FORCE_INLINE value_type &operator[](index_type iValue) {
    return reinterpret_cast<value_type &>(Values_[iValue]);
  }

  const value_type &Front() const {
    return reinterpret_cast<const value_type &>(Values_.front());
  }
  value_type &Front() {
    return reinterpret_cast<value_type &>(Values_.front());
  }

  const value_type &Back() const {
    return reinterpret_cast<const value_type &>(Values_.back());
  }
  value_type &Back() {
    return reinterpret_cast<value_type &>(Values_.back());
  }

  OVK_FORCE_INLINE const value_type *Data() const {
    return reinterpret_cast<const value_type *>(Values_.data());
  }
  OVK_FORCE_INLINE value_type *Data() {
    return reinterpret_cast<value_type *>(Values_.data());
  }

  OVK_FORCE_INLINE const_iterator Begin() const {
    return const_iterator(reinterpret_cast<const value_type *>(Values_.data()));
  }
  OVK_FORCE_INLINE iterator Begin() {
    return iterator(reinterpret_cast<value_type *>(Values_.data()));
  }

  OVK_FORCE_INLINE const_iterator End() const {
    return const_iterator(reinterpret_cast<const value_type *>(Values_.data() + Values_.size()));
  }
  OVK_FORCE_INLINE iterator End() {
    return iterator(reinterpret_cast<value_type *>(Values_.data() + Values_.size()));
  }

  // Google Test doesn't use free begin/end functions and instead expects container to have
  // lowercase begin/end methods
  OVK_FORCE_INLINE const_iterator begin() const { return Begin(); }
  OVK_FORCE_INLINE iterator begin() { return Begin(); }
  OVK_FORCE_INLINE const_iterator end() const { return End(); }
  OVK_FORCE_INLINE iterator end() { return End(); }

private:

  storage_type Values_;

  friend class test_helper<vector>;

};

template <typename T, typename Allocator> OVK_FORCE_INLINE typename vector<T, Allocator>::iterator
  begin(vector<T, Allocator> &Vector) {
  return Vector.Begin();
}

template <typename T, typename Allocator> OVK_FORCE_INLINE typename vector<T, Allocator>::
  const_iterator begin(const vector<T, Allocator> &Vector) {
  return Vector.Begin();
}

template <typename T, typename Allocator> OVK_FORCE_INLINE typename vector<T, Allocator>::iterator
  end(vector<T, Allocator> &Vector) {
  return Vector.End();
}

template <typename T, typename Allocator> OVK_FORCE_INLINE typename vector<T, Allocator>::
  const_iterator end(const vector<T, Allocator> &Vector) {
  return Vector.End();
}

}

template <typename T, typename Allocator> struct array_traits<core::vector<T, Allocator>> {
  using value_type = T;
  static constexpr int Rank = 1;
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static long long ExtentBegin(const core::vector<T, Allocator> &) { return 0; }
  template <int> static long long ExtentEnd(const core::vector<T, Allocator> &Vec) {
    return Vec.Count();
  }
  static const T *Data(const core::vector<T, Allocator> &Vec) { return Vec.Data(); }
  static T *Data(core::vector<T, Allocator> &Vec) { return Vec.Data(); }
};

}

#endif
