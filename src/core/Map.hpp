// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MAP_HPP_INCLUDED
#define OVK_CORE_MAP_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/PointerIterator.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/Set.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {

namespace map_internal {

template <typename T, bool Contiguous> class value_wrapper;

template <typename T> class value_wrapper<T, true> {

public:

  value_wrapper(const T &Value):
    Value_(Value)
  {}

  value_wrapper(T &&Value):
    Value_(std::move(Value))
  {}

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<T, Args &&...>::value &&
    !core::IsCopyOrMoveArgument<T, Args &&...>())> value_wrapper(Args &&... Arguments):
    Value_(std::forward<Args>(Arguments)...)
  {}

  value_wrapper(const value_wrapper &Other) = default;
  value_wrapper(value_wrapper &&Other) noexcept = default;
  value_wrapper &operator=(const value_wrapper &Other) = default;
  value_wrapper &operator=(value_wrapper &&Other) noexcept = default;

  const T &Get() const { return Value_; }
  T &Get() { return Value_; }

private:

  T Value_;

};

template <typename T> class value_wrapper<T, false> {

public:

  value_wrapper(const T &Value):
    Value_(new T(Value))
  {}

  value_wrapper(T &&Value):
    Value_(new T(std::move(Value)))
  {}

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<T, Args &&...>::value &&
    !core::IsCopyOrMoveArgument<T, Args &&...>())> value_wrapper(Args &&... Arguments):
    Value_(new T(std::forward<Args>(Arguments)...))
  {}

  value_wrapper(const value_wrapper &Other):
    Value_(new T(*Other.Value_))
  {}

  value_wrapper(value_wrapper &&Other) noexcept:
    Value_(std::move(Other.Value_))
  {}

  value_wrapper &operator=(value_wrapper Other) noexcept {
    using std::swap;
    swap(Value_, Other.Value_);
    return *this;
  }

  const T &Get() const { return *Value_; }
  T &Get() { return *Value_; }

private:

  std::unique_ptr<T> Value_;

};

}

template <typename T> constexpr bool MapContiguousDefault() {
  // 16 seems like a good number *shrug*
  return sizeof(T) <= 16;
}

template <typename KeyType, typename ValueType, bool Contiguous_=MapContiguousDefault<ValueType>()>
  class map_entry {

public:

  using key_type = KeyType;
  using value_type = ValueType;
  static constexpr bool Contiguous = Contiguous_;

  map_entry(const key_type &Key, const value_type &Value):
    Key_(Key),
    Value_(Value)
  {}

  map_entry(const key_type &Key, value_type &&Value):
    Key_(Key),
    Value_(std::move(Value))
  {}

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> map_entry(const key_type
    &Key, Args &&... Arguments):
    Key_(Key),
    Value_(std::forward<Args>(Arguments)...)
  {}

  const key_type &Key() const { return Key_; }
  int Key(int iDim) const { return Key_(iDim); }

  const value_type &Value() const { return Value_.Get(); }
  value_type &Value() { return Value_.Get(); }

private:

  using value_wrapper_type = map_internal::value_wrapper<value_type, Contiguous>;

  key_type Key_;
  value_wrapper_type Value_;

  friend class core::test_helper<map_entry>;

};

template <typename KeyType, typename ValueType, typename KeyCompareType=std::less<KeyType>,
  bool Contiguous_=MapContiguousDefault<ValueType>()> class map {

public:

  using key_type = KeyType;
  using value_type = ValueType;
  using key_compare_type = KeyCompareType;
  using key_set_type = set<KeyType, KeyCompareType>;
  static constexpr bool Contiguous = Contiguous_;
  using index_type = long long;
  using entry = map_entry<KeyType,ValueType,Contiguous>;
  using iterator = core::pointer_iterator<map, entry *>;
  using const_iterator = core::pointer_iterator<map, const entry *>;

  map() = default;

  map(key_compare_type KeyCompare):
    Keys_(std::move(KeyCompare))
  {}

  map(std::initializer_list<entry> EntriesList) {
    Keys_.Reserve(EntriesList.size());
    Entries_.Reserve(EntriesList.size());
    for (auto &Entry : EntriesList) {
      Insert(Entry);
    }
  }

  map(std::initializer_list<entry> EntriesList, key_compare_type KeyCompare):
    Keys_(std::move(KeyCompare))
  {
    Keys_.Reserve(EntriesList.size());
    Entries_.Reserve(EntriesList.size());
    for (auto &Entry : EntriesList) {
      Insert(Entry);
    }
  }

  // Don't care about input iterators enough to implement a second overload
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>() &&
    std::is_convertible<core::iterator_reference_type<IterType>, entry>::value)>
    map(IterType Begin, IterType End) {
    index_type NumEntries = index_type(std::distance(Begin, End));
    Keys_.Reserve(NumEntries);
    Entries_.Reserve(NumEntries);
    IterType Iter = Begin;
    while (Iter != End) {
      Insert(*Iter);
      ++Iter;
    }
  }

  // Don't care about input iterators enough to implement a second overload
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>() &&
    std::is_convertible<core::iterator_reference_type<IterType>, entry>::value)>
    map(IterType Begin, IterType End, key_compare_type KeyCompare):
    Keys_(std::move(KeyCompare))
  {
    index_type NumEntries = index_type(std::distance(Begin, End));
    Keys_.Reserve(NumEntries);
    Entries_.Reserve(NumEntries);
    IterType Iter = Begin;
    while (Iter != End) {
      Insert(*Iter);
      ++Iter;
    }
  }

  map &operator=(std::initializer_list<entry> EntriesList) {
    return Assign(EntriesList);
  }

  map &Assign(std::initializer_list<entry> EntriesList) {
    Keys_.Clear();
    Entries_.Clear();
    Keys_.Reserve(EntriesList.size());
    Entries_.Reserve(EntriesList.size());
    for (auto &Entry : EntriesList) {
      Insert(Entry);
    }
    return *this;
  }

  // Don't care about input iterators enough to implement a second overload
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>() &&
    std::is_convertible<core::iterator_reference_type<IterType>, entry>::value)>
    map &Assign(IterType Begin, IterType End) {
    index_type NumEntries = index_type(std::distance(Begin, End));
    Keys_.Clear();
    Entries_.Clear();
    Keys_.Reserve(NumEntries);
    Entries_.Reserve(NumEntries);
    IterType Iter = Begin;
    while (Iter != End) {
      Insert(*Iter);
      ++Iter;
    }
    return *this;
  }

  value_type &Insert(const entry &Entry) {
    auto KeysIter = Keys_.LowerBound(Entry.Key());
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Compare(Entry.Key(), *KeysIter)) {
      Keys_.Insert(KeysIter, Entry.Key());
      EntriesIter = Entries_.Insert(EntriesIter, Entry);
    } else {
      *EntriesIter = Entry;
    }
    return EntriesIter->Value();
  }

  value_type &Insert(entry &&Entry) {
    auto KeysIter = Keys_.LowerBound(Entry.Key());
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Compare(Entry.Key(), *KeysIter)) {
      Keys_.Insert(KeysIter, Entry.Key());
      EntriesIter = Entries_.Insert(EntriesIter, std::move(Entry));
    } else {
      *EntriesIter = std::move(Entry);
    }
    return EntriesIter->Value();
  }

  value_type &Insert(const key_type &Key, const value_type &Value) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Compare(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, Value);
    } else {
      *EntriesIter = {Key, Value};
    }
    return EntriesIter->Value();
  }

  value_type &Insert(const key_type &Key, value_type &&Value) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Compare(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::move(Value));
    } else {
      *EntriesIter = {Key, std::move(Value)};
    }
    return EntriesIter->Value();
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> value_type &Insert(const
    key_type &Key, Args &&... Arguments) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Compare(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::forward<Args>(Arguments)...);
    } else {
      *EntriesIter = {Key, std::forward<Args>(Arguments)...};
    }
    return EntriesIter->Value();
  }

  iterator Insert(const_iterator LowerBoundIter, const key_type &Key, const value_type &Value) {
    index_type iEntry = index_type(LowerBoundIter - Begin());
    auto KeysIter = Keys_.Begin() + iEntry;
    auto EntriesIter = Entries_.Begin() + iEntry;
    if (KeysIter == Keys_.End() || Keys_.Compare(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, Value);
    } else {
      *EntriesIter = {Key, Value};
    }
    return iterator(Entries_.Data() + (EntriesIter - Entries_.Begin()));
  }

  iterator Insert(const_iterator LowerBoundIter, const key_type &Key, value_type &&Value) {
    index_type iEntry = index_type(LowerBoundIter - Begin());
    auto KeysIter = Keys_.Begin() + iEntry;
    auto EntriesIter = Entries_.Begin() + iEntry;
    if (KeysIter == Keys_.End() || Keys_.Compare(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::move(Value));
    } else {
      *EntriesIter = {Key, std::move(Value)};
    }
    return iterator(Entries_.Data() + (EntriesIter - Entries_.Begin()));
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> iterator Insert(
    const_iterator LowerBoundIter, const key_type &Key, Args &&...  Arguments) {
    index_type iEntry = index_type(LowerBoundIter - Begin());
    auto KeysIter = Keys_.Begin() + iEntry;
    auto EntriesIter = Entries_.Begin() + iEntry;
    if (KeysIter == Keys_.End() || Keys_.Compare(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::forward<Args>(Arguments)...);
    } else {
      *EntriesIter = {Key, std::forward<Args>(Arguments)...};
    }
    return iterator(Entries_.Data() + (EntriesIter - Entries_.Begin()));
  }

  void Erase(const key_type &Key) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter != Keys_.End() && !Keys_.Compare(Key, *KeysIter)) {
      Keys_.Erase(KeysIter);
      Entries_.Erase(EntriesIter);
    }
  }

  iterator Erase(const_iterator Pos) {
    index_type iEntry = index_type(Pos - Begin());
    Keys_.Erase(Keys_.Begin()+iEntry);
    auto EntriesIter = Entries_.Erase(Entries_.Begin()+iEntry);
    return iterator(Entries_.Data() + (EntriesIter - Entries_.Begin()));
  }

  template <typename F, OVK_FUNCTION_REQUIRES(core::IsCallableAs<F &&, bool(const entry &)>())> void
    EraseIf(F &&Predicate) {
    auto KeysIter = Keys_.Begin();
    auto EntriesIter = Entries_.Begin();
    while (EntriesIter != Entries_.End()) {
      if (std::forward<F>(Predicate)(*EntriesIter)) {
        KeysIter = Keys_.Erase(KeysIter);
        EntriesIter = Entries_.Erase(EntriesIter);
      } else {
        ++KeysIter;
        ++EntriesIter;
      }
    }
  }

  template <typename F, OVK_FUNCTION_REQUIRES(!core::IsCallableAs<F &&, bool(const entry &)>() &&
    core::IsCallableAs<F &&, bool(const key_type &)>())> void EraseIf(F &&Predicate) {
    auto KeysIter = Keys_.Begin();
    auto EntriesIter = Entries_.Begin();
    while (KeysIter != Keys_.End()) {
      if (std::forward<F>(Predicate)(*KeysIter)) {
        KeysIter = Keys_.Erase(KeysIter);
        EntriesIter = Entries_.Erase(EntriesIter);
      } else {
        ++KeysIter;
        ++EntriesIter;
      }
    }
  }

  void Clear() {
    Keys_.Clear();
    Entries_.Clear();
  }

  void Reserve(index_type Count) {
    Keys_.Reserve(Count);
    Entries_.Reserve(Count);
  }

  bool Contains(const key_type &Key) const {
    return Keys_.Contains(Key);
  }

  const_iterator Find(const key_type &Key) const {
    auto KeysIter = Keys_.Find(Key);
    return const_iterator(Entries_.Data() + (KeysIter - Keys_.Begin()));
  }
  iterator Find(const key_type &Key) {
    auto KeysIter = Keys_.Find(Key);
    return iterator(Entries_.Data() + (KeysIter - Keys_.Begin()));
  }

  const_iterator LowerBound(const key_type &Key) const {
    auto KeysIter = Keys_.LowerBound(Key);
    return const_iterator(Entries_.Data() + (KeysIter - Keys_.Begin()));
  }
  iterator LowerBound(const key_type &Key) {
    auto KeysIter = Keys_.LowerBound(Key);
    return iterator(Entries_.Data() + (KeysIter - Keys_.Begin()));
  }

  const_iterator UpperBound(const key_type &Key) const {
    auto KeysIter = Keys_.UpperBound(Key);
    return const_iterator(Entries_.Data() + (KeysIter - Keys_.Begin()));
  }
  iterator UpperBound(const key_type &Key) {
    auto KeysIter = Keys_.UpperBound(Key);
    return iterator(Entries_.Data() + (KeysIter - Keys_.Begin()));
  }

  index_type Count() const { return Entries_.Count(); }

  bool Empty() const { return Entries_.Empty(); }

  index_type Capacity() const { return Entries_.Capacity(); }

  const key_set_type &Keys() const { return Keys_; }

  const value_type &operator()(const key_type &Key) const {
    auto KeysIter = Keys_.LowerBound(Key);
    return Entries_(KeysIter - Keys_.Begin()).Value();
  }

  value_type &operator()(const key_type &Key) {
    auto KeysIter = Keys_.LowerBound(Key);
    return Entries_(KeysIter - Keys_.Begin()).Value();
  }

  value_type &Fetch(const key_type &Key, const value_type &InsertValue) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Compare(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, InsertValue);
    }
    return EntriesIter->Value();
  }

  value_type &Fetch(const key_type &Key, value_type &&InsertValue) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Compare(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::move(InsertValue));
    }
    return EntriesIter->Value();
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> value_type &Fetch(const
    key_type &Key, Args &&... Arguments) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Compare(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::forward<Args>(Arguments)...);
    }
    return EntriesIter->Value();
  }

  const entry &operator[](index_type Index) const { return Entries_(Index); }
  entry &operator[](index_type Index) { return Entries_(Index); }

  const entry *Data() const { return Entries_.Data(); }
  entry *Data() { return Entries_.Data(); }

  const_iterator Begin() const { return const_iterator(Entries_.Data()); }
  iterator Begin() { return iterator(Entries_.Data()); }
  const_iterator CBegin() const { return const_iterator(Entries_.Data()); }

  const_iterator End() const { return const_iterator(Entries_.Data() + Entries_.Count()); }
  iterator End() { return iterator(Entries_.Data() + Entries_.Count()); }
  const_iterator CEnd() const { return const_iterator(Entries_.Data() + Entries_.Count()); }

  const key_compare_type &KeyCompare() const { return Keys_.Compare(); }

  bool KeyCompare(const value_type &Left, const value_type &Right) const {
    return Keys_.Compare(Left, Right);
  }

private:

  key_set_type Keys_;
  array<entry> Entries_;

  friend class core::test_helper<map>;

};

template <typename KeyType, typename ValueType, typename KeyCompareType, bool Contiguous> typename
  map<KeyType, ValueType, KeyCompareType, Contiguous>::iterator begin(map<KeyType, ValueType,
  KeyCompareType, Contiguous> &Map) {
  return Map.Begin();
}

template <typename KeyType, typename ValueType, typename KeyCompareType, bool Contiguous> typename
  map<KeyType, ValueType, KeyCompareType, Contiguous>::const_iterator begin(const map<KeyType,
  ValueType, KeyCompareType, Contiguous> &Map) {
  return Map.Begin();
}

template <typename KeyType, typename ValueType, typename KeyCompareType, bool Contiguous> typename
  map<KeyType, ValueType, KeyCompareType, Contiguous>::iterator end(map<KeyType, ValueType,
  KeyCompareType, Contiguous> &Map) {
  return Map.End();
}

template <typename KeyType, typename ValueType, typename KeyCompareType, bool Contiguous> typename
  map<KeyType, ValueType, KeyCompareType, Contiguous>::const_iterator end(const map<KeyType,
  ValueType, KeyCompareType, Contiguous> &Map) {
  return Map.End();
}

template <typename KeyType, typename ValueType, typename KeyCompareType, bool Contiguous> struct
  array_traits<map<KeyType, ValueType, KeyCompareType, Contiguous>> {
  using map_type = map<KeyType, ValueType, KeyCompareType, Contiguous>;
  using value_type = typename map_type::entry;
  static constexpr int Rank = 1;
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  template <int> static long long ExtentBegin(const map_type &) { return 0; }
  template <int> static long long ExtentEnd(const map_type &Map) {
    return Map.Count();
  }
  static const value_type *Data(const map_type &Map) { return Map.Data(); }
  static value_type *Data(map_type &Map) { return Map.Data(); }
};

template <typename KeyType, typename ValueType> using map_entry_contig = map_entry<KeyType,
  ValueType, true>;
template <typename KeyType, typename ValueType> using map_entry_noncontig = map_entry<KeyType,
  ValueType, false>;

template <typename KeyType, typename ValueType, typename KeyCompareType=std::less<KeyType>> using
  map_contig = map<KeyType, ValueType, KeyCompareType, true>;
template <typename KeyType, typename ValueType, typename KeyCompareType=std::less<KeyType>> using
  map_noncontig = map<KeyType, ValueType, KeyCompareType, false>;

}

#endif
