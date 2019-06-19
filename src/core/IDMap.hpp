// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ID_MAP_HPP_INCLUDED
#define OVK_CORE_ID_MAP_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IDSet.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeSequence.hpp>

#include <algorithm>
#include <initializer_list>
#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {

namespace id_map_internal {

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

template <int Rank_, typename T, bool Contiguous_> class id_map_entry {

public:

  static constexpr int Rank = Rank_;
  using key_type = elem<int,Rank>;
  using value_type = T;
  static constexpr bool Contiguous = Contiguous_;

  id_map_entry(const key_type &Key, const value_type &Value):
    Key_(Key),
    Value_(Value)
  {}

  id_map_entry(const key_type &Key, value_type &&Value):
    Key_(Key),
    Value_(std::move(Value))
  {}

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> id_map_entry(const key_type
    &Key, Args &&... Arguments):
    Key_(Key),
    Value_(std::forward<Args>(Arguments)...)
  {}

  const key_type &Key() const { return Key_; }
  int Key(int iDim) const { return Key_(iDim); }

  const value_type &Value() const { return Value_.Get(); }
  value_type &Value() { return Value_.Get(); }

private:

  using value_wrapper_type = value_wrapper<value_type, Contiguous>;

  key_type Key_;
  value_wrapper_type Value_;

  friend class core::test_helper<id_map_entry>;

};

template <int Rank, typename T, bool Contiguous, typename IDTypeSequence> class id_map_base;

template <int Rank, typename T, bool Contiguous_, typename... IDTypes> class id_map_base<Rank, T,
  Contiguous_, core::type_sequence<IDTypes...>> {

public:

  using key_type = elem<int,Rank>;
  using value_type = T;
  static constexpr bool Contiguous = Contiguous_;
  using id_set_type = id_set<Rank>;
  using entry = id_map_entry<Rank,T,Contiguous>;
  using iterator = entry *;
  using const_iterator = const entry *;

  const value_type &operator()(IDTypes... IDs) const {
    auto KeysIter = Keys_.LowerBound(IDs...);
    return Entries_(KeysIter - Keys_.Begin()).Value();
  }

  value_type &operator()(IDTypes... IDs) {
    auto KeysIter = Keys_.LowerBound(IDs...);
    return Entries_(KeysIter - Keys_.Begin()).Value();
  }

  value_type &Insert(IDTypes... IDs) {
    key_type Key(IDs...);
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key);
    } else {
      *EntriesIter = {Key};
    }
    return EntriesIter->Value();
  }

  void Erase(IDTypes... IDs) {
    key_type Key(IDs...);
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter != Keys_.End() && !Keys_.Less(Key, *KeysIter)) {
      Keys_.Erase(KeysIter);
      Entries_.Erase(EntriesIter);
    }
  }

  bool Contains(IDTypes... IDs) const {
    return Keys_.Contains(IDs...);
  }

  const_iterator Find(IDTypes... IDs) const {
    auto KeysIter = Keys_.Find(IDs...);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }
  iterator Find(IDTypes... IDs) {
    auto KeysIter = Keys_.Find(IDs...);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }

  const_iterator LowerBound(IDTypes... IDs) const {
    auto KeysIter = Keys_.LowerBound(IDs...);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }
  iterator LowerBound(IDTypes... IDs) {
    auto KeysIter = Keys_.LowerBound(IDs...);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }

  const_iterator UpperBound(IDTypes... IDs) const {
    auto KeysIter = Keys_.UpperBound(IDs...);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }
  iterator UpperBound(IDTypes... IDs) {
    auto KeysIter = Keys_.UpperBound(IDs...);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }

  value_type &Get(IDTypes... IDs) {
    key_type Key(IDs...);
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key);
    }
    return EntriesIter->Value();
  }

protected:

  id_set_type Keys_;
  array<entry> Entries_;

};

}

template <int Rank_, typename T, bool Contiguous_=sizeof(T) <= 16> class id_map : public
  id_map_internal::id_map_base<Rank_, T, Contiguous_, core::repeated_type_sequence_of_size<
  int, Rank_>> {

private:

  using parent_type = id_map_internal::id_map_base<Rank_, T, Contiguous_,
    core::repeated_type_sequence_of_size<int, Rank_>>;

  using parent_type::Keys_;
  using parent_type::Entries_;

public:

  static constexpr int Rank = Rank_;
  using key_type = typename parent_type::key_type;
  using value_type = typename parent_type::value_type;
  static constexpr bool Contiguous = Contiguous_;
  using id_set_type = typename parent_type::id_set_type;
  using entry = typename parent_type::entry;
  using iterator = typename parent_type::iterator;
  using const_iterator = typename parent_type::const_iterator;

  using parent_type::operator();
  using parent_type::Insert;
  using parent_type::Erase;
  using parent_type::Contains;
  using parent_type::Find;
  using parent_type::LowerBound;
  using parent_type::UpperBound;
  using parent_type::Get;

  id_map() = default;

  id_map(std::initializer_list<entry> EntriesList) {
    Keys_.Reserve(EntriesList.size());
    Entries_.Reserve(EntriesList.size());
    for (auto &Entry : EntriesList) {
      Insert(Entry);
    }
  }

  // Don't care about input iterators enough to implement a second overload
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>())>
    id_map(IterType Begin, IterType End) {
    long long NumEntries = (long long)(std::distance(Begin, End));
    Keys_.Reserve(NumEntries);
    Entries_.Reserve(NumEntries);
    IterType Iter = Begin;
    while (Iter != End) {
      Insert(*Iter);
      ++Iter;
    }
  }

  id_map &operator=(std::initializer_list<entry> EntriesList) {
    return Assign(EntriesList);
  }

  id_map &Assign(std::initializer_list<entry> EntriesList) {
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
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>())>
    id_map &Assign(IterType Begin, IterType End) {
    long long NumEntries = (long long)(std::distance(Begin, End));
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
    if (KeysIter == Keys_.End() || Keys_.Less(Entry.Key(), *KeysIter)) {
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
    if (KeysIter == Keys_.End() || Keys_.Less(Entry.Key(), *KeysIter)) {
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
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
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
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
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
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::forward<Args>(Arguments)...);
    } else {
      *EntriesIter = {Key, std::forward<Args>(Arguments)...};
    }
    return EntriesIter->Value();
  }

  value_type &Insert(const_iterator LowerBoundIter, const key_type &Key, const value_type &Value) {
    auto KeysIter = Keys_.Begin() + (LowerBoundIter - Entries_.Data());
    auto EntriesIter = Entries_.Begin() + (LowerBoundIter - Entries_.Data());
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, Value);
    } else {
      *EntriesIter = {Key, Value};
    }
    return EntriesIter->Value();
  }

  value_type &Insert(const_iterator LowerBoundIter, const key_type &Key, value_type &&Value) {
    auto KeysIter = Keys_.Begin() + (LowerBoundIter - Entries_.Data());
    auto EntriesIter = Entries_.Begin() + (LowerBoundIter - Entries_.Data());
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::move(Value));
    } else {
      *EntriesIter = {Key, std::move(Value)};
    }
    return EntriesIter->Value();
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> value_type &Insert(
    const_iterator LowerBoundIter, const key_type &Key, Args &&...  Arguments) {
    auto KeysIter = Keys_.Begin() + (LowerBoundIter - Entries_.Data());
    auto EntriesIter = Entries_.Begin() + (LowerBoundIter - Entries_.Data());
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::forward<Args>(Arguments)...);
    } else {
      *EntriesIter = {Key, std::forward<Args>(Arguments)...};
    }
    return EntriesIter->Value();
  }

  void Erase(const key_type &Key) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter != Keys_.End() && !Keys_.Less(Key, *KeysIter)) {
      Keys_.Erase(KeysIter);
      Entries_.Erase(EntriesIter);
    }
  }

  iterator Erase(const_iterator Pos) {
    auto KeysIter = Keys_.Begin() + (Pos - Entries_.Data());
    auto EntriesIter = Entries_.Begin() + (Pos - Entries_.Data());
    Keys_.Erase(KeysIter);
    EntriesIter = Entries_.Erase(EntriesIter);
    return Entries_.Data(EntriesIter - Entries_.Begin());
  }

  void Clear() {
    Keys_.Clear();
    Entries_.Clear();
  }

  void Reserve(long long Count) {
    Keys_.Reserve(Count);
    Entries_.Reserve(Count);
  }

  bool Contains(const key_type &Key) const {
    return Keys_.Contains(Key);
  }

  const_iterator Find(const key_type &Key) const {
    auto KeysIter = Keys_.Find(Key);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }
  iterator Find(const key_type &Key) {
    auto KeysIter = Keys_.Find(Key);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }

  const_iterator LowerBound(const key_type &Key) const {
    auto KeysIter = Keys_.LowerBound(Key);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }
  iterator LowerBound(const key_type &Key) {
    auto KeysIter = Keys_.LowerBound(Key);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }

  const_iterator UpperBound(const key_type &Key) const {
    auto KeysIter = Keys_.UpperBound(Key);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }
  iterator UpperBound(const key_type &Key) {
    auto KeysIter = Keys_.UpperBound(Key);
    return Entries_.Data() + (KeysIter - Keys_.Begin());
  }

  key_type NextAvailableKey(int iDim=0) const {
    return Keys_.NextAvailableValue(iDim);
  }

  long long Count() const { return Entries_.Count(); }

  bool Empty() const { return Entries_.Empty(); }

  long long Capacity() const { return Entries_.Capacity(); }

  const id_set_type &Keys() const { return Keys_; }

  const value_type &operator()(const key_type &Key) const {
    auto KeysIter = Keys_.LowerBound(Key);
    return Entries_(KeysIter - Keys_.Begin()).Value();
  }

  value_type &operator()(const key_type &Key) {
    auto KeysIter = Keys_.LowerBound(Key);
    return Entries_(KeysIter - Keys_.Begin()).Value();
  }

  value_type &Get(const key_type &Key, const value_type &InsertValue) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, InsertValue);
    }
    return EntriesIter->Value();
  }

  value_type &Get(const key_type &Key, value_type &&InsertValue) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::move(InsertValue));
    }
    return EntriesIter->Value();
  }

  template <typename... Args, OVK_FUNCTION_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> value_type &Get(const
    key_type &Key, Args &&... Arguments) {
    auto KeysIter = Keys_.LowerBound(Key);
    auto EntriesIter = Entries_.Begin() + (KeysIter - Keys_.Begin());
    if (KeysIter == Keys_.End() || Keys_.Less(Key, *KeysIter)) {
      Keys_.Insert(KeysIter, Key);
      EntriesIter = Entries_.Insert(EntriesIter, Key, std::forward<Args>(Arguments)...);
    }
    return EntriesIter->Value();
  }

  const entry &operator[](long long Index) const { return Entries_(Index); }
  entry &operator[](long long Index) { return Entries_(Index); }

  const entry *Data() const { return Entries_.Data(); }
  entry *Data() { return Entries_.Data(); }

  const_iterator Begin() const { return Entries_.Data(); }
  iterator Begin() { return Entries_.Data(); }
  const_iterator CBegin() const { return Entries_.Data(); }

  const_iterator End() const { return Entries_.Data() + Entries_.Count(); }
  iterator End() { return Entries_.Data() + Entries_.Count(); }
  const_iterator CEnd() const { return Entries_.Data() + Entries_.Count(); }

  static constexpr bool KeyLess(const key_type &Left, const key_type &Right) {
    return id_set_type::Less(Left, Right);
  }

private:

  friend class core::test_helper<id_map>;

};

template <typename T, int Rank, bool Contiguous> typename id_map<Rank, T>::iterator begin(
  id_map<Rank, T, Contiguous> &Map) {
  return Map.Begin();
}

template <typename T, int Rank, bool Contiguous> typename id_map<Rank, T>::const_iterator begin(
  const id_map<Rank, T, Contiguous> &Map) {
  return Map.Begin();
}

template <typename T, int Rank, bool Contiguous> typename id_map<Rank, T>::iterator end(id_map<Rank,
  T, Contiguous> &Map) {
  return Map.End();
}

template <typename T, int Rank, bool Contiguous> typename id_map<Rank, T>::const_iterator end(const
  id_map<Rank, T, Contiguous> &Map) {
  return Map.End();
}

}

#endif