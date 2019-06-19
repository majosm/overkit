// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ID_SET_HPP_INCLUDED
#define OVK_CORE_ID_SET_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/Requires.hpp>

#include <algorithm>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <utility>

namespace ovk {

namespace id_set_internal {

template <int Rank, typename IDTypeSequence> class id_set_base;

template <int Rank, typename... IDTypes> class id_set_base<Rank, core::type_sequence<IDTypes...>> {

public:

  using value_type = elem<int,Rank>;
  // Values aren't mutable, so iterator is const
  using iterator = const value_type *;
  using const_iterator = iterator;

  void Insert(IDTypes... IDs) {
    value_type Value(IDs...);
    auto ValuesIter = LowerBound_(Value);
    if (ValuesIter == Values_.End() || Less(Value, *ValuesIter)) {
      ValuesIter = Values_.Insert(ValuesIter, Value);
    }
  }

  void Erase(IDTypes... IDs) {
    value_type Value(IDs...);
    auto ValuesIter = LowerBound_(Value);
    if (ValuesIter != Values_.End() && !Less(Value, *ValuesIter)) {
      Values_.Erase(ValuesIter);
    }
  }

  bool Contains(IDTypes... IDs) const {
    value_type Value(IDs...);
    auto ValuesIter = LowerBound_(Value);
    return ValuesIter != Values_.End() && !Less(Value, *ValuesIter);
  }

  iterator Find(IDTypes... IDs) const {
    value_type Value(IDs...);
    auto ValuesIter = LowerBound_(Value);
    if (ValuesIter != Values_.End() && Less(Value, *ValuesIter)) {
      ValuesIter = Values_.End();
    }
    return Values_.Data() + (ValuesIter - Values_.Begin());
  }

  iterator LowerBound(IDTypes... IDs) const {
    value_type Value(IDs...);
    auto ValuesIter = LowerBound_(Value);
    return Values_.Data() + (ValuesIter - Values_.Begin());
  }

  iterator UpperBound(IDTypes... IDs) const {
    value_type Value(IDs...);
    auto ValuesIter = UpperBound_(Value);
    return Values_.Data() + (ValuesIter - Values_.Begin());
  }

  static constexpr bool Less(const value_type &Left, const value_type &Right) {
    return Less_(Left, Right, index_tag<0>());
  }

protected:

  template <int N> struct index_tag {};

  array<value_type> Values_;

  typename array<value_type>::const_iterator LowerBound_(const value_type &Value) const {
    auto Compare = [](const value_type &Left, const value_type &Right) -> bool {
      return Less(Left, Right);
    };
    return std::lower_bound(Values_.Begin(), Values_.End(), Value, Compare);
  }

  typename array<value_type>::iterator LowerBound_(const value_type &Value) {
    auto Compare = [](const value_type &Left, const value_type &Right) -> bool {
      return Less(Left, Right);
    };
    return std::lower_bound(Values_.Begin(), Values_.End(), Value, Compare);
  }

  typename array<value_type>::const_iterator UpperBound_(const value_type &Value) const {
    auto Compare = [](const value_type &Left, const value_type &Right) -> bool {
      return Less(Left, Right);
    };
    return std::upper_bound(Values_.Begin(), Values_.End(), Value, Compare);
  }

  typename array<value_type>::iterator UpperBound_(const value_type &Value) {
    auto Compare = [](const value_type &Left, const value_type &Right) -> bool {
      return Less(Left, Right);
    };
    return std::upper_bound(Values_.Begin(), Values_.End(), Value, Compare);
  }

  template <int N> static constexpr OVK_FORCE_INLINE bool Less_(const value_type &Left, const
    value_type &Right, index_tag<N>) {
    return Left(N) < Right(N) || (Left(N) == Right(N) && Less_(Left, Right, index_tag<N+1>()));
  }

  static constexpr OVK_FORCE_INLINE bool Less_(const value_type &Left, const value_type &Right,
    index_tag<Rank-1>) {
    return Left(Rank-1) < Right(Rank-1);
  }

};

}

template <int Rank_> class id_set : public id_set_internal::id_set_base<Rank_,
  core::repeated_type_sequence_of_size<int, Rank_>> {

private:

  using parent_type = id_set_internal::id_set_base<Rank_, core::repeated_type_sequence_of_size<int,
    Rank_>>;

  using parent_type::Values_;
  using parent_type::LowerBound_;
  using parent_type::UpperBound_;

public:

  static constexpr int Rank = Rank_;
  using value_type = typename parent_type::value_type;
  using iterator = typename parent_type::iterator;
  using const_iterator = typename parent_type::const_iterator;

  using parent_type::Insert;
  using parent_type::Erase;
  using parent_type::Contains;
  using parent_type::Find;
  using parent_type::LowerBound;
  using parent_type::UpperBound;
  using parent_type::Less;

  id_set() = default;

  id_set(std::initializer_list<value_type> ValuesList) {
    Values_.Reserve(ValuesList.size());
    for (auto &Value : ValuesList) {
      Insert(Value);
    }
  }

  // Don't care about input iterators enough to implement a second overload
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>())>
    id_set(IterType Begin, IterType End) {
    Values_.Reserve(std::distance(Begin, End));
    IterType Iter = Begin;
    while (Iter != End) {
      Insert(*Iter);
      ++Iter;
    }
  }

  id_set &operator=(std::initializer_list<value_type> ValuesList) {
    return Assign(ValuesList);
  }

  id_set &Assign(std::initializer_list<value_type> ValuesList) {
    Values_.Clear();
    Values_.Reserve(ValuesList.size());
    for (auto &Value : ValuesList) {
      Insert(Value);
    }
    return *this;
  }

  // Don't care about input iterators enough to implement a second overload
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsForwardIterator<IterType>())>
    id_set &Assign(IterType Begin, IterType End) {
    Values_.Clear();
    Values_.Reserve(std::distance(Begin, End));
    IterType Iter = Begin;
    while (Iter != End) {
      Insert(*Iter);
      ++Iter;
    }
    return *this;
  }

  void Reserve(long long Count) {
    Values_.Reserve(Count);
  }

  void Insert(const value_type &Value) {
    auto ValuesIter = LowerBound_(Value);
    if (ValuesIter == Values_.End() || Less(Value, *ValuesIter)) {
      ValuesIter = Values_.Insert(ValuesIter, Value);
    }
  }

  void Insert(iterator LowerBoundIter, const value_type &Value) {
    auto ValuesIter = Values_.Begin() + (LowerBoundIter - Values_.Data());
    if (ValuesIter == Values_.End() || Less(Value, *ValuesIter)) {
      ValuesIter = Values_.Insert(ValuesIter, Value);
    }
  }

  void Erase(const value_type &Value) {
    auto ValuesIter = LowerBound_(Value);
    if (ValuesIter != Values_.End() && !Less(Value, *ValuesIter)) {
      Values_.Erase(ValuesIter);
    }
  }

  iterator Erase(iterator Pos) {
    auto ValuesIter = Values_.Begin() + (Pos - Values_.Data());
    ValuesIter = Values_.Erase(ValuesIter);
    return Values_.Data(ValuesIter - Values_.Begin());
  }

  void Clear() {
    Values_.Clear();
  }

  bool Contains(const value_type &Value) const {
    auto ValuesIter = LowerBound_(Value);
    return ValuesIter != Values_.End() && !Less(Value, *ValuesIter);
  }

  iterator Find(const value_type &Value) const {
    auto ValuesIter = LowerBound_(Value);
    if (ValuesIter != Values_.End() && Less(Value, *ValuesIter)) {
      ValuesIter = Values_.End();
    }
    return Values_.Data() + (ValuesIter - Values_.Begin());
  }

  iterator LowerBound(const value_type &Value) const {
    auto ValuesIter = LowerBound_(Value);
    return Values_.Data() + (ValuesIter - Values_.Begin());
  }

  iterator UpperBound(const value_type &Value) const {
    auto ValuesIter = UpperBound_(Value);
    return Values_.Data() + (ValuesIter - Values_.Begin());
  }

  // Not sure if this is actually useful for Rank > 1
  value_type NextAvailableValue(int iDim=0) const {
    value_type Result = MakeUniformElem<int,Rank>(0);
    if (Values_.Count() > 0 && Values_(0)(iDim) == 0) {
      int PrevID = Values_(0)(iDim);
      for (int iValue = 1; iValue < Values_.Count(); ++iValue) {
        int ID = Values_(iValue)(iDim);
        if (ID - PrevID > 1) {
          break;
        }
        PrevID = ID;
      }
      Result(iDim) = PrevID + 1;
    }
    return Result;
  }

  long long Count() const { return Values_.Count(); }

  bool Empty() const { return Values_.Empty(); }

  long long Capacity() const { return Values_.Capacity(); }

  const value_type &operator[](long long Index) const { return Values_(Index); }

  const value_type *Data() const { return Values_.Data(); }

  iterator Begin() const { return Values_.Data(); }

  iterator End() const { return Values_.Data() + Values_.Count(); }

private:

  friend class core::test_helper<id_set>;

};

template <int Rank> typename id_set<Rank>::iterator begin(const id_set<Rank> &Set) {
  return Set.Begin();
}

template <int Rank> typename id_set<Rank>::iterator end(const id_set<Rank> &Set) {
  return Set.End();
}

}

#endif
