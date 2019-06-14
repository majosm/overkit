// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_INTERVAL_HPP_INCLUDED
#define OVK_CORE_INTERVAL_HPP_INCLUDED

#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/ScalarOps.hpp>

#include <type_traits>

namespace ovk {

template <typename T, int N> class interval_base_1 {

public:

  using value_type = T;
  static constexpr int Rank = N;
  using tuple_type = elem<T,N>;

protected:

  tuple_type Begin_;
  tuple_type End_;

public:

  constexpr interval_base_1() = default;

  constexpr interval_base_1(const tuple_type &Size):
    Begin_(MakeUniformElem<value_type,Rank>(value_type(0))),
    End_(Size)
  {}

  constexpr interval_base_1(const tuple_type &Begin, const tuple_type &End):
    Begin_(Begin),
    End_(End)
  {}

  template <typename U, OVK_FUNCTION_REQUIRES(!std::is_same<U, value_type>::value &&
    std::is_convertible<U, value_type>::value)> constexpr interval_base_1(const interval_base_1<U,N>
    &Other):
    Begin_(Other.Begin_),
    End_(Other.End_)
  {}

  constexpr OVK_FORCE_INLINE const tuple_type &Begin() const { return Begin_; }
  OVK_FORCE_INLINE tuple_type &Begin() { return Begin_; }
  constexpr OVK_FORCE_INLINE value_type Begin(int iDim) const { return Begin_(iDim); }
  OVK_FORCE_INLINE value_type &Begin(int iDim) { return Begin_(iDim); }

  constexpr OVK_FORCE_INLINE const tuple_type &End() const { return End_; }
  OVK_FORCE_INLINE tuple_type &End() { return End_; }
  constexpr OVK_FORCE_INLINE value_type End(int iDim) const { return End_(iDim); }
  OVK_FORCE_INLINE value_type &End(int iDim) { return End_(iDim); }

  constexpr OVK_FORCE_INLINE tuple_type Size() const {
    return Max(End_ - Begin_, MakeUniformElem<value_type,N>(value_type(0)));
  }
  constexpr OVK_FORCE_INLINE value_type Size(int iDim) const {
    return Max(End_(iDim) - Begin_(iDim), value_type(0));
  }

protected:

  template <typename ResultType=value_type, std::size_t Index1, std::size_t Index2, std::size_t...
    RemainingIndices> constexpr OVK_FORCE_INLINE ResultType VolumeCountHelper_(core::index_sequence<
    Index1, Index2, RemainingIndices...>) const {
    return ResultType(Size(Index1)) * VolumeCountHelper_<ResultType>(core::index_sequence<Index2,
      RemainingIndices...>());
  }
  template <typename ResultType=value_type, std::size_t Index> constexpr OVK_FORCE_INLINE
    value_type VolumeCountHelper_(core::index_sequence<Index>) const {
    return ResultType(Size(Index));
  }

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices>
    constexpr OVK_FORCE_INLINE bool IncludesHelper_(core::index_sequence<Index1, Index2,
    RemainingIndices...>, const interval_base_1 &Other) const {
    return (Begin_(Index1) <= Other.Begin_(Index1) && End_(Index1) >= Other.End_(Index1)) &&
      IncludesHelper_(core::index_sequence<Index2, RemainingIndices...>(), Other);
  }
  template <std::size_t Index> constexpr OVK_FORCE_INLINE bool IncludesHelper_(core::index_sequence<
    Index>, const interval_base_1 &Other) const {
    return Begin_(Index) <= Other.Begin_(Index) && End_(Index) >= Other.End_(Index);
  }

  template <typename U, int M> friend class interval_base_1;

};

template <typename T, int N, typename=void> class interval_base_2;

template <typename T, int N> class interval_base_2<T, N, OVK_SPECIALIZATION_REQUIRES(
  std::is_integral<T>::value)> : public interval_base_1<T,N> {

private:

  using parent_type = interval_base_1<T,N>;
  using parent_type::VolumeCountHelper_;
  using parent_type::IncludesHelper_;

protected:

  using parent_type::Begin_;
  using parent_type::End_;

public:

  using typename parent_type::value_type;
  using parent_type::Rank;
  using typename parent_type::tuple_type;
  using parent_type::parent_type;

  template <typename IndexType=long long> constexpr OVK_FORCE_INLINE IndexType Count() const {
    return parent_type::template VolumeCountHelper_<IndexType>(core::index_sequence_of_size<N>());
  }

  template <typename... Args> value_type Volume(Args &&...) const {
    // Assert on Args instead of just 'false' so it doesn't trigger unless method is
    // instantiated
    static_assert(int(sizeof...(Args)) < 0, "Cannot use interval::Volume for integral value types.");
    return value_type(0);
  }

  constexpr OVK_FORCE_INLINE bool Empty() const {
    return Empty_(core::index_sequence_of_size<N>());
  }

  constexpr OVK_FORCE_INLINE bool Contains(const tuple_type &Tuple) const {
    return Contains_(core::index_sequence_of_size<N>(), Tuple);
  }

  constexpr OVK_FORCE_INLINE bool Includes(const interval_base_2 &Other) const {
    return Other.Empty() || IncludesHelper_(core::index_sequence_of_size<N>(), Other);
  }

private:

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices>
    constexpr OVK_FORCE_INLINE bool Empty_(core::index_sequence<Index1, Index2,
    RemainingIndices...>) const {
    return (End_(Index1) <= Begin_(Index1)) || Empty_(core::index_sequence<Index2,
      RemainingIndices...>());
  }
  template <std::size_t Index> constexpr OVK_FORCE_INLINE bool Empty_(core::index_sequence<Index>)
    const {
    return End_(Index) <= Begin_(Index);
  }

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices>
    constexpr OVK_FORCE_INLINE bool Contains_(core::index_sequence<Index1, Index2,
    RemainingIndices...>, const tuple_type &Tuple) const {
    return (Tuple(Index1) >= Begin_(Index1) && Tuple(Index1) < End_(Index1)) && Contains_(
      core::index_sequence<Index2, RemainingIndices...>(), Tuple);
  }
  template <std::size_t Index> constexpr OVK_FORCE_INLINE bool Contains_(core::index_sequence<Index>,
    const tuple_type &Tuple) const {
    return Tuple(Index) >= Begin_(Index) && Tuple(Index) < End_(Index);
  }

};

template <typename T, int N> class interval_base_2<T, N, OVK_SPECIALIZATION_REQUIRES(
  std::is_floating_point<T>::value)> : public interval_base_1<T,N> {

private:

  using parent_type = interval_base_1<T,N>;
  using parent_type::VolumeCountHelper_;
  using parent_type::IncludesHelper_;

protected:

  using parent_type::Begin_;
  using parent_type::End_;

public:

  using typename parent_type::value_type;
  using parent_type::Rank;
  using typename parent_type::tuple_type;
  using parent_type::parent_type;

  template <typename IndexType=long long> IndexType Count() const {
    // Assert on IndexType instead of just 'false' so it doesn't trigger unless method is
    // instantiated
    static_assert(std::is_same<IndexType,void>::value, "Cannot use interval::Count for floating "
      "point value types.");
    return IndexType(0);
  }

  constexpr OVK_FORCE_INLINE value_type Volume() const {
    return VolumeCountHelper_(core::index_sequence_of_size<N>());
  }

  constexpr OVK_FORCE_INLINE bool Empty() const {
    return Empty_(core::index_sequence_of_size<N>());
  }

  constexpr OVK_FORCE_INLINE bool Contains(const tuple_type &Tuple) const {
    return Contains_(core::index_sequence_of_size<N>(), Tuple);
  }

  constexpr OVK_FORCE_INLINE bool Includes(const interval_base_2 &Other) const {
    return Other.Empty() || IncludesHelper_(core::index_sequence_of_size<N>(), Other);
  }

private:

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices>
    constexpr OVK_FORCE_INLINE bool Empty_(core::index_sequence<Index1, Index2,
    RemainingIndices...>) const {
    return (End_(Index1) < Begin_(Index1)) || Empty_(core::index_sequence<Index2,
      RemainingIndices...>());
  }
  template <std::size_t Index> constexpr OVK_FORCE_INLINE bool Empty_(core::index_sequence<Index>)
    const {
    return End_(Index) < Begin_(Index);
  }

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices>
    constexpr OVK_FORCE_INLINE bool Contains_(core::index_sequence<Index1, Index2,
    RemainingIndices...>, const tuple_type &Tuple) const {
    return (Tuple(Index1) >= Begin_(Index1) && Tuple(Index1) <= End_(Index1)) && Contains_(
      core::index_sequence<Index2, RemainingIndices...>(), Tuple);
  }
  template <std::size_t Index> constexpr OVK_FORCE_INLINE bool Contains_(core::index_sequence<Index>,
    const tuple_type &Tuple) const {
    return Tuple(Index) >= Begin_(Index) && Tuple(Index) <= End_(Index);
  }

};

template <typename T, int N=1> class interval : public interval_base_2<T,N> {

  static_assert(std::is_arithmetic<T>::value, "Interval value type must be an arithmetic type.");

private:

  using parent_type = interval_base_2<T,N>;
  using parent_type::Begin_;
  using parent_type::End_;

public:

  using typename parent_type::value_type;
  using parent_type::Rank;
  using typename parent_type::tuple_type;
  using parent_type::parent_type;
  using parent_type::Count;
  using parent_type::Volume;
  using parent_type::Empty;
  using parent_type::Contains;

private:

  friend class core::test_helper<interval>;

};

template <typename T, int N, OVK_FUNCTION_REQUIRES(std::is_integral<T>::value)> constexpr
  interval<T,N> MakeEmptyInterval() {
  return {MakeUniformElem<T,N>(T(0))};
}

template <typename T, int N, OVK_FUNCTION_REQUIRES(std::is_floating_point<T>::value)> constexpr
  interval<T,N> MakeEmptyInterval() {
  return {MakeUniformElem<T,N>(T(-1))};
}

template <typename T, int N> constexpr bool operator==(const interval<T,N> &Left, const interval<T,
  N> &Right) {
  return Left.Begin() == Right.Begin() && Left.End() == Right.End();
}
template <typename T, int N> constexpr bool operator!=(const interval<T,N> &Left, const interval<T,
  N> &Right) {
  return !(Left == Right);
}

}

#endif
