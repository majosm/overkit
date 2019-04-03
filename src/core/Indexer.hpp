// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_INDEXER_HPP_INCLUDED
#define OVK_CORE_INDEXER_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IntegerSequence.hpp>
#include <ovk/core/Interval.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/TypeSequence.hpp>

#include <cstddef>
#include <type_traits>

namespace ovk {

namespace indexer_internal {

// Can't specialize inside class, so abuse overloads instead
template <array_layout Layout> struct layout_tag {};
template <std::size_t Index> struct index_tag {};

template <typename IndexType, typename TupleElementType, int Rank, array_layout Layout,
  typename TupleElementTypeSequence> class indexer_base;

template <typename IndexType, typename TupleElementType, int Rank_, array_layout Layout_,
  typename... TupleElementTypes> class indexer_base<IndexType, TupleElementType, Rank_, Layout_,
  core::type_sequence<TupleElementTypes...>> {

public:

  static_assert(std::is_integral<IndexType>::value, "Index type must be an integer type.");
  static_assert(std::is_integral<TupleElementType>::value, "Tuple element type must be an integer "
    "type.");

  using index_type = IndexType;
  using tuple_element_type = TupleElementType;
  static constexpr const int Rank = Rank_;
  static constexpr const array_layout Layout = Layout_;
  using tuple_type = elem<tuple_element_type,Rank>;
  using stride_type = elem<index_type,Rank>;
  using interval_type = interval<tuple_element_type,Rank>;

protected:

  tuple_type Begin_;
  stride_type Stride_;

public:

  constexpr OVK_FORCE_INLINE indexer_base():
    Begin_(MakeUniformElem<tuple_element_type,Rank>(0)),
    Stride_(MakeStride_(core::index_sequence_of_size<Rank>(), Begin_, MakeUniformElem<
      tuple_element_type,Rank>(0)))
  {}

  constexpr OVK_FORCE_INLINE explicit indexer_base(const tuple_type &Size):
    Begin_(MakeUniformElem<tuple_element_type,Rank>(0)),
    Stride_(MakeStride_(core::index_sequence_of_size<Rank>(), Begin_, Size))
  {}

  constexpr OVK_FORCE_INLINE explicit indexer_base(const interval_type &Extents):
    Begin_(Extents.Begin()),
    Stride_(MakeStride_(core::index_sequence_of_size<Rank>(), Extents.Begin(), Extents.End()))
  {}

  const tuple_type &Begin() const { return Begin_; }
  const stride_type &Stride() const { return Stride_; }

  // Want to use iterator directly here instead of constructing intermediate elem type
  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsRandomAccessIterator<IterType>() &&
    std::is_convertible<core::iterator_deref_type<IterType>, tuple_element_type>::value)> constexpr
    OVK_FORCE_INLINE index_type ToIndex(IterType First) const {
    return TupleToIndex_(core::index_sequence_of_size<Rank>(), First);
  }

  // Want to use array directly here instead of constructing intermediate elem type
  template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
    !core::IsIterator<ArrayType>() && std::is_convertible<core::array_access_type<const ArrayType
    &>, tuple_element_type>::value && core::ArrayRank<ArrayType>() == 1 &&
    (core::ArrayHasRuntimeExtents<ArrayType>() || (core::StaticArrayHasBegin<ArrayType,0>() &&
    core::StaticArrayHasEnd<ArrayType,Rank>())))> constexpr OVK_FORCE_INLINE index_type ToIndex(
    const ArrayType &Array) const {
    return TupleToIndex_(core::index_sequence_of_size<Rank>(), Array);
  }

  constexpr OVK_FORCE_INLINE index_type ToIndex(TupleElementTypes... TupleElements) const {
    return TupleElementsToIndex_(core::index_sequence_of_size<Rank>(), TupleElements...);
  }

  OVK_FORCE_INLINE tuple_type ToTuple(index_type Index) const {
    return IndexToTuple_(layout_tag<Layout>(), Index);
  }

private:

  template <std::size_t Index, typename BeginType, typename EndType> constexpr OVK_FORCE_INLINE
    index_type MakeStrideElement_(layout_tag<array_layout::ROW_MAJOR> LayoutTag,
    index_tag<Index>, const BeginType &Begin, const EndType &End) const {
    return index_type(End[Index+1]-Begin[Index+1])*MakeStrideElement_(LayoutTag,
      index_tag<Index+1>(), Begin, End);
  }
  template <typename BeginType, typename EndType> constexpr OVK_FORCE_INLINE index_type
    MakeStrideElement_(layout_tag<array_layout::ROW_MAJOR>, index_tag<Rank-1>, const BeginType &,
    const EndType &) const {
    return 1;
  }

  template <std::size_t Index, typename BeginType, typename EndType> constexpr OVK_FORCE_INLINE
    index_type MakeStrideElement_(layout_tag<array_layout::COLUMN_MAJOR> LayoutTag,
    index_tag<Index>, const BeginType &Begin, const EndType &End) const {
    return index_type(End[Index-1]-Begin[Index-1])*MakeStrideElement_(LayoutTag,
      index_tag<Index-1>(), Begin, End);
  }
  template <typename BeginType, typename EndType> constexpr OVK_FORCE_INLINE index_type
    MakeStrideElement_(layout_tag<array_layout::COLUMN_MAJOR>, index_tag<0>, const BeginType &,
    const EndType &) const {
    return 1;
  }

  template <std::size_t... Indices, typename BeginType, typename EndType> constexpr OVK_FORCE_INLINE
    stride_type MakeStride_(core::index_sequence<Indices...>, const BeginType &Begin, const EndType
    &End) const {
    return {MakeStrideElement_(layout_tag<Layout>(), index_tag<Indices>(), Begin, End)...};
  }

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices, typename
    TupleType> constexpr index_type OVK_FORCE_INLINE TupleToIndex_(core::index_sequence<Index1,
    Index2, RemainingIndices...>, const TupleType &Tuple) const {
    return Stride_[Index1]*(Tuple[Index1]-Begin_[Index1]) + TupleToIndex_(core::index_sequence<
      Index2, RemainingIndices...>(), Tuple);
  }

  template <std::size_t Index, typename TupleType> constexpr index_type OVK_FORCE_INLINE
    TupleToIndex_(core::index_sequence<Index>, const TupleType &Tuple) const {
    return Stride_[Index]*(Tuple[Index]-Begin_[Index]);
  }

  template <std::size_t Index1, std::size_t Index2, std::size_t... RemainingIndices, typename...
    RemainingTupleElementTypes> constexpr index_type OVK_FORCE_INLINE TupleElementsToIndex_(
    core::index_sequence<Index1, Index2, RemainingIndices...>, TupleElementType TupleElement1,
    RemainingTupleElementTypes... RemainingTupleElements) const {
    return Stride_[Index1]*(TupleElement1-Begin_[Index1]) + TupleElementsToIndex_(
      core::index_sequence<Index2, RemainingIndices...>(), RemainingTupleElements...);
  }

  template <std::size_t Index> constexpr index_type OVK_FORCE_INLINE TupleElementsToIndex_(
    core::index_sequence<Index>, TupleElementType TupleElement) const {
    return Stride_[Index]*(TupleElement-Begin_[Index]);
  }

  tuple_type OVK_FORCE_INLINE IndexToTuple_(layout_tag<array_layout::ROW_MAJOR>, index_type Index)
    const {

    tuple_type Tuple;

    index_type ReducedIndex = Index;
    for (int i = 0; i < Rank; ++i) {
      tuple_element_type Offset = tuple_element_type(ReducedIndex/Stride_[i]);
      Tuple[i] = Begin_[i] + Offset;
      ReducedIndex -= Offset*Stride_[i];
    }

    return Tuple;

  }

  tuple_type OVK_FORCE_INLINE IndexToTuple_(layout_tag<array_layout::COLUMN_MAJOR>, index_type
    Index) const {

    tuple_type Tuple;

    index_type ReducedIndex = Index;
    for (int i = Rank-1; i >= 0; --i) {
      tuple_element_type Offset = tuple_element_type(ReducedIndex/Stride_[i]);
      Tuple[i] = Begin_[i] + Offset;
      ReducedIndex -= Offset*Stride_[i];
    }

    return Tuple;

  }

};

}

template <typename IndexType, typename TupleElementType, int Rank_, array_layout Layout_=
  array_layout::ROW_MAJOR> class indexer : public indexer_internal::indexer_base<IndexType,
  TupleElementType, Rank_, Layout_, core::repeated_type_sequence_of_size<TupleElementType,Rank_>> {

private:

  using parent_type = indexer_internal::indexer_base<IndexType, TupleElementType, Rank_, Layout_,
    core::repeated_type_sequence_of_size<TupleElementType,Rank_>>;

public:

  using typename parent_type::index_type;
  using typename parent_type::tuple_element_type;
  using parent_type::Rank;
  using parent_type::Layout;
  using typename parent_type::tuple_type;
  using typename parent_type::interval_type;
  using parent_type::parent_type;
  using parent_type::ToIndex;
  using parent_type::ToTuple;

  friend class core::test_helper<indexer>;

};

template <typename IndexType, typename TupleElementType, int Rank, array_layout Layout> constexpr
  OVK_FORCE_INLINE bool operator==(const indexer<IndexType, TupleElementType, Rank, Layout> &Left,
  const indexer<IndexType, TupleElementType, Rank, Layout> &Right) {
  return Left.Begin() == Right.Begin() && Left.Stride() == Right.Stride();
}

template <typename IndexType, typename TupleElementType, int Rank, array_layout Layout> constexpr
  OVK_FORCE_INLINE bool operator!=(const indexer<IndexType, TupleElementType, Rank, Layout> &Left,
  const indexer<IndexType, TupleElementType, Rank, Layout> &Right) {
  return !(Left == Right);
}

}

#endif
