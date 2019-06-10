// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_INDEXER_HPP_INCLUDED
#define OVK_CORE_INDEXER_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
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

// Slightly faster than using constexpr recursion (compiles a bit faster too)
// May be able to replace with C++17 fold expressions once they are widely supported
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,1> &Begin, const elem<
  IndexType,1> &, const ArrayOrIterType &Tuple) {
  return IndexType(Tuple[0] - Begin(0));
}
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,2> &Begin, const elem<
  IndexType,2> &Stride, const ArrayOrIterType &Tuple) {
  return
    Stride(0)*IndexType(Tuple[0] - Begin(0)) +
    Stride(1)*IndexType(Tuple[1] - Begin(1));
}
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,3> &Begin, const elem<
  IndexType,3> &Stride, const ArrayOrIterType &Tuple) {
  return
    Stride(0)*IndexType(Tuple[0] - Begin(0)) +
    Stride(1)*IndexType(Tuple[1] - Begin(1)) +
    Stride(2)*IndexType(Tuple[2] - Begin(2));
}
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,4> &Begin,
  const elem<IndexType,4> &Stride, const ArrayOrIterType &Tuple) {
  return
    Stride(0)*IndexType(Tuple[0] - Begin(0)) +
    Stride(1)*IndexType(Tuple[1] - Begin(1)) +
    Stride(2)*IndexType(Tuple[2] - Begin(2)) +
    Stride(3)*IndexType(Tuple[3] - Begin(3));
}
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,5> &Begin,
  const elem<IndexType,5> &Stride, const ArrayOrIterType &Tuple) {
  return
    Stride(0)*IndexType(Tuple[0] - Begin(0)) +
    Stride(1)*IndexType(Tuple[1] - Begin(1)) +
    Stride(2)*IndexType(Tuple[2] - Begin(2)) +
    Stride(3)*IndexType(Tuple[3] - Begin(3)) +
    Stride(4)*IndexType(Tuple[4] - Begin(4));
}
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,6> &Begin,
  const elem<IndexType,6> &Stride, const ArrayOrIterType &Tuple) {
  return
    Stride(0)*IndexType(Tuple[0] - Begin(0)) +
    Stride(1)*IndexType(Tuple[1] - Begin(1)) +
    Stride(2)*IndexType(Tuple[2] - Begin(2)) +
    Stride(3)*IndexType(Tuple[3] - Begin(3)) +
    Stride(4)*IndexType(Tuple[4] - Begin(4)) +
    Stride(5)*IndexType(Tuple[5] - Begin(5));
}
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,7> &Begin,
  const elem<IndexType,7> &Stride, const ArrayOrIterType &Tuple) {
  return
    Stride(0)*IndexType(Tuple[0] - Begin(0)) +
    Stride(1)*IndexType(Tuple[1] - Begin(1)) +
    Stride(2)*IndexType(Tuple[2] - Begin(2)) +
    Stride(3)*IndexType(Tuple[3] - Begin(3)) +
    Stride(4)*IndexType(Tuple[4] - Begin(4)) +
    Stride(5)*IndexType(Tuple[5] - Begin(5)) +
    Stride(6)*IndexType(Tuple[6] - Begin(6));
}
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,8> &Begin,
  const elem<IndexType,8> &Stride, const ArrayOrIterType &Tuple) {
  return
    Stride(0)*IndexType(Tuple[0] - Begin(0)) +
    Stride(1)*IndexType(Tuple[1] - Begin(1)) +
    Stride(2)*IndexType(Tuple[2] - Begin(2)) +
    Stride(3)*IndexType(Tuple[3] - Begin(3)) +
    Stride(4)*IndexType(Tuple[4] - Begin(4)) +
    Stride(5)*IndexType(Tuple[5] - Begin(5)) +
    Stride(6)*IndexType(Tuple[6] - Begin(6)) +
    Stride(7)*IndexType(Tuple[7] - Begin(7));
}
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,9> &Begin,
  const elem<IndexType,9> &Stride, const ArrayOrIterType &Tuple) {
  return
    Stride(0)*IndexType(Tuple[0] - Begin(0)) +
    Stride(1)*IndexType(Tuple[1] - Begin(1)) +
    Stride(2)*IndexType(Tuple[2] - Begin(2)) +
    Stride(3)*IndexType(Tuple[3] - Begin(3)) +
    Stride(4)*IndexType(Tuple[4] - Begin(4)) +
    Stride(5)*IndexType(Tuple[5] - Begin(5)) +
    Stride(6)*IndexType(Tuple[6] - Begin(6)) +
    Stride(7)*IndexType(Tuple[7] - Begin(7)) +
    Stride(8)*IndexType(Tuple[8] - Begin(8));
}
template <typename IndexType, typename TupleElementType, typename ArrayOrIterType> constexpr
  OVK_FORCE_INLINE IndexType TupleToIndexHelper(const elem<TupleElementType,10> &Begin,
  const elem<IndexType,10> &Stride, const ArrayOrIterType &Tuple) {
  return
    Stride(0)*IndexType(Tuple[0] - Begin(0)) +
    Stride(1)*IndexType(Tuple[1] - Begin(1)) +
    Stride(2)*IndexType(Tuple[2] - Begin(2)) +
    Stride(3)*IndexType(Tuple[3] - Begin(3)) +
    Stride(4)*IndexType(Tuple[4] - Begin(4)) +
    Stride(5)*IndexType(Tuple[5] - Begin(5)) +
    Stride(6)*IndexType(Tuple[6] - Begin(6)) +
    Stride(7)*IndexType(Tuple[7] - Begin(7)) +
    Stride(8)*IndexType(Tuple[8] - Begin(8)) +
    Stride(9)*IndexType(Tuple[9] - Begin(9));
}

template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,1> &Begin, const elem<IndexType,1> &,
  TupleElementType Element0) {
  return IndexType(Element0 - Begin(0));
}
template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,2> &Begin, const elem<IndexType,2> &Stride,
  TupleElementType Element0, TupleElementType Element1) {
  return
    Stride(0)*IndexType(Element0 - Begin(0)) +
    Stride(1)*IndexType(Element1 - Begin(1));
}
template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,3> &Begin, const elem<IndexType,3> &Stride,
  TupleElementType Element0, TupleElementType Element1, TupleElementType Element2) {
  return
    Stride(0)*IndexType(Element0 - Begin(0)) +
    Stride(1)*IndexType(Element1 - Begin(1)) +
    Stride(2)*IndexType(Element2 - Begin(2));
}
template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,4> &Begin, const elem<IndexType,4> &Stride,
  TupleElementType Element0, TupleElementType Element1, TupleElementType Element2, TupleElementType
  Element3) {
  return
    Stride(0)*IndexType(Element0 - Begin(0)) +
    Stride(1)*IndexType(Element1 - Begin(1)) +
    Stride(2)*IndexType(Element2 - Begin(2)) +
    Stride(3)*IndexType(Element3 - Begin(3));
}
template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,5> &Begin, const elem<IndexType,5> &Stride,
  TupleElementType Element0, TupleElementType Element1, TupleElementType Element2, TupleElementType
  Element3, TupleElementType Element4) {
  return
    Stride(0)*IndexType(Element0 - Begin(0)) +
    Stride(1)*IndexType(Element1 - Begin(1)) +
    Stride(2)*IndexType(Element2 - Begin(2)) +
    Stride(3)*IndexType(Element3 - Begin(3)) +
    Stride(4)*IndexType(Element4 - Begin(4));
}
template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,6> &Begin, const elem<IndexType,6> &Stride,
  TupleElementType Element0, TupleElementType Element1, TupleElementType Element2, TupleElementType
  Element3, TupleElementType Element4, TupleElementType Element5) {
  return
    Stride(0)*IndexType(Element0 - Begin(0)) +
    Stride(1)*IndexType(Element1 - Begin(1)) +
    Stride(2)*IndexType(Element2 - Begin(2)) +
    Stride(3)*IndexType(Element3 - Begin(3)) +
    Stride(4)*IndexType(Element4 - Begin(4)) +
    Stride(5)*IndexType(Element5 - Begin(5));
}
template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,7> &Begin, const elem<IndexType,7> &Stride,
  TupleElementType Element0, TupleElementType Element1, TupleElementType Element2, TupleElementType
  Element3, TupleElementType Element4, TupleElementType Element5, TupleElementType Element6) {
  return
    Stride(0)*IndexType(Element0 - Begin(0)) +
    Stride(1)*IndexType(Element1 - Begin(1)) +
    Stride(2)*IndexType(Element2 - Begin(2)) +
    Stride(3)*IndexType(Element3 - Begin(3)) +
    Stride(4)*IndexType(Element4 - Begin(4)) +
    Stride(5)*IndexType(Element5 - Begin(5)) +
    Stride(6)*IndexType(Element6 - Begin(6));
}
template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,8> &Begin, const elem<IndexType,8> &Stride,
  TupleElementType Element0, TupleElementType Element1, TupleElementType Element2, TupleElementType
  Element3, TupleElementType Element4, TupleElementType Element5, TupleElementType Element6,
  TupleElementType Element7) {
  return
    Stride(0)*IndexType(Element0 - Begin(0)) +
    Stride(1)*IndexType(Element1 - Begin(1)) +
    Stride(2)*IndexType(Element2 - Begin(2)) +
    Stride(3)*IndexType(Element3 - Begin(3)) +
    Stride(4)*IndexType(Element4 - Begin(4)) +
    Stride(5)*IndexType(Element5 - Begin(5)) +
    Stride(6)*IndexType(Element6 - Begin(6)) +
    Stride(7)*IndexType(Element7 - Begin(7));
}
template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,9> &Begin, const elem<IndexType,9> &Stride,
  TupleElementType Element0, TupleElementType Element1, TupleElementType Element2, TupleElementType
  Element3, TupleElementType Element4, TupleElementType Element5, TupleElementType Element6,
  TupleElementType Element7, TupleElementType Element8) {
  return
    Stride(0)*IndexType(Element0 - Begin(0)) +
    Stride(1)*IndexType(Element1 - Begin(1)) +
    Stride(2)*IndexType(Element2 - Begin(2)) +
    Stride(3)*IndexType(Element3 - Begin(3)) +
    Stride(4)*IndexType(Element4 - Begin(4)) +
    Stride(5)*IndexType(Element5 - Begin(5)) +
    Stride(6)*IndexType(Element6 - Begin(6)) +
    Stride(7)*IndexType(Element7 - Begin(7)) +
    Stride(8)*IndexType(Element8 - Begin(8));
}
template <typename IndexType, typename TupleElementType> constexpr OVK_FORCE_INLINE IndexType
  TupleElementsToIndexHelper(const elem<TupleElementType,10> &Begin, const elem<IndexType,10> &Stride,
  TupleElementType Element0, TupleElementType Element1, TupleElementType Element2, TupleElementType
  Element3, TupleElementType Element4, TupleElementType Element5, TupleElementType Element6,
  TupleElementType Element7, TupleElementType Element8, TupleElementType Element9) {
  return
    Stride(0)*IndexType(Element0 - Begin(0)) +
    Stride(1)*IndexType(Element1 - Begin(1)) +
    Stride(2)*IndexType(Element2 - Begin(2)) +
    Stride(3)*IndexType(Element3 - Begin(3)) +
    Stride(4)*IndexType(Element4 - Begin(4)) +
    Stride(5)*IndexType(Element5 - Begin(5)) +
    Stride(6)*IndexType(Element6 - Begin(6)) +
    Stride(7)*IndexType(Element7 - Begin(7)) +
    Stride(8)*IndexType(Element8 - Begin(8)) +
    Stride(9)*IndexType(Element9 - Begin(9));
}

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
    std::is_convertible<core::iterator_reference_type<IterType>, tuple_element_type>::value)>
    constexpr OVK_FORCE_INLINE index_type ToIndex(IterType First) const {
    return TupleToIndexHelper<index_type, tuple_element_type>(Begin_, Stride_, First);
  }

  // Want to use array directly here instead of constructing intermediate elem type
  template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
    !core::IsIterator<ArrayType>() && std::is_convertible<core::array_access_type<const ArrayType
    &>, tuple_element_type>::value && core::ArrayRank<ArrayType>() == 1 &&
    (core::ArrayHasRuntimeExtents<ArrayType>() || (core::StaticArrayHasBegin<ArrayType,0>() &&
    core::StaticArrayHasEnd<ArrayType,Rank>())))> constexpr OVK_FORCE_INLINE index_type ToIndex(
    const ArrayType &Array) const {
    return TupleToIndexHelper<index_type, tuple_element_type>(Begin_, Stride_, Array);
  }

  constexpr OVK_FORCE_INLINE index_type ToIndex(TupleElementTypes... TupleElements) const {
    return TupleElementsToIndexHelper<index_type, tuple_element_type>(Begin_, Stride_,
      TupleElements...);
  }

  tuple_type ToTuple(index_type Index) const {
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

  tuple_type IndexToTuple_(layout_tag<array_layout::ROW_MAJOR>, index_type Index) const {

    tuple_type Tuple;

    index_type ReducedIndex = Index;
    for (int i = 0; i < Rank; ++i) {
      tuple_element_type Offset = tuple_element_type(ReducedIndex/Stride_(i));
      Tuple(i) = Begin_(i) + Offset;
      ReducedIndex -= Offset*Stride_(i);
    }

    return Tuple;

  }

  tuple_type IndexToTuple_(layout_tag<array_layout::COLUMN_MAJOR>, index_type Index) const {

    tuple_type Tuple;

    index_type ReducedIndex = Index;
    for (int i = Rank-1; i >= 0; --i) {
      tuple_element_type Offset = tuple_element_type(ReducedIndex/Stride_(i));
      Tuple(i) = Begin_(i) + Offset;
      ReducedIndex -= Offset*Stride_(i);
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
  bool operator==(const indexer<IndexType, TupleElementType, Rank, Layout> &Left, const indexer<
  IndexType, TupleElementType, Rank, Layout> &Right) {
  return Left.Begin() == Right.Begin() && Left.Stride() == Right.Stride();
}

template <typename IndexType, typename TupleElementType, int Rank, array_layout Layout> constexpr
  bool operator!=(const indexer<IndexType, TupleElementType, Rank, Layout> &Left, const indexer<
  IndexType, TupleElementType, Rank, Layout> &Right) {
  return !(Left == Right);
}

}

#endif
