// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ARRAY_OPS_HPP_INCLUDED
#define OVK_CORE_ARRAY_OPS_HPP_INCLUDED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/ForEach.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/ScalarTraits.hpp>

#include <type_traits>
#include <utility>

namespace ovk {

template <typename ArrayRefType, typename FRef, OVK_FUNCTION_REQUIRES(core::IsArray<core::
  remove_cvref<ArrayRefType>>() && core::IsCallableWith<FRef &&, core::array_access_type<
  ArrayRefType &&>>())> void ArrayForEach(ArrayRefType &&Array, FRef &&Func) {

  long long NumValues = core::ArrayCount(Array);

  for (long long i = 0; i < NumValues; ++i) {
    std::forward<FRef>(Func)(static_cast<core::array_access_type<ArrayRefType &&>>(Array[i]));
  }

}

template <typename T, int Rank, array_layout Layout, typename FRef, OVK_FUNCTION_REQUIRES(
  core::IsCallableWith<FRef &&, T &>())> void ArrayForEach(const array_view<T, Rank, Layout> &View,
  FRef &&Func) {

  for (auto &Value : View) {
    std::forward<FRef>(Func)(Value);
  }

}

template <typename ArrayRefType, typename TransformedArrayType, typename FRef,
  OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<ArrayRefType>>() && core::IsArray<
  TransformedArrayType>() && core::ArraysAreCongruent<core::remove_cvref<ArrayRefType>,
  TransformedArrayType>() && core::IsCallableAs<FRef &&, core::array_value_type<
  TransformedArrayType>(core::array_access_type<ArrayRefType &&>)>())> void ArrayTransform(
  ArrayRefType &&Array, TransformedArrayType &TransformedArray, FRef &&Func) {

  long long NumValues = core::ArrayCount(Array);

  for (long long i = 0; i < NumValues; ++i) {
    TransformedArray[i] = std::forward<FRef>(Func)(static_cast<core::array_access_type<ArrayRefType
      &&>>(Array[i]));
  }

}

template <typename ArrayRefType, typename T, int Rank, array_layout Layout, typename FRef,
  OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<ArrayRefType>>() && !std::is_const<T>::
  value && core::ArrayHasFootprint<core::remove_cvref<ArrayRefType>, Rank, Layout>() &&
  core::IsCallableAs<FRef &&, T(core::array_access_type<ArrayRefType &&>)>())> void ArrayTransform(
  ArrayRefType &&Array, const array_view<T, Rank, Layout> &TransformedView, FRef &&Func) {

  long long NumValues = core::ArrayCount(Array);

  for (long long i = 0; i < NumValues; ++i) {
    TransformedView[i] = std::forward<FRef>(Func)(static_cast<core::array_access_type<ArrayRefType
      &&>>(Array[i]));
  }

}

template <typename T, int Rank, array_layout Layout, typename TransformedArrayType, typename FRef,
  OVK_FUNCTION_REQUIRES(core::IsArray<TransformedArrayType>() && core::ArrayHasFootprint<
  TransformedArrayType, Rank, Layout>() && core::IsCallableAs<FRef &&, core::array_value_type<
  TransformedArrayType>(T &)>())> void ArrayTransform(const array_view<T, Rank, Layout> &View,
  TransformedArrayType &TransformedArray, FRef &&Func) {

  long long NumValues = View.Count();

  for (long long i = 0; i < NumValues; ++i) {
    TransformedArray[i] = std::forward<FRef>(Func)(View[i]);
  }

}

template <typename T, typename U, int Rank, array_layout Layout, typename FRef,
  OVK_FUNCTION_REQUIRES(!std::is_const<U>::value && core::IsCallableAs<FRef &&, U(T &)>())> void
  ArrayTransform(const array_view<T, Rank, Layout> &View, const array_view<U, Rank, Layout>
  &TransformedView, FRef &&Func) {

  long long NumValues = View.Count();

  for (long long i = 0; i < NumValues; ++i) {
    TransformedView[i] = std::forward<FRef>(Func)(View[i]);
  }

}

// Not sure if I like these overloads (awkward due to different possible index and tuple element
// types; revisit later. May want to consider taking separate arguments instead of (or in addition
// to) tuple, like this:
//   ovk::ArrayForEach(InterpCoefs, [](int &Coef, int iDim, int iPoint, long long iDonor) {
//     Coef = ...;
//   });
// (though this still has index type conversion issues).

// template <typename IndexType=long long, typename ArrayRefType, typename FRef, OVK_FUNCTION_REQUIRES(
//   core::IsArray<core::remove_cvref<ArrayRefType>>() && !core::IsCallableWith<FRef &&,
//   core::array_access_type<ArrayRefType &&>>() && core::IsCallableWith<FRef &&,
//   core::array_access_type<ArrayRefType &&>, IndexType>())> void ArrayForEach(ArrayRefType &&Array,
//   FRef &&Func) {

//   long long NumValues = core::ArrayCount(Array);

//   for (long long i = 0; i < NumValues; ++i) {
//     std::forward<FRef>(Func)(static_cast<core::array_access_type<ArrayRefType &&>>(Array[i]),
//       IndexType(i));
//   }

// }

// template <typename IndexType=long long, typename T, int Rank, array_layout Layout, typename FRef,
//   OVK_FUNCTION_REQUIRES(!core::IsCallableWith<FRef &&, T &>() && core::IsCallableWith<FRef &&,
//   T &, IndexType>())> void ArrayForEach(const array_view<T, Rank, Layout> &View, FRef &&Func) {

//   for (long long i = 0; i < View.Count(); ++i) {
//     std::forward<FRef>(Func)(View[i], IndexType(i));
//   }

// }

// template <typename TupleElementType=long long, typename ArrayRefType, typename FRef,
//   OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<ArrayRefType>>() && !core::IsCallableWith<
//   FRef &&, core::array_access_type<ArrayRefType &&>>() && !core::IsCallableWith<FRef &&,
//   core::array_access_type<ArrayRefType &&>, long long>() && core::IsCallableWith<FRef &&,
//   core::array_access_type<ArrayRefType &&>, const elem<TupleElementType,core::ArrayRank<
//   core::remove_cvref<ArrayRefType>>()> &>())> void ArrayForEach(ArrayRefType &&Array, FRef &&Func) {

//   using array_type = core::remove_cvref<ArrayRefType>;
//   constexpr int Rank = core::ArrayRank<array_type>();
//   constexpr array_layout Layout = core::ArrayLayout<array_type>();

//   interval<TupleElementType,Rank> Extents = core::ArrayExtents<TupleElementType>(Array);

//   indexer<long long, TupleElementType, Rank, Layout> Indexer(Extents);

//   core::ForEach<Layout>(Extents, [&](const elem<TupleElementType,Rank> &Tuple) {
//     long long i = Indexer.ToIndex(Tuple);
//     std::forward<FRef>(Func)(static_cast<core::array_access_type<ArrayRefType &&>>(Array[i]),
//       Tuple);
//   });

// }

// template <typename TupleElementType=long long, typename T, int Rank, array_layout Layout, typename
//   FRef, OVK_FUNCTION_REQUIRES(!core::IsCallableWith<FRef &&, T &>() && !core::IsCallableWith<FRef
//   &&, T &, long long>() && core::IsCallableWith<FRef &&, T &, const elem<TupleElementType,Rank>
//   &>())> void ArrayForEach(const array_view<T, Rank, Layout> &View, FRef &&Func) {

//   // View extents value type might not be the same as TupleElementType
//   interval<TupleElementType,Rank> Extents(View.Extents());

//   // Can't use View.Indexer() because its tuple element type may not be the same as TupleElementType
//   indexer<long long, TupleElementType, Rank, Layout> Indexer(Extents);

//   core::ForEach<Layout>(Extents, [&](const elem<TupleElementType,Rank> &Tuple) {
//     long long i = Indexer.ToIndex(Tuple);
//     std::forward<FRef>(Func)(View[i], Tuple);
//   });

// }

template <typename ArrayRefType, typename U, typename FRef, OVK_FUNCTION_REQUIRES(core::IsArray<
  core::remove_cvref<ArrayRefType>>() && core::IsCallableWith<FRef &&, U &, core::array_access_type<
  ArrayRefType &&>>())> U ArrayReduce(ArrayRefType &&Array, U BaseValue, FRef &&Func) {

  U Result = std::move(BaseValue);

  long long NumValues = core::ArrayCount(Array);

  for (long long i = 0; i < NumValues; ++i) {
    std::forward<FRef>(Func)(Result, static_cast<core::array_access_type<ArrayRefType &&>>(
      Array[i]));
  }

  return Result;

}

template <typename T, int Rank, array_layout Layout, typename U, typename FRef,
  OVK_FUNCTION_REQUIRES(core::IsCallableWith<FRef &&, U &, T &>())> U ArrayReduce(const
  array_view<T, Rank, Layout> &View, U BaseValue, FRef &&Func) {

  U Result = std::move(BaseValue);

  for (auto &Value : View) {
    std::forward<FRef>(Func)(Result, Value);
  }

  return Result;

}

// See note above

// template <typename IndexType=long long, typename ArrayRefType, typename U, typename FRef,
//   OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<ArrayRefType>>() &&
//   !core::IsCallableWith<FRef &&, U &, core::array_access_type<ArrayRefType &&>>() &&
//   core::IsCallableWith<FRef &&, U &, core::array_access_type<ArrayRefType &&>, IndexType>())> U
//   ArrayReduce(ArrayRefType &&Array, U BaseValue, FRef &&Func) {

//   U Result = std::move(BaseValue);

//   long long NumValues = core::ArrayCount(Array);

//   for (long long i = 0; i < NumValues; ++i) {
//     std::forward<FRef>(Func)(Result, static_cast<core::array_access_type<ArrayRefType &&>>(
//       Array[i]), IndexType(i));
//   }

//   return Result;

// }

// template <typename IndexType=long long, typename T, int Rank, array_layout Layout, typename U,
//   typename FRef, OVK_FUNCTION_REQUIRES(!core::IsCallableWith<FRef &&, U &, T &>() &&
//   core::IsCallableWith<FRef &&, U &, T &, IndexType>())> void ArrayReduce(const array_view<T, Rank,
//   Layout> &View, U BaseValue, FRef &&Func) {

//   U Result = std::move(BaseValue);

//   for (long long i = 0; i < View.Count(); ++i) {
//     std::forward<FRef>(Func)(Result, View[i], IndexType(i));
//   }

//   return Result;

// }

// template <typename TupleElementType=long long, typename ArrayRefType, typename U, typename FRef,
//   OVK_FUNCTION_REQUIRES(core::IsArray<core::remove_cvref<ArrayRefType>>() && !core::IsCallableWith<
//   FRef &&, U &, core::array_access_type<ArrayRefType &&>>() && !core::IsCallableWith<FRef &&, U &,
//   core::array_access_type<ArrayRefType &&>, long long>() && core::IsCallableWith<FRef &&, U &,
//   core::array_access_type<ArrayRefType &&>, const elem<TupleElementType,core::ArrayRank<
//   core::remove_cvref<ArrayRefType>>()> &>())> void ArrayReduce(ArrayRefType &&Array, U BaseValue,
//   FRef &&Func) {

//   using array_type = core::remove_cvref<ArrayRefType>;
//   constexpr int Rank = core::ArrayRank<array_type>();
//   constexpr array_layout Layout = core::ArrayLayout<array_type>();

//   interval<TupleElementType,Rank> Extents = core::ArrayExtents<TupleElementType>(Array);

//   indexer<long long, TupleElementType, Rank, Layout> Indexer(Extents);

//   U Result = std::move(BaseValue);

//   core::ForEach<Layout>(Extents, [&](const elem<TupleElementType,Rank> &Tuple) {
//     long long i = Indexer.ToIndex(Tuple);
//     std::forward<FRef>(Func)(Result, static_cast<core::array_access_type<ArrayRefType &&>>(
//       Array[i]), Tuple);
//   });

//   return Result;

// }

// template <typename TupleElementType=long long, typename T, int Rank, array_layout Layout, typename
//   U, typename FRef, OVK_FUNCTION_REQUIRES(!core::IsCallableWith<FRef &&, U &, T &>() &&
//   !core::IsCallableWith<FRef &&, U &, T &, long long>() && core::IsCallableWith<FRef &&, U &, T &,
//   const elem<TupleElementType,Rank> &>())> void ArrayReduce(const array_view<T, Rank, Layout>
//   &View, U BaseValue, FRef &&Func) {

//   // View extents value type might not be the same as TupleElementType
//   interval<TupleElementType,Rank> Extents(View.Extents());

//   // Can't use View.Indexer() because its tuple element type may not be the same as TupleElementType
//   indexer<long long, TupleElementType, Rank, Layout> Indexer(Extents);

//   U Result = std::move(BaseValue);

//   core::ForEach<Layout>(Extents, [&](const elem<TupleElementType,Rank> &Tuple) {
//     long long i = Indexer.ToIndex(Tuple);
//     std::forward<FRef>(Func)(Result, View[i], Tuple);
//   });

//   return Result;

// }

template <typename ArrayRefType, typename FRef, OVK_FUNCTION_REQUIRES(core::IsArray<
  core::remove_cvref<ArrayRefType>>() && core::IsCallableWith<FRef &&, core::array_value_type<
  core::remove_cvref<ArrayRefType>> &, core::array_access_type<ArrayRefType &&>>())>
  core::array_value_type<core::remove_cvref<ArrayRefType>> ArrayCollapse(ArrayRefType &&Array, FRef
  &&Func) {

  using value_type = core::array_value_type<core::remove_cvref<ArrayRefType>>;

  long long NumValues = core::ArrayCount(Array);

  OVK_DEBUG_ASSERT(NumValues > 0, "Cannot collapse zero-size array.");

  value_type Result = static_cast<core::array_access_type<ArrayRefType &&>>(Array[0]);

  for (long long i = 1; i < NumValues; ++i) {
    std::forward<FRef>(Func)(Result, static_cast<core::array_access_type<ArrayRefType &&>>(
      Array[i]));
  }

  return Result;

}

template <typename T, int Rank, array_layout Layout, typename FRef, OVK_FUNCTION_REQUIRES(
  core::IsCallableWith<FRef &&, typename std::remove_const<T>::type &, T &>())> typename
  std::remove_const<T>::type ArrayCollapse(const array_view<T, Rank, Layout> &View, FRef &&Func) {

  using value_type = typename std::remove_const<T>::type;

  OVK_DEBUG_ASSERT(View.Count() > 0, "Cannot collapse zero-size array.");

  value_type Result = View[0];

  for (int i = 1; i < View.Count(); ++i) {
    std::forward<FRef>(Func)(Result, View[i]);
  }

  return Result;

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>())> void ArrayFill(
  ArrayType &Array, const core::array_value_type<ArrayType> &Value_) {

  using value_type = core::array_value_type<ArrayType>;

  ArrayForEach(Array, [&](value_type &Value) { Value = Value_; });

}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(!std::is_const<T>::value
  )> void ArrayFill(const array_view<T, Rank, Layout> &View, const T &Value) {

  View.Fill(Value);

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>())> void ArrayFill(
  ArrayType &Array, std::initializer_list<core::array_value_type<ArrayType>> ValuesList) {

  using value_type = core::array_value_type<ArrayType>;

  auto Iter = ValuesList.begin();
  ArrayForEach(Array, [&](value_type &Value) { Value = *Iter++; });

}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(!std::is_const<T>::value
  )> void ArrayFill(const array_view<T, Rank, Layout> &View, std::initializer_list<T> ValuesList) {

  View.Fill(ValuesList);

}

template <typename ArrayType, typename IterType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  core::IsInputIterator<IterType>() && std::is_convertible<core::iterator_reference_type<IterType>,
  core::array_value_type<ArrayType>>::value)> void ArrayFill(ArrayType &Array, IterType First) {

  using value_type = core::array_value_type<ArrayType>;

  IterType Iter = First;
  ArrayForEach(Array, [&](value_type &Value) { Value = *Iter++; });

}

template <typename T, int Rank, array_layout Layout, typename IterType, OVK_FUNCTION_REQUIRES(
  !std::is_const<T>::value && core::IsInputIterator<IterType>() && std::is_convertible<
  core::iterator_reference_type<IterType>, T>::value)> void ArrayFill(const array_view<T, Rank,
  Layout> &View, IterType First) {

  View.Fill(First);

}

template <typename ArrayType, typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(
  core::IsArray<ArrayType>() && core::ArrayHasFootprint<ArrayType, Rank, Layout>() &&
  std::is_convertible<typename std::remove_const<T>::type, core::array_value_type<ArrayType>>::
  value)> void ArrayFill(ArrayType &Array, const array_view<T, Rank, Layout> &SourceView) {

  using value_type = core::array_value_type<ArrayType>;

  long long i = 0;
  ArrayForEach(Array, [&](value_type &Value) { Value = SourceView[i++]; });

}

template <typename T, typename U, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(
  !std::is_const<T>::value && std::is_convertible<typename std::remove_const<U>::type, T>::value)>
  void ArrayFill(const array_view<T, Rank, Layout> &View, const array_view<U, Rank, Layout>
  &SourceView) {

  View.Fill(SourceView);

}

template <typename ArrayType, typename SourceArrayRefType, OVK_FUNCTION_REQUIRES(core::IsArray<
  ArrayType>() && core::IsArray<core::remove_cvref<SourceArrayRefType>>() && !core::IsIterator<
  typename std::decay<SourceArrayRefType>::type>() && core::ArraysAreCongruent<ArrayType,
  core::remove_cvref<SourceArrayRefType>>() && std::is_convertible<core::array_access_type<
  SourceArrayRefType &&>, core::array_value_type<ArrayType>>::value)> void ArrayFill(ArrayType
  &Array, SourceArrayRefType &&SourceArray) {

  using value_type = core::array_value_type<ArrayType>;

  long long i = 0;
  ArrayForEach(Array, [&](value_type &Value) {
    Value = static_cast<core::array_access_type<SourceArrayRefType &&>>(SourceArray[i++]);
  });

}

template <typename T, int Rank, array_layout Layout, typename SourceArrayRefType,
  OVK_FUNCTION_REQUIRES(!std::is_const<T>::value && core::IsArray<core::remove_cvref<
  SourceArrayRefType>>() && !core::IsIterator<typename std::decay<SourceArrayRefType>::type>()
  && core::ArrayHasFootprint<core::remove_cvref<SourceArrayRefType>, Rank, Layout>() &&
  std::is_convertible<core::array_access_type<SourceArrayRefType &&>, T>::value)> void ArrayFill(
  const array_view<T, Rank, Layout> &View, SourceArrayRefType &&SourceArray) {

  View.Fill(std::forward<SourceArrayRefType>(SourceArray));

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  std::is_convertible<core::array_access_type<ArrayType>, bool>::value)> bool ArrayNone(const
  ArrayType &Array) {

  using value_type = core::array_value_type<ArrayType>;

  return ArrayReduce(Array, true, [](bool &Partial, const value_type &Value) {
    Partial = Partial && !static_cast<bool>(Value);
  });

}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(std::is_convertible<T,
  bool>::value)> bool ArrayNone(const array_view<T, Rank, Layout> &View) {

  return ArrayReduce(View, true, [](bool &Partial, T &Value) {
    Partial = Partial && !static_cast<bool>(Value);
  });

}

template <typename ArrayType, typename F, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  core::IsCallableAs<F, bool(core::array_access_type<const ArrayType &>)>())> bool ArrayNone(const
  ArrayType &Array, F Condition) {

  using value_type = core::array_value_type<ArrayType>;

  return ArrayReduce(Array, true, [&](bool &Partial, const value_type &Value) {
    Partial = Partial && !Condition(Value);
  });

}

template <typename T, int Rank, array_layout Layout, typename F, OVK_FUNCTION_REQUIRES(
  core::IsCallableAs<F, bool(T &)>())> bool ArrayNone(const array_view<T, Rank, Layout> &View, F
  Condition) {

  return ArrayReduce(View, true, [&](bool &Partial, T &Value) {
    Partial = Partial && !Condition(Value);
  });

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  std::is_convertible<core::array_access_type<const ArrayType &>, bool>::value)> bool ArrayAny(const
  ArrayType &Array) {

  using value_type = core::array_value_type<ArrayType>;

  return ArrayReduce(Array, false, [](bool &Partial, const value_type &Value) {
    Partial = Partial || static_cast<bool>(Value);
  });

}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(std::is_convertible<T,
  bool>::value)> bool ArrayAny(const array_view<T, Rank, Layout> &View) {

  return ArrayReduce(View, false, [](bool &Partial, T &Value) {
    Partial = Partial || static_cast<bool>(Value);
  });

}

template <typename ArrayType, typename F, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  core::IsCallableAs<F, bool(core::array_access_type<const ArrayType &>)>())> bool ArrayAny(const
  ArrayType &Array, F Condition) {

  using value_type = core::array_value_type<ArrayType>;

  return ArrayReduce(Array, false, [&](bool &Partial, const value_type &Value) {
    Partial = Partial || Condition(Value);
  });

}

template <typename T, int Rank, array_layout Layout, typename F, OVK_FUNCTION_REQUIRES(
  core::IsCallableAs<F, bool(T &)>())> bool ArrayAny(const array_view<T, Rank, Layout> &View, F
  Condition) {

  return ArrayReduce(View, false, [&](bool &Partial, T &Value) {
    Partial = Partial || Condition(Value);
  });

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  std::is_convertible<core::array_access_type<const ArrayType &>, bool>::value)> bool ArrayNotAll(
  const ArrayType &Array) {

  using value_type = core::array_value_type<ArrayType>;

  return ArrayReduce(Array, false, [](bool &Partial, const value_type &Value) {
    Partial = Partial || !static_cast<bool>(Value);
  });

}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(std::is_convertible<T,
  bool>::value)> bool ArrayNotAll(const array_view<T, Rank, Layout> &View) {

  return ArrayReduce(View, false, [](bool &Partial, T &Value) {
    Partial = Partial || !static_cast<bool>(Value);
  });

}

template <typename ArrayType, typename F, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  core::IsCallableAs<F, bool(core::array_access_type<const ArrayType &>)>())> bool ArrayNotAll(const
  ArrayType &Array, F Condition) {

  using value_type = core::array_value_type<ArrayType>;

  return ArrayReduce(Array, false, [&](bool &Partial, const value_type &Value) {
    Partial = Partial || !Condition(Value);
  });

}

template <typename T, int Rank, array_layout Layout, typename F, OVK_FUNCTION_REQUIRES(
  core::IsCallableAs<F, bool(T &)>())> bool ArrayNotAll(const array_view<T, Rank, Layout> &View, F
  Condition) {

  return ArrayReduce(View, false, [&](bool &Partial, T &Value) {
    Partial = Partial || !Condition(Value);
  });

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  std::is_convertible<core::array_access_type<const ArrayType &>, bool>::value)> bool ArrayAll(const
  ArrayType &Array) {

  using value_type = core::array_value_type<ArrayType>;

  return ArrayReduce(Array, true, [](bool &Partial, const value_type &Value) {
    Partial = Partial && static_cast<bool>(Value);
  });

}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(std::is_convertible<T,
  bool>::value)> bool ArrayAll(const array_view<T, Rank, Layout> &View) {

  return ArrayReduce(View, true, [](bool &Partial, T &Value) {
    Partial = Partial && static_cast<bool>(Value);
  });

}

template <typename ArrayType, typename F, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() &&
  core::IsCallableAs<F, bool(core::array_access_type<const ArrayType &>)>())> bool ArrayAll(const
  ArrayType &Array, F Condition) {

  using value_type = core::array_value_type<ArrayType>;

  return ArrayReduce(Array, true, [&](bool &Partial, const value_type &Value) {
    Partial = Partial && Condition(Value);
  });

}

template <typename T, int Rank, array_layout Layout, typename F, OVK_FUNCTION_REQUIRES(
  core::IsCallableAs<F, bool(T &)>())> bool ArrayAll(const array_view<T, Rank, Layout> &View, F
  Condition) {

  return ArrayReduce(View, true, [&](bool &Partial, T &Value) {
    Partial = Partial && Condition(Value);
  });

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() && core::IsScalar<
  core::array_value_type<ArrayType>>())> core::array_value_type<ArrayType> ArrayMin(const ArrayType
  &Array) {

  OVK_DEBUG_ASSERT(core::ArrayCount(Array) > 0, "Cannot take min of zero-size array.");

  using value_type = core::array_value_type<ArrayType>;

  return ArrayCollapse(Array, [](value_type &Partial, const value_type &Value) {
    Partial = Min(Partial, Value);
  });

}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(core::IsScalar<typename
  std::remove_const<T>::type>())> typename std::remove_const<T>::type ArrayMin(const array_view<T,
  Rank, Layout> &View) {

  OVK_DEBUG_ASSERT(View.Count() > 0, "Cannot take min of zero-size array.");

  using value_type = typename std::remove_const<T>::type;

  return ArrayCollapse(View, [](value_type &Partial, const value_type &Value) {
    Partial = Min(Partial, Value);
  });

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>() && core::IsScalar<
  core::array_value_type<ArrayType>>())> core::array_value_type<ArrayType> ArrayMax(const ArrayType
  &Array) {

  OVK_DEBUG_ASSERT(core::ArrayCount(Array) > 0, "Cannot take max of zero-size array.");

  using value_type = core::array_value_type<ArrayType>;

  return ArrayCollapse(Array, [](value_type &Partial, const value_type &Value) {
    Partial = Max(Partial, Value);
  });

}

template <typename T, int Rank, array_layout Layout, OVK_FUNCTION_REQUIRES(core::IsScalar<typename
  std::remove_const<T>::type>())> typename std::remove_const<T>::type ArrayMax(const array_view<T,
  Rank, Layout> &View) {

  OVK_DEBUG_ASSERT(View.Count() > 0, "Cannot take max of zero-size array.");

  using value_type = typename std::remove_const<T>::type;

  return ArrayCollapse(View, [](value_type &Partial, const value_type &Value) {
    Partial = Max(Partial, Value);
  });

}

template <typename ArrayType, OVK_FUNCTION_REQUIRES(core::IsArray<ArrayType>())>
  core::array_value_type<ArrayType> ArraySum(const ArrayType &Array) {

  OVK_DEBUG_ASSERT(core::ArrayCount(Array) > 0, "Cannot sum zero-size array.");

  using value_type = core::array_value_type<ArrayType>;

  return ArrayCollapse(Array, [](value_type &Partial, const value_type &Value) {
    Partial += Value;
  });


}

template <typename T, int Rank, array_layout Layout> T ArraySum(const array_view<T, Rank, Layout>
  &View) {

  OVK_DEBUG_ASSERT(View.Count() > 0, "Cannot sum zero-size array.");

  using value_type = typename std::remove_const<T>::type;

  return ArrayCollapse(View, [](value_type &Partial, const value_type &Value) {
    Partial += Value;
  });

}

}

#endif
