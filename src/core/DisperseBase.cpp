// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/DisperseBase.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"

namespace ovk {
namespace core {

template <array_layout Layout> disperse_base<Layout>::disperse_base(const array<int,2> &Points,
  int Count, const range &FieldValuesRange, profiler &Profiler):
  Points_(Points),
  Profiler_(&Profiler),
  Count_(Count),
  FieldValuesRange_(FieldValuesRange),
  FieldValuesIndexer_(FieldValuesRange)
{}

template class disperse_base<array_layout::ROW_MAJOR>;
template class disperse_base<array_layout::COLUMN_MAJOR>;

template <typename T, array_layout Layout> disperse_base_for_type<T, Layout>::
  disperse_base_for_type(const array<int,2> &Points, int Count, const range &FieldValuesRange,
  profiler &Profiler):
  parent_type(Points, Count, FieldValuesRange, Profiler)
{
  PackedValues_.Resize({Count_});
  FieldValues_.Resize({Count_});
}

template <typename T, array_layout Layout> void disperse_base_for_type<T, Layout>::SetBufferViews(
  const void * const *PackedValuesVoid, void **FieldValuesVoid) {

  long long NumPoints = Points_.Size(1);

  OVK_DEBUG_ASSERT(PackedValuesVoid || Count_ == 0, "Invalid packed values pointer.");
  OVK_DEBUG_ASSERT(FieldValuesVoid || Count_ == 0, "Invalid field values pointer.");

  for (int iCount = 0; iCount < Count_; ++iCount) {
    OVK_DEBUG_ASSERT(PackedValuesVoid[iCount] || NumPoints == 0, "Invalid packed values "
      "pointer.");
    PackedValues_(iCount) = {static_cast<const value_type *>(PackedValuesVoid[iCount]),
      {NumPoints}};
  }

  for (int iCount = 0; iCount < Count_; ++iCount) {
    OVK_DEBUG_ASSERT(FieldValuesVoid[iCount] || NumPoints == 0, "Invalid field values pointer.");
    FieldValues_(iCount) = {static_cast<value_type *>(FieldValuesVoid[iCount]),
      {FieldValuesRange_.Count()}};
  }

}

template class disperse_base_for_type<bool, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<bool, array_layout::COLUMN_MAJOR>;
template class disperse_base_for_type<unsigned char, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<unsigned char, array_layout::COLUMN_MAJOR>;
template class disperse_base_for_type<int, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<int, array_layout::COLUMN_MAJOR>;
template class disperse_base_for_type<long, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<long, array_layout::COLUMN_MAJOR>;
template class disperse_base_for_type<long long, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<long long, array_layout::COLUMN_MAJOR>;
template class disperse_base_for_type<unsigned int, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<unsigned int, array_layout::COLUMN_MAJOR>;
template class disperse_base_for_type<unsigned long, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<unsigned long, array_layout::COLUMN_MAJOR>;
template class disperse_base_for_type<unsigned long long, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<unsigned long long, array_layout::COLUMN_MAJOR>;
template class disperse_base_for_type<float, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<float, array_layout::COLUMN_MAJOR>;
template class disperse_base_for_type<double, array_layout::ROW_MAJOR>;
template class disperse_base_for_type<double, array_layout::COLUMN_MAJOR>;

}}
