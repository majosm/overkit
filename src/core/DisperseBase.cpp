// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/DisperseBase.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/ConnectivityD.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"

namespace ovk {
namespace core {

template <array_layout Layout> void disperse_base<Layout>::Initialize(const exchange &Exchange,
  int Count, const range &GridValuesRange) {

  Count_ = Count;
  GridValuesRange_ = GridValuesRange;
  GridValuesIndexer_ = range_indexer(GridValuesRange);

  const connectivity &Connectivity = *Exchange.Connectivity_;

  GetConnectivityReceiverSide(Connectivity, Receivers_);
  const connectivity_r &Receivers = *Receivers_;

  GetConnectivityReceiverSideCount(Receivers, NumReceivers_);

  GetConnectivityReceiverSideGrid(Receivers, Grid_);
  const grid &Grid = *Grid_;

  if (OVK_DEBUG) {
    const range &LocalRange = Grid.LocalRange();
    OVK_DEBUG_ASSERT(GridValuesRange.Includes(LocalRange), "Invalid grid values range.");
  }

  Points_ = Receivers.Points_;

}

template class disperse_base<array_layout::ROW_MAJOR>;
template class disperse_base<array_layout::COLUMN_MAJOR>;

template <typename T, array_layout Layout> void disperse_base_for_type<T, Layout>::Initialize(const
  exchange &Exchange, int Count, const range &GridValuesRange) {

  parent_type::Initialize(Exchange, Count, GridValuesRange);

  ReceiverValues_.Resize({Count_});
  GridValues_.Resize({Count_});

}

template <typename T, array_layout Layout> void disperse_base_for_type<T, Layout>::SetBufferViews(
  const void * const *ReceiverValuesVoid, void **GridValuesVoid) {

  OVK_DEBUG_ASSERT(ReceiverValuesVoid || Count_ == 0, "Invalid receiver values pointer.");
  OVK_DEBUG_ASSERT(GridValuesVoid || Count_ == 0, "Invalid grid values pointer.");

  for (int iCount = 0; iCount < Count_; ++iCount) {
    OVK_DEBUG_ASSERT(ReceiverValuesVoid[iCount] || NumReceivers_ == 0, "Invalid receiver values "
      "pointer.");
    ReceiverValues_(iCount) = {static_cast<const value_type *>(ReceiverValuesVoid[iCount]),
      {NumReceivers_}};
  }

  for (int iCount = 0; iCount < Count_; ++iCount) {
    OVK_DEBUG_ASSERT(GridValuesVoid[iCount] || NumReceivers_ == 0, "Invalid grid values pointer.");
    GridValues_(iCount) = {static_cast<value_type *>(GridValuesVoid[iCount]),
      {GridValuesRange_.Count()}};
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
