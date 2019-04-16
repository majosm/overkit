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

template <array_layout Layout> disperse_base<Layout>::disperse_base() = default;
template <array_layout Layout> disperse_base<Layout>::disperse_base(disperse_base &&Other) noexcept
  = default;
template <array_layout Layout> disperse_base<Layout> &disperse_base<Layout>::operator=(disperse_base
  &&Other) noexcept = default;

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
    range LocalRange;
    GetGridLocalRange(Grid, LocalRange);
    OVK_DEBUG_ASSERT(GridValuesRange.Includes(LocalRange), "Invalid grid values range.");
  }

  Points_ = Receivers.Points_;

}

template class disperse_base<array_layout::ROW_MAJOR>;
template class disperse_base<array_layout::COLUMN_MAJOR>;

}}
