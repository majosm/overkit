// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISPERSE_BASE_HPP_INCLUDED
#define OVK_CORE_DISPERSE_BASE_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/ConnectivityD.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Exchange.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

template <array_layout Layout> class disperse_base {

public:

  disperse_base();
  disperse_base(const disperse_base &Other) = delete;
  disperse_base(disperse_base &&Other) noexcept;

  disperse_base &operator=(const disperse_base &Other);
  disperse_base &operator=(disperse_base &&Other) noexcept;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange);

protected:

  using range_indexer = indexer<long long, int, MAX_DIMS, Layout>;

  const connectivity_r *Receivers_;
  const grid *Grid_;
  int Count_;
  long long NumReceivers_;
  range GridValuesRange_;
  range_indexer GridValuesIndexer_;
  array_view<const int,2> Points_;

};

extern template class disperse_base<array_layout::ROW_MAJOR>;
extern template class disperse_base<array_layout::COLUMN_MAJOR>;

}}

#endif
