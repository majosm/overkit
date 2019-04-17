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

template <typename T, array_layout Layout> class disperse_base_for_type : public disperse_base<
  Layout> {

private:

  using parent_type = disperse_base<Layout>;

public:

  using value_type = T;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange);

protected:

  using typename parent_type::range_indexer;

  using parent_type::Receivers_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::NumReceivers_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesIndexer_;
  using parent_type::Points_;

  array<array_view<const value_type>> ReceiverValues_;
  array<array_view<value_type>> GridValues_;

  void SetBufferViews(const void * const *ReceiverValuesVoid, void **GridValuesVoid);

};

extern template class disperse_base_for_type<bool, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<bool, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<unsigned char, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<unsigned char, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<int, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<int, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<long, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<long, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<long long, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<long long, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<unsigned int, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<unsigned int, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<unsigned long, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<unsigned long, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<unsigned long long, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<unsigned long long, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<float, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<float, array_layout::COLUMN_MAJOR>;
extern template class disperse_base_for_type<double, array_layout::ROW_MAJOR>;
extern template class disperse_base_for_type<double, array_layout::COLUMN_MAJOR>;

}}

#endif
