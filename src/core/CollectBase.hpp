// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_BASE_HPP_INCLUDED
#define OVK_CORE_COLLECT_BASE_HPP_INCLUDED

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

// Put as much as possible in non-type-specific base class to reduce compile times
template <array_layout Layout> class collect_base {

public:

  collect_base();
  collect_base(const collect_base &Other) = delete;
  collect_base(collect_base &&Other) noexcept;

  collect_base &operator=(const collect_base &Other) = delete;
  collect_base &operator=(collect_base &&Other) noexcept;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange);

  void GetLocalDonorPointInfo_(long long iDonor, elem<int,MAX_DIMS> &DonorSize, int
    &NumLocalDonorPoints, array_view<int> LocalDonorPointIndices, array_view<long long>
    LocalDonorPointGridValuesIndices);

  using range_indexer = indexer<long long, int, MAX_DIMS, Layout>;
  using donor_indexer = indexer<int, int, MAX_DIMS, Layout>;

  const connectivity_d *Donors_;
  const grid *Grid_;
  cart Cart_;
  range GlobalRange_;
  range LocalRange_;
  core::profiler *Profiler_;
  int Count_;
  long long NumDonors_;
  int MaxPointsInCell_;
  range GridValuesRange_;
  range_indexer GridValuesIndexer_;
  array_view<const exchange::collect_send> Sends_;
  array_view<const exchange::collect_recv> Recvs_;
  array<MPI_Request> Requests_;
  array_view<const int> NumRemoteDonorPoints_;
  array_view<const long long * const> RemoteDonorPoints_;
  array_view<const int * const> RemoteDonorPointCollectRecvs_;
  array_view<const long long * const> RemoteDonorPointCollectRecvBufferIndices_;
  array<int> LocalDonorPointIndices_;
  array<long long> LocalDonorPointGridValuesIndices_;

};

extern template class collect_base<array_layout::ROW_MAJOR>;
extern template class collect_base<array_layout::COLUMN_MAJOR>;

template <typename T, array_layout Layout> class collect_base_for_type : private collect_base<
  Layout> {

private:

  using parent_type = collect_base<Layout>;

  using mpi_value_type = core::mpi_compatible_type<T>;

public:

  using value_type = T;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange);

protected:

  using typename parent_type::range_indexer;
  using typename parent_type::donor_indexer;
  using parent_type::Donors_;
  using parent_type::Grid_;
  using parent_type::Cart_;
  using parent_type::GlobalRange_;
  using parent_type::LocalRange_;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::NumDonors_;
  using parent_type::MaxPointsInCell_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesIndexer_;

  array<array_view<const value_type>> GridValues_;
  array<array_view<value_type>> DonorValues_;

  void AllocateRemoteDonorValues(array<array<value_type,2>> &RemoteDonorValues);

  void SetBufferViews(const void * const *GridValuesVoid, void **DonorValuesVoid);

  void RetrieveRemoteDonorValues(array_view<array_view<const value_type>> GridValues, array<array<
    value_type,2>> &RemoteDonorValues);

  void AssembleDonorPointValues(array_view<array_view<const value_type>> GridValues, const
    array<array<value_type,2>> &RemoteDonorValues, long long iDonor, elem<int,MAX_DIMS>
    &DonorSize, array_view<value_type,2> DonorPointValues);

private:

  using parent_type::Sends_;
  using parent_type::Recvs_;
  using parent_type::Requests_;
  using parent_type::NumRemoteDonorPoints_;
  using parent_type::RemoteDonorPoints_;
  using parent_type::RemoteDonorPointCollectRecvs_;
  using parent_type::RemoteDonorPointCollectRecvBufferIndices_;
  using parent_type::LocalDonorPointIndices_;
  using parent_type::LocalDonorPointGridValuesIndices_;

  array<array<mpi_value_type,2>> SendBuffers_;
  array<array<mpi_value_type,2>> RecvBuffers_;

};

extern template class collect_base_for_type<bool, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<bool, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<unsigned char, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<unsigned char, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<int, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<int, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<long, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<long, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<long long, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<long long, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<unsigned int, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<unsigned int, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<unsigned long, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<unsigned long, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<unsigned long long, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<unsigned long long, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<float, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<float, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<double, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<double, array_layout::COLUMN_MAJOR>;

}}

#endif
