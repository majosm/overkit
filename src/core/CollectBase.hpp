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

template <typename T, array_layout Layout> class collect_base {

  using mpi_value_type = core::mpi_compatible_type<T>;

public:

  using value_type = T;

  collect_base();
  collect_base(const collect_base &Other) = delete;
  collect_base(collect_base &&Other) noexcept;

  collect_base &operator=(const collect_base &Other) = delete;
  collect_base &operator=(collect_base &&Other) noexcept;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange);

protected:

  using range_indexer = indexer<long long, int, MAX_DIMS, Layout>;
  using donor_indexer = indexer<int, int, MAX_DIMS, Layout>;

  const connectivity_d *Donors_;
  const grid *Grid_;
  core::profiler *Profiler_;
  int Count_;
  long long NumDonors_;
  int MaxPointsInCell_;
  range_indexer GridValuesIndexer_;

  void AllocateRemoteDonorValues(array<array<value_type,2>> &RemoteDonorValues);

  void RetrieveRemoteDonorValues(array_view<array_view<const value_type>> GridValues, array<array<
    value_type,2>> &RemoteDonorValues);

  void AssembleDonorPointValues(array_view<array_view<const value_type>> GridValues, const
    array<array<value_type,2>> &RemoteDonorValues, long long iDonor, elem<int,MAX_DIMS>
    &DonorSize, array_view<value_type,2> DonorPointValues);

private:

  cart Cart_;
  range GlobalRange_;
  range LocalRange_;
  array_view<const exchange::collect_send> Sends_;
  array_view<const exchange::collect_recv> Recvs_;
  array<array<mpi_value_type,2>> SendBuffers_;
  array<array<mpi_value_type,2>> RecvBuffers_;
  array<MPI_Request> Requests_;
  array_view<const int> NumRemoteDonorPoints_;
  array_view<const long long * const> RemoteDonorPoints_;
  array_view<const int * const> RemoteDonorPointCollectRecvs_;
  array_view<const long long * const> RemoteDonorPointCollectRecvBufferIndices_;
  array<int> LocalDonorPointIndices_;
  array<long long> LocalDonorPointGridValuesIndices_;

};

extern template class collect_base<bool, array_layout::ROW_MAJOR>;
extern template class collect_base<bool, array_layout::COLUMN_MAJOR>;
extern template class collect_base<unsigned char, array_layout::ROW_MAJOR>;
extern template class collect_base<unsigned char, array_layout::COLUMN_MAJOR>;
extern template class collect_base<int, array_layout::ROW_MAJOR>;
extern template class collect_base<int, array_layout::COLUMN_MAJOR>;
extern template class collect_base<long, array_layout::ROW_MAJOR>;
extern template class collect_base<long, array_layout::COLUMN_MAJOR>;
extern template class collect_base<long long, array_layout::ROW_MAJOR>;
extern template class collect_base<long long, array_layout::COLUMN_MAJOR>;
extern template class collect_base<unsigned int, array_layout::ROW_MAJOR>;
extern template class collect_base<unsigned int, array_layout::COLUMN_MAJOR>;
extern template class collect_base<unsigned long, array_layout::ROW_MAJOR>;
extern template class collect_base<unsigned long, array_layout::COLUMN_MAJOR>;
extern template class collect_base<unsigned long long, array_layout::ROW_MAJOR>;
extern template class collect_base<unsigned long long, array_layout::COLUMN_MAJOR>;
extern template class collect_base<float, array_layout::ROW_MAJOR>;
extern template class collect_base<float, array_layout::COLUMN_MAJOR>;
extern template class collect_base<double, array_layout::ROW_MAJOR>;
extern template class collect_base<double, array_layout::COLUMN_MAJOR>;

}}

#endif
