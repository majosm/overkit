// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_BASE_HPP_INCLUDED
#define OVK_CORE_COLLECT_BASE_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {
namespace collect_internal {

// Put as much as possible in non-type-specific base class to reduce compile times
template <array_layout Layout> class collect_base {

protected:

  collect_base(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart, const range
    &LocalRange, const collect_map &CollectMap, int Count, const range &FieldValuesRange, int
    NumThreads);

  // Can't define these here due to issues with GCC < 6.3 and Intel < 17
  // implementations of extern template
//   collect_base(const collect_base &Other) = delete;
//   collect_base(collect_base &&Other) noexcept = default;

//   collect_base &operator=(const collect_base &Other) = delete;
//   collect_base &operator=(collect_base &&Other) noexcept = default;

  std::shared_ptr<context> Context_;

  comm_view Comm_;

  cart Cart_;
  range LocalRange_;

  floating_ref<const collect_map> CollectMap_;

  int Count_;
  int MaxPointsInCell_;

  range FieldValuesRange_;
  range_indexer<long long,Layout> FieldValuesIndexer_;

  array<MPI_Request> Requests_;
  array<array<int>> LocalVertexCellIndices_;
  array<array<long long>> LocalVertexFieldValuesIndices_;

  range GetCellRange_(long long iCell) const;

  void GetLocalCellInfo_(const range &CellRange, const range_indexer<int,Layout> &CellIndexer, int
    &NumLocalVertices, array_view<int> LocalCellIndices, array_view<long long>
    LocalFieldValuesIndices) const;

  static constexpr int PACK_TIME = profiler::EXCHANGER_COLLECT_PACK_TIME;
  static constexpr int MPI_TIME = profiler::EXCHANGER_COLLECT_MPI_TIME;
  static constexpr int REDUCE_TIME = profiler::EXCHANGER_COLLECT_REDUCE_TIME;

};

extern template class collect_base<array_layout::ROW_MAJOR>;
extern template class collect_base<array_layout::COLUMN_MAJOR>;

template <typename T, array_layout Layout> class collect_base_for_type : public collect_base<
  Layout> {

private:

  using parent_type = collect_base<Layout>;

  using mpi_value_type = mpi_compatible_type<T>;

public:

  using value_type = T;

  collect_base_for_type(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart, const
    range &LocalRange, const collect_map &CollectMap, int Count, const range &FieldValuesRange, int
    NumThreads=1);

protected:

  using parent_type::Context_;
  using parent_type::Comm_;
  using parent_type::Cart_;
  using parent_type::LocalRange_;
  using parent_type::CollectMap_;
  using parent_type::Count_;
  using parent_type::FieldValuesRange_;
  using parent_type::FieldValuesIndexer_;
  using parent_type::PACK_TIME;
  using parent_type::MPI_TIME;
  using parent_type::REDUCE_TIME;

  array<array_view<const value_type>> FieldValues_;
  array<array_view<value_type>> PackedValues_;

  void AllocateRemoteValues_(array<array<value_type,2>> &RemoteValues) const;

  void SetBufferViews_(const void *FieldValuesVoid, void *PackedValuesVoid);

  void RetrieveRemoteValues_(array_view<array_view<const value_type>> FieldValues, array<array<
    value_type,2>> &RemoteValues);

  void AssembleVertexValues_(array_view<array_view<const value_type>> FieldValues, const
    array<array<value_type,2>> &RemoteValues, long long iCell, const range &CellRange, const
    range_indexer<int,Layout> &CellIndexer, array_view<value_type,2> VertexValues, int iThread=0);

private:

  using parent_type::Requests_;
  using parent_type::LocalVertexCellIndices_;
  using parent_type::LocalVertexFieldValuesIndices_;

  array<array<mpi_value_type,2>> SendBuffers_;
  array<array<mpi_value_type,2>> RecvBuffers_;

};

extern template class collect_base_for_type<bool, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<bool, array_layout::COLUMN_MAJOR>;
extern template class collect_base_for_type<byte, array_layout::ROW_MAJOR>;
extern template class collect_base_for_type<byte, array_layout::COLUMN_MAJOR>;
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

}}}

#endif
