// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_BASE_HPP_INCLUDED
#define OVK_CORE_COLLECT_BASE_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Indexer.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

#include <type_traits>
#include <utility>

namespace ovk {
namespace core {
namespace collect_internal {

// Put as much as possible in non-type-specific base class to reduce compile times
template <array_layout Layout> class collect_base {

protected:

  collect_base(comm_view Comm, const cart &Cart, const range &LocalRange, const collect_map
    &CollectMap, int Count, const range &FieldValuesRange, profiler &Profiler);

  // Can't define these here due to issues with GCC < 6.3 and Intel < 17
  // implementations of extern template
//   collect_base(const collect_base &Other) = delete;
//   collect_base(collect_base &&Other) noexcept = default;

//   collect_base &operator=(const collect_base &Other) = delete;
//   collect_base &operator=(collect_base &&Other) noexcept = default;

  using field_indexer = indexer<long long, int, MAX_DIMS, Layout>;
  using cell_indexer = indexer<int, int, MAX_DIMS, Layout>;

  comm_view Comm_;
  cart Cart_;
  range LocalRange_;
  const collect_map *CollectMap_;
  mutable profiler *Profiler_;
  int Count_;
  int MaxPointsInCell_;
  range FieldValuesRange_;
  field_indexer FieldValuesIndexer_;
  array<MPI_Request> Requests_;
  array<int> LocalVertexCellIndices_;
  array<long long> LocalVertexFieldValuesIndices_;

  range GetCellRange_(long long iCell) const;

  void GetLocalCellInfo_(const range &CellRange, const cell_indexer &CellIndexer, int
    &NumLocalVertices, array_view<int> LocalCellIndices, array_view<long long>
    LocalFieldValuesIndices) const;

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

  collect_base_for_type(comm_view Comm, const cart &Cart, const range &LocalRange, const collect_map
    &CollectMap, int Count, const range &FieldValuesRange, profiler &Profiler);

protected:

  using typename parent_type::field_indexer;
  using typename parent_type::cell_indexer;
  using parent_type::Comm_;
  using parent_type::Cart_;
  using parent_type::LocalRange_;
  using parent_type::CollectMap_;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::FieldValuesRange_;
  using parent_type::FieldValuesIndexer_;

  array<array_view<const value_type>> FieldValues_;
  array<array_view<value_type>> PackedValues_;

  void AllocateRemoteValues_(array<array<value_type,2>> &RemoteValues) const;

  void SetBufferViews_(const void * const *FieldValuesVoid, void **PackedValuesVoid);

  void RetrieveRemoteValues_(array_view<array_view<const value_type>> FieldValues, array<array<
    value_type,2>> &RemoteValues);

  void AssembleVertexValues_(array_view<array_view<const value_type>> FieldValues, const
    array<array<value_type,2>> &RemoteValues, long long iCell, const range &CellRange, const
    cell_indexer &CellIndexer, array_view<value_type,2> VertexValues);

private:

  using parent_type::Requests_;
  using parent_type::LocalVertexCellIndices_;
  using parent_type::LocalVertexFieldValuesIndices_;

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

}}}

#endif
