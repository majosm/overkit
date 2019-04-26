// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_NOT_ALL_HPP_INCLUDED
#define OVK_CORE_COLLECT_NOT_ALL_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/CollectBase.hpp>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

namespace ovk {
namespace core {
namespace collect_internal {

template <typename T, array_layout Layout> class collect_not_all : public collect_base_for_type<T,
  Layout> {

protected:

  using parent_type = collect_base_for_type<T, Layout>;

  using typename parent_type::cell_indexer;
  using parent_type::CollectMap_;
  using parent_type::Profiler_;
  using parent_type::Count_;
  using parent_type::FieldValues_;
  using parent_type::PackedValues_;

public:

  using typename parent_type::value_type;

  collect_not_all(comm_view Comm, const cart &Cart, const range &LocalRange, const collect_map
    &CollectMap, int Count, const range &FieldValuesRange, profiler &Profiler):
    parent_type(Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange, Profiler),
    MemAllocTime_(GetProfilerTimerID(*Profiler_, "Collect::MemAlloc")),
    ReduceTime_(GetProfilerTimerID(*Profiler_, "Collect::Reduce"))
  {

    StartProfile(*Profiler_, MemAllocTime_);

    parent_type::AllocateRemoteValues_(RemoteValues_);

    VertexValues_.Resize({{Count_,CollectMap_->MaxVertices()}});

    EndProfile(*Profiler_, MemAllocTime_);

  }

  collect_not_all(const collect_not_all &Other) = delete;
  collect_not_all(collect_not_all &&Other) noexcept = default;

  collect_not_all &operator=(const collect_not_all &Other) = delete;
  collect_not_all &operator=(collect_not_all &&Other) noexcept = default;

  void Collect(const void * const *FieldValuesVoid, void **PackedValuesVoid) {

    parent_type::SetBufferViews_(FieldValuesVoid, PackedValuesVoid);
    parent_type::RetrieveRemoteValues_(FieldValues_, RemoteValues_);

    StartProfile(*Profiler_, ReduceTime_);

    for (long long iCell = 0; iCell < CollectMap_->Count(); ++iCell) {

      range CellRange = parent_type::GetCellRange_(iCell);
      cell_indexer CellIndexer(CellRange);
      int NumVertices = CellRange.Count<int>();

      parent_type::AssembleVertexValues_(FieldValues_, RemoteValues_, iCell, CellRange, CellIndexer,
        VertexValues_);

      for (int iCount = 0; iCount < Count_; ++iCount) {
        PackedValues_(iCount)(iCell) = value_type(false);
        for (int iVertex = 0; iVertex < NumVertices; ++iVertex) {
          PackedValues_(iCount)(iCell) = PackedValues_(iCount)(iVertex) || !VertexValues_(iCount,
            iVertex);
        }
      }

    }

    EndProfile(*Profiler_, ReduceTime_);

  }

private:

  int MemAllocTime_;
  int ReduceTime_;
  array<array<value_type,2>> RemoteValues_;
  array<value_type,2> VertexValues_;

};

}}}

#endif
