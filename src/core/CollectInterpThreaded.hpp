// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_INTERP_THREADED_HPP_INCLUDED
#define OVK_CORE_COLLECT_INTERP_THREADED_HPP_INCLUDED

#ifdef OVK_HAVE_OPENMP

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/CollectBase.hpp>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>
#include <omp.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {
namespace collect_internal {

template <typename T, array_layout Layout> class collect_interp_threaded : public
  collect_base_for_type<T, Layout> {

protected:

  using parent_type = collect_base_for_type<T, Layout>;

  using parent_type::Context_;
  using parent_type::CollectMap_;
  using parent_type::Count_;
  using parent_type::FieldValues_;
  using parent_type::PackedValues_;
  using parent_type::REDUCE_TIME;

public:

  using typename parent_type::value_type;

  collect_interp_threaded(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
    const range &LocalRange, const collect_map &CollectMap, int Count, const range
    &FieldValuesRange, floating_ref<const array<double,3>> InterpCoefs):
    collect_interp_threaded(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count,
      FieldValuesRange, InterpCoefs, GetThreadCount_())
  {}

  collect_interp_threaded(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart,
    const range &LocalRange, const collect_map &CollectMap, int Count, const range
    &FieldValuesRange, floating_ref<const array<double,3>> InterpCoefs, int NumThreads):
    parent_type(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
      NumThreads),
    InterpCoefs_(InterpCoefs)
  {

    parent_type::AllocateRemoteValues_(RemoteValues_);

    VertexValues_.Resize({NumThreads});
    VertexCoefs_.Resize({NumThreads});
    for (int iThread = 0; iThread < NumThreads; ++iThread) {
      VertexValues_(iThread).Resize({{Count_,CollectMap_->MaxVertices()}});
      VertexCoefs_(iThread).Resize({CollectMap_->MaxVertices()});
    }

  }

  collect_interp_threaded(const collect_interp_threaded &Other) = delete;
  collect_interp_threaded(collect_interp_threaded &&Other) noexcept = default;

  collect_interp_threaded &operator=(const collect_interp_threaded &Other) = delete;
  collect_interp_threaded &operator=(collect_interp_threaded &&Other) noexcept = default;

  void Collect(const void *FieldValuesVoid, void *PackedValuesVoid) {

    profiler &Profiler = Context_->core_Profiler();

    const array<double,3> &InterpCoefs = *InterpCoefs_;

    parent_type::SetBufferViews_(FieldValuesVoid, PackedValuesVoid);
    parent_type::RetrieveRemoteValues_(FieldValues_, RemoteValues_);

    Profiler.Start(REDUCE_TIME);

    long long NumCells = CollectMap_->Count();

    #pragma omp parallel firstprivate(NumCells)
    {

      int iThread = omp_get_thread_num();

      array<value_type,2> &VertexValues = VertexValues_(iThread);
      array<double> &VertexCoefs = VertexCoefs_(iThread);

      #pragma omp for
      for (long long iCell = 0; iCell < NumCells; ++iCell) {

        range CellRange = parent_type::GetCellRange_(iCell);
        range_indexer<int,Layout> CellIndexer(CellRange);
        int NumVertices = CellRange.Count<int>();

        parent_type::AssembleVertexValues_(FieldValues_, RemoteValues_, iCell, CellRange,
          CellIndexer, VertexValues, iThread);

        for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
          for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
            for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
              int iVertex = CellIndexer.ToIndex(i,j,k);
              VertexCoefs(iVertex) =
                InterpCoefs(0,i-CellRange.Begin(0),iCell) *
                InterpCoefs(1,j-CellRange.Begin(1),iCell) *
                InterpCoefs(2,k-CellRange.Begin(2),iCell);
            }
          }
        }

        for (int iCount = 0; iCount < Count_; ++iCount) {
          PackedValues_(iCount)(iCell) = value_type(0);
          for (int iVertex = 0; iVertex < NumVertices; ++iVertex) {
            PackedValues_(iCount)(iCell) += VertexCoefs(iVertex)*VertexValues(iCount,iVertex);
          }
        }

      }

    }

    Profiler.Stop(REDUCE_TIME);

  }

private:

  floating_ref<const array<double,3>> InterpCoefs_;
  array<array<value_type,2>> RemoteValues_;
  array<array<value_type,2>> VertexValues_;
  array<array<double>> VertexCoefs_;

  static int GetThreadCount_() {
    int NumThreads;
#pragma omp parallel
    {
      NumThreads = omp_get_num_threads();
    }
    return NumThreads;
  }

};

}}}

#endif

#endif
