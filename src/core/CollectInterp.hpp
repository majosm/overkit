// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_INTERP_HPP_INCLUDED
#define OVK_CORE_COLLECT_INTERP_HPP_INCLUDED

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

#include <memory>
#include <utility>

namespace ovk {
namespace core {
namespace collect_internal {

template <typename T, array_layout Layout> class collect_interp : public collect_base_for_type<T,
  Layout> {

protected:

  using parent_type = collect_base_for_type<T, Layout>;

  using typename parent_type::cell_indexer;
  using parent_type::Context_;
  using parent_type::CollectMap_;
  using parent_type::Count_;
  using parent_type::FieldValues_;
  using parent_type::PackedValues_;
  using parent_type::REDUCE_TIME;

public:

  using typename parent_type::value_type;

  collect_interp(std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart, const range
    &LocalRange, const collect_map &CollectMap, int Count, const range &FieldValuesRange,
    floating_ref<const array<double,3>> InterpCoefs):
    parent_type(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange),
    InterpCoefs_(InterpCoefs)
  {

    parent_type::AllocateRemoteValues_(RemoteValues_);

    VertexValues_.Resize({{Count_,CollectMap_->MaxVertices()}});
    VertexCoefs_.Resize({CollectMap_->MaxVertices()});

  }

  collect_interp(const collect_interp &Other) = delete;
  collect_interp(collect_interp &&Other) noexcept = default;

  collect_interp &operator=(const collect_interp &Other) = delete;
  collect_interp &operator=(collect_interp &&Other) noexcept = default;

  void Collect(const void *FieldValuesVoid, void *PackedValuesVoid) {

    profiler &Profiler = Context_->core_Profiler();

    const array<double,3> &InterpCoefs = *InterpCoefs_;

    parent_type::SetBufferViews_(FieldValuesVoid, PackedValuesVoid);
    parent_type::RetrieveRemoteValues_(FieldValues_, RemoteValues_);

    Profiler.Start(REDUCE_TIME);

    for (long long iCell = 0; iCell < CollectMap_->Count(); ++iCell) {

      range CellRange = parent_type::GetCellRange_(iCell);
      cell_indexer CellIndexer(CellRange);
      int NumVertices = CellRange.Count<int>();

      parent_type::AssembleVertexValues_(FieldValues_, RemoteValues_, iCell, CellRange, CellIndexer,
        VertexValues_);

      for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
        for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
          for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
            int iVertex = CellIndexer.ToIndex(i,j,k);
            VertexCoefs_(iVertex) =
              InterpCoefs(0,i-CellRange.Begin(0),iCell) *
              InterpCoefs(1,j-CellRange.Begin(1),iCell) *
              InterpCoefs(2,k-CellRange.Begin(2),iCell);
          }
        }
      }

      for (int iCount = 0; iCount < Count_; ++iCount) {
        PackedValues_(iCount)(iCell) = value_type(0);
        for (int iVertex = 0; iVertex < NumVertices; ++iVertex) {
          PackedValues_(iCount)(iCell) += VertexCoefs_(iVertex)*VertexValues_(iCount,iVertex);
        }
      }

    }

    Profiler.Stop(REDUCE_TIME);

  }

private:

  floating_ref<const array<double,3>> InterpCoefs_;
  array<array<value_type,2>> RemoteValues_;
  array<value_type,2> VertexValues_;
  array<double> VertexCoefs_;

};

}}}

#endif
