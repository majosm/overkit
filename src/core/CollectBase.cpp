// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/CollectBase.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/CollectMap.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {
namespace collect_internal {

template <array_layout Layout> collect_base<Layout>::collect_base(std::shared_ptr<context>
  &&Context, comm_view Comm, const cart &Cart, const range &LocalRange, const collect_map
  &CollectMap, int Count, const range &FieldValuesRange, int NumThreads):
  Context_(std::move(Context)),
  Comm_(Comm),
  Cart_(Cart),
  LocalRange_(LocalRange),
  CollectMap_(CollectMap.GetFloatingRef()),
  Count_(Count),
  FieldValuesRange_(FieldValuesRange),
  FieldValuesIndexer_(FieldValuesRange)
{

  const array<collect_map::send> &Sends = CollectMap_->Sends();
  const array<collect_map::recv> &Recvs = CollectMap_->Recvs();

  Requests_.Reserve(Sends.Count()+Recvs.Count());

  LocalVertexCellIndices_.Resize({NumThreads});
  LocalVertexFieldValuesIndices_.Resize({NumThreads});
  for (int iThread = 0; iThread < NumThreads; ++iThread) {
    LocalVertexCellIndices_(iThread).Resize({CollectMap_->MaxVertices()});
    LocalVertexFieldValuesIndices_(iThread).Resize({CollectMap_->MaxVertices()});
  }

}

template <array_layout Layout> range collect_base<Layout>::GetCellRange_(long long iCell) const {

  const array_view<const int,3> &CellExtents = CollectMap_->CellExtents();

  range CellRange;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    CellRange.Begin(iDim) = CellExtents(0,iDim,iCell);
    CellRange.End(iDim) = CellExtents(1,iDim,iCell);
  }

  return CellRange;

}

template <array_layout Layout> void collect_base<Layout>::GetLocalCellInfo_(const range &CellRange,
  const range_indexer<int,Layout> &CellIndexer, int &NumLocalVertices, array_view<int>
  LocalVertexCellIndices, array_view<long long> LocalVertexFieldValuesIndices) const {

  bool AwayFromEdge = Cart_.Range().Includes(CellRange);

  if (AwayFromEdge) {
    range LocalCellRange = IntersectRanges(LocalRange_, CellRange);
    int iLocalVertex = 0;
    for (int k = LocalCellRange.Begin(2); k < LocalCellRange.End(2); ++k) {
      for (int j = LocalCellRange.Begin(1); j < LocalCellRange.End(1); ++j) {
        for (int i = LocalCellRange.Begin(0); i < LocalCellRange.End(0); ++i) {
          LocalVertexCellIndices(iLocalVertex) = CellIndexer.ToIndex(i,j,k);
          LocalVertexFieldValuesIndices(iLocalVertex) = FieldValuesIndexer_.ToIndex(i,j,k);
          ++iLocalVertex;
        }
      }
    }
    NumLocalVertices = iLocalVertex;
  } else {
    int iLocalVertex = 0;
    for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
      for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
        for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
          tuple<int> Vertex = {i,j,k};
          tuple<int> AdjustedVertex = Cart_.PeriodicAdjust(Vertex);
          if (LocalRange_.Contains(AdjustedVertex)) {
            LocalVertexCellIndices(iLocalVertex) = CellIndexer.ToIndex(Vertex);
            LocalVertexFieldValuesIndices(iLocalVertex) = FieldValuesIndexer_.ToIndex(AdjustedVertex);
            ++iLocalVertex;
          }
        }
      }
    }
    NumLocalVertices = iLocalVertex;
  }

}

template class collect_base<array_layout::ROW_MAJOR>;
template class collect_base<array_layout::COLUMN_MAJOR>;

template <typename T, array_layout Layout> collect_base_for_type<T, Layout>::collect_base_for_type(
  std::shared_ptr<context> &&Context, comm_view Comm, const cart &Cart, const range &LocalRange,
  const collect_map &CollectMap, int Count, const range &FieldValuesRange, int NumThreads):
  parent_type(std::move(Context), Comm, Cart, LocalRange, CollectMap, Count, FieldValuesRange,
    NumThreads)
{

  const array<collect_map::send> &Sends = CollectMap_->Sends();
  const array<collect_map::recv> &Recvs = CollectMap_->Recvs();

  SendBuffers_.Resize({Sends.Count()});
  for (int iSend = 0; iSend < Sends.Count(); ++iSend) {
    SendBuffers_(iSend).Resize({{Count_,Sends(iSend).NumPoints}});
  }

  if (!std::is_same<value_type, mpi_value_type>::value) {
    RecvBuffers_.Resize({Recvs.Count()});
    for (int iRecv = 0; iRecv < Recvs.Count(); ++iRecv) {
      RecvBuffers_(iRecv).Resize({{Count_,Recvs(iRecv).NumPoints}});
    }
  }

  FieldValues_.Resize({Count_});
  PackedValues_.Resize({Count_});

}

template <typename T, array_layout Layout> void collect_base_for_type<T, Layout>::
  AllocateRemoteValues_(array<array<value_type,2>> &RemoteValues) const {

  const array<collect_map::recv> &Recvs = CollectMap_->Recvs();

  RemoteValues.Resize({Recvs.Count()});
  for (int iRecv = 0; iRecv < Recvs.Count(); ++iRecv) {
    RemoteValues(iRecv).Resize({{Count_,Recvs(iRecv).NumPoints}});
  }

}

template <typename T, array_layout Layout> void collect_base_for_type<T, Layout>::
  SetBufferViews_(const void *FieldValuesVoid, void *PackedValuesVoid) {

  auto FieldValuesRaw = static_cast<const value_type * const *>(FieldValuesVoid);
  auto PackedValuesRaw = static_cast<value_type **>(PackedValuesVoid);

  OVK_DEBUG_ASSERT(FieldValuesRaw || Count_ == 0, "Invalid field values pointer.");
  OVK_DEBUG_ASSERT(PackedValuesRaw || Count_ == 0, "Invalid packed values pointer.");

  long long NumCells = CollectMap_->Count();

  for (int iCount = 0; iCount < Count_; ++iCount) {
    OVK_DEBUG_ASSERT(FieldValuesRaw[iCount] || NumCells == 0, "Invalid field values pointer.");
    FieldValues_(iCount) = {FieldValuesRaw[iCount], {FieldValuesRange_.Count()}};
  }

  for (int iCount = 0; iCount < Count_; ++iCount) {
    OVK_DEBUG_ASSERT(PackedValuesRaw[iCount] || NumCells == 0, "Invalid packed values pointer.");
    PackedValues_(iCount) = {PackedValuesRaw[iCount], {NumCells}};
  }

}

template <typename T, array_layout Layout> void collect_base_for_type<T, Layout>::
  RetrieveRemoteValues_(array_view<array_view<const value_type>> FieldValues, array<array<
  value_type,2>> &RemoteValues) {

  core::profiler &Profiler = Context_->core_Profiler();

  MPI_Datatype MPIDataType = core::GetMPIDataType<mpi_value_type>();

  const array<collect_map::send> &Sends = CollectMap_->Sends();
  const array<collect_map::recv> &Recvs = CollectMap_->Recvs();

  Profiler.Start(MPI_TIME);

  if (std::is_same<value_type, mpi_value_type>::value) {
    for (int iRecv = 0; iRecv < Recvs.Count(); ++iRecv) {
      const collect_map::recv &Recv = Recvs(iRecv);
      MPI_Irecv(RemoteValues(iRecv).Data(), Count_*Recv.NumPoints, MPIDataType, Recv.Rank, 0,
        Comm_, &Requests_.Append());
    }
  } else {
    for (int iRecv = 0; iRecv < Recvs.Count(); ++iRecv) {
      const collect_map::recv &Recv = Recvs(iRecv);
      MPI_Irecv(RecvBuffers_(iRecv).Data(), Count_*Recv.NumPoints, MPIDataType, Recv.Rank, 0,
        Comm_, &Requests_.Append());
    }
  }

  Profiler.Stop(MPI_TIME);
  Profiler.Start(PACK_TIME);

  for (int iSend = 0; iSend < Sends.Count(); ++iSend) {
    const collect_map::send &Send = Sends(iSend);
    for (long long iSendPoint = 0; iSendPoint < Send.NumPoints; ++iSendPoint) {
      tuple<int> Point = {
        Send.Points(0,iSendPoint),
        Send.Points(1,iSendPoint),
        Send.Points(2,iSendPoint)
      };
      long long iFieldPoint = FieldValuesIndexer_.ToIndex(Point);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        SendBuffers_(iSend)(iCount,iSendPoint) = mpi_value_type(FieldValues(iCount)(iFieldPoint));
      }
    }
  }

  Profiler.Stop(PACK_TIME);
  Profiler.Start(MPI_TIME);

  for (int iSend = 0; iSend < Sends.Count(); ++iSend) {
    const collect_map::send &Send = Sends(iSend);
    MPI_Isend(SendBuffers_(iSend).Data(), Count_*Send.NumPoints, MPIDataType, Send.Rank, 0, Comm_,
      &Requests_.Append());
  }

  MPI_Waitall(Requests_.Count(), Requests_.Data(), MPI_STATUSES_IGNORE);
  Requests_.Clear();

  Profiler.Stop(MPI_TIME);

  if (!std::is_same<value_type, mpi_value_type>::value) {
    for (int iRecv = 0; iRecv < RecvBuffers_.Count(); ++iRecv) {
      for (long long iBuffer = 0; iBuffer < RecvBuffers_(iRecv).Count(); ++iBuffer) {
        RemoteValues(iRecv)[iBuffer] = value_type(RecvBuffers_(iRecv)[iBuffer]);
      }
    }
  }

}

template <typename T, array_layout Layout> void collect_base_for_type<T, Layout>::
  AssembleVertexValues_(array_view<array_view<const value_type>> FieldValues, const array<
  array<value_type,2>> &RemoteValues, long long iCell, const range &CellRange, const range_indexer<
  int,Layout> &CellIndexer, array_view<value_type,2> VertexValues, int iThread) {

  array<int> &LocalVertexCellIndices = LocalVertexCellIndices_(iThread);
  array<long long> &LocalVertexFieldValuesIndices = LocalVertexFieldValuesIndices_(iThread);

  int NumLocalVertices;
  parent_type::GetLocalCellInfo_(CellRange, CellIndexer, NumLocalVertices, LocalVertexCellIndices,
    LocalVertexFieldValuesIndices);

  // Fill in the local data
  for (int iCount = 0; iCount < Count_; ++iCount) {
    for (int iLocalVertex = 0; iLocalVertex < NumLocalVertices; ++iLocalVertex) {
      VertexValues(iCount,LocalVertexCellIndices(iLocalVertex)) = FieldValues(iCount)(
        LocalVertexFieldValuesIndices(iLocalVertex));
    }
  }

  // Fill in the remote data
  const array<int> &NumRemoteVertices = CollectMap_->RemoteVertexCounts();
  const array<long long *> &RemoteVertices = CollectMap_->RemoteVertices();
  const array<int *> &RemoteVertexRecvs = CollectMap_->RemoteVertexRecvs();
  const array<long long *> &RemoteVertexRecvBufferIndices =
    CollectMap_->RemoteVertexRecvBufferIndices();

  for (int iRemoteVertex = 0; iRemoteVertex < NumRemoteVertices(iCell); ++iRemoteVertex) {
    int iVertex = RemoteVertices(iCell)[iRemoteVertex];
    int iRecv = RemoteVertexRecvs(iCell)[iRemoteVertex];
    long long iValue = RemoteVertexRecvBufferIndices(iCell)[iRemoteVertex];
    for (int iCount = 0; iCount < Count_; ++iCount) {
      VertexValues(iCount,iVertex) = RemoteValues(iRecv)(iCount,iValue);
    }
  }

}

template class collect_base_for_type<bool, array_layout::ROW_MAJOR>;
template class collect_base_for_type<bool, array_layout::COLUMN_MAJOR>;
template class collect_base_for_type<unsigned char, array_layout::ROW_MAJOR>;
template class collect_base_for_type<unsigned char, array_layout::COLUMN_MAJOR>;
template class collect_base_for_type<int, array_layout::ROW_MAJOR>;
template class collect_base_for_type<int, array_layout::COLUMN_MAJOR>;
template class collect_base_for_type<long, array_layout::ROW_MAJOR>;
template class collect_base_for_type<long, array_layout::COLUMN_MAJOR>;
template class collect_base_for_type<long long, array_layout::ROW_MAJOR>;
template class collect_base_for_type<long long, array_layout::COLUMN_MAJOR>;
template class collect_base_for_type<unsigned int, array_layout::ROW_MAJOR>;
template class collect_base_for_type<unsigned int, array_layout::COLUMN_MAJOR>;
template class collect_base_for_type<unsigned long, array_layout::ROW_MAJOR>;
template class collect_base_for_type<unsigned long, array_layout::COLUMN_MAJOR>;
template class collect_base_for_type<unsigned long long, array_layout::ROW_MAJOR>;
template class collect_base_for_type<unsigned long long, array_layout::COLUMN_MAJOR>;
template class collect_base_for_type<float, array_layout::ROW_MAJOR>;
template class collect_base_for_type<float, array_layout::COLUMN_MAJOR>;
template class collect_base_for_type<double, array_layout::ROW_MAJOR>;
template class collect_base_for_type<double, array_layout::COLUMN_MAJOR>;

}}}
