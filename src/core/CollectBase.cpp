// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/CollectBase.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/ConnectivityD.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"

namespace ovk {
namespace core {

template <array_layout Layout> collect_base<Layout>::collect_base() = default;
template <array_layout Layout> collect_base<Layout>::collect_base(collect_base &&Other) noexcept
  = default;
template <array_layout Layout> collect_base<Layout> &collect_base<Layout>::operator=(collect_base
  &&Other) noexcept = default;

template <array_layout Layout> void collect_base<Layout>::Initialize(const exchange &Exchange, int
  Count, const range &GridValuesRange) {

  int NumDims = Exchange.NumDims_;
  Count_ = Count;

  const connectivity &Connectivity = *Exchange.Connectivity_;

  GetConnectivityDonorSide(Connectivity, Donors_);
  const connectivity_d &Donors = *Donors_;

  GetConnectivityDonorSideGrid(Donors, Grid_);
  const grid &Grid = *Grid_;

  Profiler_ = Exchange.Profiler_;
  int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");

  GetGridCart(Grid, Cart_);
  GetGridGlobalRange(Grid, GlobalRange_);
  GetGridLocalRange(Grid, LocalRange_);

  OVK_DEBUG_ASSERT(GridValuesRange.Includes(LocalRange_), "Invalid grid values range.");

  int MaxSize;
  GetConnectivityDonorSideCount(Donors, NumDonors_);
  GetConnectivityDonorSideMaxSize(Donors, MaxSize);

  MaxPointsInCell_ = 1;
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    MaxPointsInCell_ *= MaxSize;
  }

  GridValuesIndexer_ = range_indexer(GridValuesRange);

  Sends_ = Exchange.CollectSends_;
  Recvs_ = Exchange.CollectRecvs_;

  int NumSends = Sends_.Count();
  int NumRecvs = Recvs_.Count();

  core::StartProfile(*Profiler_, MemAllocTime);

  Requests_.Reserve(NumSends+NumRecvs);

  core::EndProfile(*Profiler_, MemAllocTime);

  NumRemoteDonorPoints_ = Exchange.NumRemoteDonorPoints_;
  RemoteDonorPoints_ = Exchange.RemoteDonorPoints_;
  RemoteDonorPointCollectRecvs_ = Exchange.RemoteDonorPointCollectRecvs_;
  RemoteDonorPointCollectRecvBufferIndices_ = Exchange.RemoteDonorPointCollectRecvBufferIndices_;

  core::StartProfile(*Profiler_, MemAllocTime);

  LocalDonorPointIndices_.Resize({MaxPointsInCell_});
  LocalDonorPointGridValuesIndices_.Resize({MaxPointsInCell_});

  core::EndProfile(*Profiler_, MemAllocTime);

}

template <array_layout Layout> void collect_base<Layout>::GetLocalDonorPointInfo_(long long iDonor,
  elem<int,MAX_DIMS> &DonorSize, int &NumLocalDonorPoints, array_view<int> LocalDonorPointIndices,
  array_view<long long> LocalDonorPointGridValuesIndices) {

  const connectivity_d &Donors = *Donors_;

  range DonorRange;
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    DonorRange.Begin(iDim) = Donors.Extents_(0,iDim,iDonor);
    DonorRange.End(iDim) = Donors.Extents_(1,iDim,iDonor);
  }

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    DonorSize[iDim] = DonorRange.Size(iDim);
  }

  using donor_indexer = indexer<int, int, MAX_DIMS, array_layout::GRID>;
  donor_indexer DonorIndexer(DonorRange);

  bool AwayFromEdge = GlobalRange_.Includes(DonorRange);

  if (AwayFromEdge) {
    range LocalDonorRange = IntersectRanges(LocalRange_, DonorRange);
    int iLocalDonorPoint = 0;
    for (int k = LocalDonorRange.Begin(2); k < LocalDonorRange.End(2); ++k) {
      for (int j = LocalDonorRange.Begin(1); j < LocalDonorRange.End(1); ++j) {
        for (int i = LocalDonorRange.Begin(0); i < LocalDonorRange.End(0); ++i) {
          LocalDonorPointIndices(iLocalDonorPoint) = DonorIndexer.ToIndex(i,j,k);
          LocalDonorPointGridValuesIndices(iLocalDonorPoint) = GridValuesIndexer_.ToIndex(i,j,k);
          ++iLocalDonorPoint;
        }
      }
    }
    NumLocalDonorPoints = iLocalDonorPoint;
  } else {
    int iLocalDonorPoint = 0;
    for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
      for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
        for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
          elem<int,MAX_DIMS> Point = {i,j,k};
          elem<int,MAX_DIMS> AdjustedPoint = Cart_.PeriodicAdjust(Point);
          if (LocalRange_.Contains(AdjustedPoint)) {
            LocalDonorPointIndices(iLocalDonorPoint) = DonorIndexer.ToIndex(Point);
            LocalDonorPointGridValuesIndices(iLocalDonorPoint) = GridValuesIndexer_.ToIndex(
              AdjustedPoint);
            ++iLocalDonorPoint;
          }
        }
      }
    }
    NumLocalDonorPoints = iLocalDonorPoint;
  }

}

template <typename T, array_layout Layout> void collect_base_for_type<T, Layout>::Initialize(const
  exchange &Exchange, int Count, const range &GridValuesRange) {

  parent_type::Initialize(Exchange, Count, GridValuesRange);

  int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "Collect::MemAlloc");

  int NumSends = Sends_.Count();
  int NumRecvs = Recvs_.Count();

  core::StartProfile(*Profiler_, MemAllocTime);

  SendBuffers_.Resize({NumSends});
  for (int iSend = 0; iSend < NumSends; ++iSend) {
    SendBuffers_(iSend).Resize({{Count_,Exchange.CollectSends_(iSend).NumPoints}});
  }

  if (!std::is_same<value_type, mpi_value_type>::value) {
    RecvBuffers_.Resize({NumRecvs});
    for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
      RecvBuffers_(iRecv).Resize({{Count_,Exchange.CollectRecvs_(iRecv).NumPoints}});
    }
  }

  core::EndProfile(*Profiler_, MemAllocTime);

}

template <typename T, array_layout Layout> void collect_base_for_type<T, Layout>::
  AllocateRemoteDonorValues(array<array<value_type,2>> &RemoteDonorValues) {

  int NumRecvs = Recvs_.Count();

  RemoteDonorValues.Resize({NumRecvs});
  for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
    RemoteDonorValues(iRecv).Resize({{Count_,Recvs_(iRecv).NumPoints}});
  }

}

template <typename T, array_layout Layout> void collect_base_for_type<T, Layout>::
  RetrieveRemoteDonorValues(array_view<array_view<const value_type>> GridValues, array<array<
  value_type,2>> &RemoteDonorValues) {

  MPI_Comm Comm;
  GetGridComm(*Grid_, Comm);

  int MPITime = core::GetProfilerTimerID(*Profiler_, "Collect::MPI");
  int PackTime = core::GetProfilerTimerID(*Profiler_, "Collect::Pack");

  MPI_Datatype MPIDataType = core::GetMPIDataType<mpi_value_type>();

  core::StartProfile(*Profiler_, MPITime);

  if (std::is_same<value_type, mpi_value_type>::value) {
    for (int iRecv = 0; iRecv < Recvs_.Count(); ++iRecv) {
      const exchange::collect_recv &Recv = Recvs_(iRecv);
      MPI_Request &Request = Requests_.Append();
      MPI_Irecv(RemoteDonorValues(iRecv).Data(), Count_*Recv.NumPoints, MPIDataType, Recv.Rank, 0,
        Comm, &Request);
    }
  } else {
    for (int iRecv = 0; iRecv < Recvs_.Count(); ++iRecv) {
      const exchange::collect_recv &Recv = Recvs_(iRecv);
      MPI_Request &Request = Requests_.Append();
      MPI_Irecv(RecvBuffers_(iRecv).Data(), Count_*Recv.NumPoints, MPIDataType, Recv.Rank, 0,
        Comm, &Request);
    }
  }

  core::EndProfile(*Profiler_, MPITime);
  core::StartProfile(*Profiler_, PackTime);

  for (int iSend = 0; iSend < Sends_.Count(); ++iSend) {
    const exchange::collect_send &Send = Sends_(iSend);
    for (long long iSendPoint = 0; iSendPoint < Send.NumPoints; ++iSendPoint) {
      elem<int,MAX_DIMS> Point = {
        Send.Points(0,iSendPoint),
        Send.Points(1,iSendPoint),
        Send.Points(2,iSendPoint)
      };
      long long iGridPoint = GridValuesIndexer_.ToIndex(Point);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        SendBuffers_(iSend)(iCount,iSendPoint) = mpi_value_type(GridValues(iCount)(iGridPoint));
      }
    }
  }

  core::EndProfile(*Profiler_, PackTime);
  core::StartProfile(*Profiler_, MPITime);

  for (int iSend = 0; iSend < Sends_.Count(); ++iSend) {
    const exchange::collect_send &Send = Sends_(iSend);
    MPI_Request &Request = Requests_.Append();
    MPI_Isend(SendBuffers_(iSend).Data(), Count_*Send.NumPoints, MPIDataType, Send.Rank, 0, Comm,
      &Request);
  }

  MPI_Waitall(Requests_.Count(), Requests_.Data(), MPI_STATUSES_IGNORE);

  core::EndProfile(*Profiler_, MPITime);

  if (!std::is_same<value_type, mpi_value_type>::value) {
    for (int iRecv = 0; iRecv < RecvBuffers_.Count(); ++iRecv) {
      for (long long iBuffer = 0; iBuffer < RecvBuffers_(iRecv).Count(); ++iBuffer) {
        RemoteDonorValues(iRecv)[iBuffer] = value_type(RecvBuffers_(iRecv)[iBuffer]);
      }
    }
  }

  Requests_.Clear();

}

template <typename T, array_layout Layout> void collect_base_for_type<T, Layout>::
  AssembleDonorPointValues(array_view<array_view<const value_type>> GridValues, const array<
  array<value_type,2>> &RemoteDonorValues, long long iDonor, elem<int,MAX_DIMS> &DonorSize,
  array_view<value_type,2> DonorPointValues) {

  int NumLocalDonorPoints;
  parent_type::GetLocalDonorPointInfo_(iDonor, DonorSize, NumLocalDonorPoints,
    LocalDonorPointIndices_, LocalDonorPointGridValuesIndices_);

  // Fill in the local data
  for (int iCount = 0; iCount < Count_; ++iCount) {
    for (int iLocalDonorPoint = 0; iLocalDonorPoint < NumLocalDonorPoints; ++iLocalDonorPoint) {
      DonorPointValues(iCount,LocalDonorPointIndices_(iLocalDonorPoint)) = GridValues(iCount)(
        LocalDonorPointGridValuesIndices_(iLocalDonorPoint));
    }
  }

  // Fill in the remote data
  for (int iRemotePoint = 0; iRemotePoint < NumRemoteDonorPoints_(iDonor); ++iRemotePoint) {
    int iPointInCell = RemoteDonorPoints_(iDonor)[iRemotePoint];
    int iRecv = RemoteDonorPointCollectRecvs_(iDonor)[iRemotePoint];
    long long iBuffer = RemoteDonorPointCollectRecvBufferIndices_(iDonor)[iRemotePoint];
    for (int iCount = 0; iCount < Count_; ++iCount) {
      DonorPointValues(iCount,iPointInCell) = RemoteDonorValues(iRecv)(iCount,iBuffer);
    }
  }

}

template class collect_base<array_layout::ROW_MAJOR>;
template class collect_base<array_layout::COLUMN_MAJOR>;

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

}}
