// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Recv.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Request.hpp"

#include <mpi.h>

#include <utility>

namespace ovk {
namespace core {

template <typename T> class recv_request {

public:

  using value_type = T;
  using mpi_value_type = core::mpi_compatible_type<value_type>;

  recv_request(const exchange &Exchange, int Count, int NumRecvs, array<array<mpi_value_type,2>>
    Buffers, array<MPI_Request> MPIRequests, array<array_view<value_type>> ReceiverValues):
    Exchange_(&Exchange),
    Count_(Count),
    NumRecvs_(NumRecvs),
    Buffers_(std::move(Buffers)),
    MPIRequests_(std::move(MPIRequests)),
    ReceiverValues_(std::move(ReceiverValues)),
    Profiler_(Exchange.Profiler_),
    MemAllocTime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc")),
    MPITime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::MPI")),
    UnpackTime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::Unpack"))
  {
    GetConnectivityReceiverSide(*Exchange.Connectivity_, Receivers_);
  }

  array_view<MPI_Request> MPIRequests() { return MPIRequests_; }

  void Finish(int) { /* Can't finish until all requests are done */ }

  void Wait();

  void StartProfileMemAlloc() const { core::StartProfile(*Profiler_, MemAllocTime_); }
  void EndProfileMemAlloc() const { core::EndProfile(*Profiler_, MemAllocTime_); }
  void StartProfileMPI() const { core::StartProfile(*Profiler_, MPITime_); }
  void EndProfileMPI() const { core::EndProfile(*Profiler_, MPITime_); }

private:

  const exchange *Exchange_;
  const connectivity_r *Receivers_;
  int Count_;
  int NumRecvs_;
  array<array<mpi_value_type,2>> Buffers_;
  array<MPI_Request> MPIRequests_;
  array<array_view<value_type>> ReceiverValues_;
  mutable core::profiler *Profiler_;
  int MemAllocTime_;
  int MPITime_;
  int UnpackTime_;

};

template <typename T> class recv_impl {

public:

  using value_type = T;
  using mpi_value_type = core::mpi_compatible_type<value_type>;
  using request_type = recv_request<T>;

  recv_impl() = default;
  recv_impl(const recv_impl &) = delete;
  recv_impl(recv_impl &&) = default;

  void Initialize(const exchange &Exchange, int Count, int Tag) {

    Exchange_ = &Exchange;
    Count_ = Count;
    Tag_ = Tag;

    const connectivity &Connectivity = *Exchange.Connectivity_;

    GetConnectivityReceiverSide(Connectivity, Receivers_);

    Profiler_ = Exchange.Profiler_;

    NumRecvs_ = Exchange.Recvs_.Count();

  }

  request Recv(array_view<array_view<value_type>> ReceiverValues) {

    const exchange &Exchange = *Exchange_;

    MPI_Datatype MPIDataType = core::GetMPIDataType<mpi_value_type>();

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc");
    int MPITime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MPI");

    core::StartProfile(*Profiler_, MemAllocTime);

    // Will be moved into request object
    array<array<mpi_value_type,2>> Buffers({NumRecvs_});
    for (int iRecv = 0; iRecv < NumRecvs_; ++iRecv) {
      const exchange::recv &Recv = Exchange.Recvs_(iRecv);
      Buffers(iRecv).Resize({{Count_,Recv.Count}});
    }

    // Will be moved into request object
    array<MPI_Request> MPIRequests({NumRecvs_});

    core::EndProfile(*Profiler_, MemAllocTime);
    core::StartProfile(*Profiler_, MPITime);

    for (int iRecv = 0; iRecv < NumRecvs_; ++iRecv) {
      const exchange::recv &Recv = Exchange.Recvs_(iRecv);
      MPI_Irecv(Buffers(iRecv).Data(), Count_*Recv.Count, MPIDataType, Recv.Rank, Tag_,
        Exchange.Comm_, &MPIRequests(iRecv));
    }

    core::EndProfile(*Profiler_, MPITime);
    core::StartProfile(*Profiler_, MemAllocTime);

    // Will be moved into request object
    array<array_view<value_type>> ReceiverValuesSaved({ReceiverValues});

    core::EndProfile(*Profiler_, MemAllocTime);

    return request_type(Exchange, Count_, NumRecvs_, std::move(Buffers), std::move(MPIRequests),
      std::move(ReceiverValuesSaved));

  }

private:

  const exchange *Exchange_;
  const connectivity_r *Receivers_;
  core::profiler *Profiler_;
  int Count_;
  int Tag_;
  int NumRecvs_;

};

template <typename T> void recv_request<T>::Wait() {

  const exchange &Exchange = *Exchange_;

  const connectivity_r &Receivers = *Receivers_;

  long long NumReceivers;
  GetConnectivityReceiverSideCount(Receivers, NumReceivers);

  core::StartProfile(*Profiler_, MemAllocTime_);

  array<long long> NextBufferEntry({NumRecvs_}, 0);

  core::EndProfile(*Profiler_, MemAllocTime_);
  core::StartProfile(*Profiler_, MPITime_);

  MPI_Waitall(NumRecvs_, MPIRequests_.Data(), MPI_STATUSES_IGNORE);

  core::EndProfile(*Profiler_, MPITime_);
  core::StartProfile(*Profiler_, UnpackTime_);

  for (long long iReceiverOrder = 0; iReceiverOrder < NumReceivers; ++iReceiverOrder) {
    long long iReceiver = Exchange.ReceiversSorted_(iReceiverOrder);
    int iRecv = Exchange.ReceiverRecvIndices_(iReceiver);
    if (iRecv >= 0) {
      long long iBuffer = NextBufferEntry(iRecv);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        ReceiverValues_(iCount)(iReceiver) = value_type(Buffers_(iRecv)(iCount,iBuffer));
      }
      ++NextBufferEntry(iRecv);
    }
  }

  core::EndProfile(*Profiler_, UnpackTime_);

  Buffers_.Clear();
  MPIRequests_.Clear();
  ReceiverValues_.Clear();

}

recv MakeRecv(data_type ValueType) {

  switch (ValueType) {
  case data_type::BOOL: return recv_impl<bool>();
  case data_type::BYTE: return recv_impl<unsigned char>();
  case data_type::INT: return recv_impl<int>();
  case data_type::LONG: return recv_impl<long>();
  case data_type::LONG_LONG: return recv_impl<long long>();
  case data_type::UNSIGNED_INT: return recv_impl<unsigned int>();
  case data_type::UNSIGNED_LONG: return recv_impl<unsigned long>();
  case data_type::UNSIGNED_LONG_LONG: return recv_impl<unsigned long long>();
  case data_type::FLOAT: return recv_impl<float>();
  case data_type::DOUBLE: return recv_impl<double>();
  }

  return {};

}

}}
