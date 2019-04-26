// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Recv.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/RecvMap.hpp"
#include "ovk/core/Request.hpp"

#include <mpi.h>

#include <utility>

namespace ovk {
namespace core {

template <typename T> class recv_request {

public:

  using value_type = T;
  using mpi_value_type = mpi_compatible_type<value_type>;

  recv_request(const recv_map &RecvMap, int Count, array<array_view<value_type>> &Values,
    array<array<mpi_value_type,2>> &Buffers, array<long long> &NextBufferEntry, array<MPI_Request>
    &MPIRequests, profiler &Profiler):
    RecvMap_(&RecvMap),
    Profiler_(&Profiler),
    MemAllocTime_(GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc")),
    MPITime_(GetProfilerTimerID(*Profiler_, "SendRecv::MPI")),
    UnpackTime_(GetProfilerTimerID(*Profiler_, "SendRecv::Unpack")),
    Count_(Count),
    Values_(Values),
    Buffers_(Buffers),
    NextBufferEntry_(NextBufferEntry),
    MPIRequests_(MPIRequests)
  {}

  array_view<MPI_Request> MPIRequests() { return MPIRequests_; }

  void Finish(int) { /* Can't finish until all requests are done */ }

  void Wait();

  void StartProfileMemAlloc() const { StartProfile(*Profiler_, MemAllocTime_); }
  void EndProfileMemAlloc() const { EndProfile(*Profiler_, MemAllocTime_); }
  void StartProfileMPI() const { StartProfile(*Profiler_, MPITime_); }
  void EndProfileMPI() const { EndProfile(*Profiler_, MPITime_); }

private:

  const recv_map *RecvMap_;
  profiler *Profiler_;
  int MemAllocTime_;
  int MPITime_;
  int UnpackTime_;
  int Count_;
  array_view<array_view<value_type>> Values_;
  array_view<array<mpi_value_type,2>> Buffers_;
  array_view<long long> NextBufferEntry_;
  array_view<MPI_Request> MPIRequests_;

};

template <typename T> class recv_impl {

public:

  using value_type = T;
  using mpi_value_type = mpi_compatible_type<value_type>;

  recv_impl(comm_view Comm, const recv_map &RecvMap, int Count, int Tag, profiler &Profiler):
    Comm_(Comm),
    RecvMap_(&RecvMap),
    Profiler_(&Profiler),
    MemAllocTime_(GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc")),
    MPITime_(GetProfilerTimerID(*Profiler_, "SendRecv::MPI")),
    Count_(Count),
    Tag_(Tag)
  {

    StartProfile(*Profiler_, MemAllocTime_);

    const array<recv_map::recv> &Recvs = RecvMap_->Recvs();

    Values_.Resize({Count_});

    Buffers_.Resize({Recvs.Count()});
    for (int iRecv = 0; iRecv < Recvs.Count(); ++iRecv) {
      const recv_map::recv &Recv = Recvs(iRecv);
      Buffers_(iRecv).Resize({{Count_,Recv.NumValues}});
    }

    NextBufferEntry_.Resize({Recvs.Count()});

    MPIRequests_.Resize({Recvs.Count()});

    EndProfile(*Profiler_, MemAllocTime_);

  }

  recv_impl(const recv_impl &Other) = delete;
  recv_impl(recv_impl &&Other) = default;

  request Recv(void **ValuesVoid) {

    long long NumValues = RecvMap_->Count();

    OVK_DEBUG_ASSERT(ValuesVoid || Count_ == 0, "Invalid values pointer.");
    for (int iCount = 0; iCount < Count_; ++iCount) {
      OVK_DEBUG_ASSERT(ValuesVoid[iCount] || NumValues == 0, "Invalid values pointer.");
      Values_(iCount) = {static_cast<value_type *>(ValuesVoid[iCount]), {NumValues}};
    }

    MPI_Datatype MPIDataType = GetMPIDataType<mpi_value_type>();

    const array<recv_map::recv> &Recvs = RecvMap_->Recvs();

    StartProfile(*Profiler_, MPITime_);

    for (int iRecv = 0; iRecv < Recvs.Count(); ++iRecv) {
      const recv_map::recv &Recv = Recvs(iRecv);
      MPI_Irecv(Buffers_(iRecv).Data(), Count_*Recv.NumValues, MPIDataType, Recv.Rank, Tag_, Comm_,
        MPIRequests_.Data(iRecv));
    }

    EndProfile(*Profiler_, MPITime_);

    return recv_request<value_type>(*RecvMap_, Count_, Values_, Buffers_, NextBufferEntry_,
      MPIRequests_, *Profiler_);

  }

private:

  comm_view Comm_;
  const recv_map *RecvMap_;
  mutable profiler *Profiler_;
  int MemAllocTime_;
  int MPITime_;
  int Count_;
  int Tag_;
  array<array_view<value_type>> Values_;
  array<array<mpi_value_type,2>> Buffers_;
  array<long long> NextBufferEntry_;
  array<MPI_Request> MPIRequests_;


};

template <typename T> void recv_request<T>::Wait() {

  long long NumValues = RecvMap_->Count();

  StartProfile(*Profiler_, MemAllocTime_);

  EndProfile(*Profiler_, MemAllocTime_);
  StartProfile(*Profiler_, MPITime_);

  MPI_Waitall(MPIRequests_.Count(), MPIRequests_.Data(), MPI_STATUSES_IGNORE);

  EndProfile(*Profiler_, MPITime_);
  StartProfile(*Profiler_, UnpackTime_);

  const array<int> &RecvOrder = RecvMap_->RecvOrder();
  const array<int> &RecvIndices = RecvMap_->RecvIndices();

  NextBufferEntry_.Fill(0);

  for (long long iOrder = 0; iOrder < NumValues; ++iOrder) {
    long long iValue = RecvOrder(iOrder);
    int iRecv = RecvIndices(iValue);
    if (iRecv >= 0) {
      long long iBuffer = NextBufferEntry_(iRecv);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        Values_(iCount)(iValue) = value_type(Buffers_(iRecv)(iCount,iBuffer));
      }
      ++NextBufferEntry_(iRecv);
    }
  }

  EndProfile(*Profiler_, UnpackTime_);

}

recv MakeRecv(comm_view Comm, const recv_map &RecvMap, data_type ValueType, int Count, int Tag,
  profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return recv_impl<bool>(Comm, RecvMap, Count, Tag, Profiler);
  case data_type::BYTE:
    return recv_impl<unsigned char>(Comm, RecvMap, Count, Tag, Profiler);
  case data_type::INT:
    return recv_impl<int>(Comm, RecvMap, Count, Tag, Profiler);
  case data_type::LONG:
    return recv_impl<long>(Comm, RecvMap, Count, Tag, Profiler);
  case data_type::LONG_LONG:
    return recv_impl<long long>(Comm, RecvMap, Count, Tag, Profiler);
  case data_type::UNSIGNED_INT:
    return recv_impl<unsigned int>(Comm, RecvMap, Count, Tag, Profiler);
  case data_type::UNSIGNED_LONG:
    return recv_impl<unsigned long>(Comm, RecvMap, Count, Tag, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return recv_impl<unsigned long long>(Comm, RecvMap, Count, Tag, Profiler);
  case data_type::FLOAT:
    return recv_impl<float>(Comm, RecvMap, Count, Tag, Profiler);
  case data_type::DOUBLE:
    return recv_impl<double>(Comm, RecvMap, Count, Tag, Profiler);
  }

  return {};

}

}}
