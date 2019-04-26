// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Send.hpp"

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
#include "ovk/core/Request.hpp"
#include "ovk/core/SendMap.hpp"

#include <mpi.h>

#include <utility>

namespace ovk {
namespace core {

namespace {

template <typename T> class send_request {

public:

  using value_type = T;
  using mpi_value_type = mpi_compatible_type<value_type>;

  send_request(const send_map &SendMap, int Count, array<MPI_Request> &MPIRequests, profiler
    &Profiler):
    SendMap_(&SendMap),
    Profiler_(&Profiler),
    MemAllocTime_(GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc")),
    MPITime_(GetProfilerTimerID(*Profiler_, "SendRecv::MPI")),
    Count_(Count),
    MPIRequests_(MPIRequests)
  {}

  array_view<MPI_Request> MPIRequests() { return MPIRequests_; }

  void Finish(int) { /* Nothing to finish */ }

  void Wait();

  void StartProfileMemAlloc() const { StartProfile(*Profiler_, MemAllocTime_); }
  void EndProfileMemAlloc() const { EndProfile(*Profiler_, MemAllocTime_); }
  void StartProfileMPI() const { StartProfile(*Profiler_, MPITime_); }
  void EndProfileMPI() const { EndProfile(*Profiler_, MPITime_); }

private:

  const send_map *SendMap_;
  profiler *Profiler_;
  int MemAllocTime_;
  int MPITime_;
  int Count_;
  array_view<MPI_Request> MPIRequests_;

};

template <typename T> class send_impl {

public:

  using value_type = T;
  using mpi_value_type = mpi_compatible_type<value_type>;

  send_impl(comm_view Comm, const send_map &SendMap, int Count, int Tag, profiler &Profiler):
    Comm_(Comm),
    SendMap_(&SendMap),
    Profiler_(&Profiler),
    MemAllocTime_(GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc")),
    PackTime_(GetProfilerTimerID(*Profiler_, "SendRecv::Pack")),
    MPITime_(GetProfilerTimerID(*Profiler_, "SendRecv::MPI")),
    Count_(Count),
    Tag_(Tag)
  {

    StartProfile(*Profiler_, MemAllocTime_);

    const array<send_map::send> &Sends = SendMap.Sends();

    Values_.Resize({Count_});

    Buffers_.Resize({Sends.Count()});
    for (int iSend = 0; iSend < Sends.Count(); ++iSend) {
      const send_map::send &Send = Sends(iSend);
      Buffers_(iSend).Resize({{Count_,Send.NumValues}});
    }

    NextBufferEntry_.Resize({Sends.Count()});

    MPIRequests_.Resize({Sends.Count()});

    EndProfile(*Profiler_, MemAllocTime_);

  }

  send_impl(const send_impl &Other) = delete;
  send_impl(send_impl &&Other) = default;

  request Send(const void * const *ValuesVoid) {

    long long NumValues = SendMap_->Count();

    OVK_DEBUG_ASSERT(ValuesVoid || Count_ == 0, "Invalid values pointer.");
    for (int iCount = 0; iCount < Count_; ++iCount) {
      OVK_DEBUG_ASSERT(ValuesVoid[iCount] || NumValues == 0, "Invalid values pointer.");
      Values_(iCount) = {static_cast<const value_type *>(ValuesVoid[iCount]), {NumValues}};
    }

    MPI_Datatype MPIDataType = GetMPIDataType<mpi_value_type>();

    const array<send_map::send> &Sends = SendMap_->Sends();

    StartProfile(*Profiler_, PackTime_);

    const array<int> &SendOrder = SendMap_->SendOrder();
    const array<int> &SendIndices = SendMap_->SendIndices();

    NextBufferEntry_.Fill(0);

    for (long long iOrder = 0; iOrder < NumValues; ++iOrder) {
      long long iValue = SendOrder(iOrder);
      int iSend = SendIndices(iValue);
      if (iSend >= 0) {
        long long iBuffer = NextBufferEntry_(iSend);
        for (int iCount = 0; iCount < Count_; ++iCount) {
          Buffers_(iSend)(iCount,iBuffer) = mpi_value_type(Values_(iCount)(iValue));
        }
        ++NextBufferEntry_(iSend);
      }
    }

    EndProfile(*Profiler_, PackTime_);
    StartProfile(*Profiler_, MPITime_);

    for (int iSend = 0; iSend < Sends.Count(); ++iSend) {
      const send_map::send &Send = Sends(iSend);
      MPI_Isend(Buffers_(iSend).Data(), Count_*Send.NumValues, MPIDataType, Send.Rank, Tag_, Comm_,
        MPIRequests_.Data(iSend));
    }

    EndProfile(*Profiler_, MPITime_);

    return send_request<value_type>(*SendMap_, Count_, MPIRequests_, *Profiler_);

  }

private:

  comm_view Comm_;
  const send_map *SendMap_;
  mutable profiler *Profiler_;
  int MemAllocTime_;
  int PackTime_;
  int MPITime_;
  int Count_;
  int Tag_;
  array<array_view<const value_type>> Values_;
  array<array<mpi_value_type,2>> Buffers_;
  array<long long> NextBufferEntry_;
  array<MPI_Request> MPIRequests_;

};

template <typename T> void send_request<T>::Wait() {

  StartProfile(*Profiler_, MPITime_);

  MPI_Waitall(MPIRequests_.Count(), MPIRequests_.Data(), MPI_STATUSES_IGNORE);

  EndProfile(*Profiler_, MPITime_);

}

}

send MakeSend(comm_view Comm, const send_map &SendMap, data_type ValueType, int Count,
  int Tag, profiler &Profiler) {

  switch (ValueType) {
  case data_type::BOOL:
    return send_impl<bool>(Comm, SendMap, Count, Tag, Profiler);
  case data_type::BYTE:
    return send_impl<unsigned char>(Comm, SendMap, Count, Tag, Profiler);
  case data_type::INT:
    return send_impl<int>(Comm, SendMap, Count, Tag, Profiler);
  case data_type::LONG:
    return send_impl<long>(Comm, SendMap, Count, Tag, Profiler);
  case data_type::LONG_LONG:
    return send_impl<long long>(Comm, SendMap, Count, Tag, Profiler);
  case data_type::UNSIGNED_INT:
    return send_impl<unsigned int>(Comm, SendMap, Count, Tag, Profiler);
  case data_type::UNSIGNED_LONG:
    return send_impl<unsigned long>(Comm, SendMap, Count, Tag, Profiler);
  case data_type::UNSIGNED_LONG_LONG:
    return send_impl<unsigned long long>(Comm, SendMap, Count, Tag, Profiler);
  case data_type::FLOAT:
    return send_impl<float>(Comm, SendMap, Count, Tag, Profiler);
  case data_type::DOUBLE:
    return send_impl<double>(Comm, SendMap, Count, Tag, Profiler);
  }

  return {};

}

}}
