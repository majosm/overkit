// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Send.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <memory>
#include <string>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

namespace {

template <typename T> class send_request {

public:

  using value_type = T;
  using mpi_value_type = core::mpi_compatible_type<value_type>;

  send_request(const exchange &Exchange, int Count, int NumSends, array<array<mpi_value_type,2>>
    Buffers, array<MPI_Request> MPIRequests):
    Exchange_(&Exchange),
    Count_(Count),
    NumSends_(NumSends),
    Buffers_(std::move(Buffers)),
    MPIRequests_(std::move(MPIRequests)),
    Profiler_(Exchange.Profiler_),
    MemAllocTime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc")),
    MPITime_(core::GetProfilerTimerID(*Profiler_, "SendRecv::MPI"))
  {
    GetConnectivityDonorSide(*Exchange.Connectivity_, Donors_);
  }

  array_view<MPI_Request> MPIRequests() { return MPIRequests_; }

  void Finish(int) { /* Nothing to finish */ }

  void Wait();

  void StartProfileMemAlloc() const { core::StartProfile(*Profiler_, MemAllocTime_); }
  void EndProfileMemAlloc() const { core::EndProfile(*Profiler_, MemAllocTime_); }
  void StartProfileMPI() const { core::StartProfile(*Profiler_, MPITime_); }
  void EndProfileMPI() const { core::EndProfile(*Profiler_, MPITime_); }

private:

  const exchange *Exchange_;
  const connectivity_d *Donors_;
  int Count_;
  int NumSends_;
  array<array<mpi_value_type,2>> Buffers_;
  array<MPI_Request> MPIRequests_;
  mutable core::profiler *Profiler_;
  int MemAllocTime_;
  int MPITime_;

};

template <typename T> class send_impl {

public:

  using value_type = T;
  using mpi_value_type = core::mpi_compatible_type<value_type>;

  send_impl() = default;
  send_impl(const send_impl &) = delete;
  send_impl(send_impl &&) = default;

  void Initialize(const exchange &Exchange, int Count, int Tag) {

    Exchange_ = &Exchange;
    Count_ = Count;
    Tag_ = Tag;

    const connectivity &Connectivity = *Exchange.Connectivity_;

    GetConnectivityDonorSide(Connectivity, Donors_);

    Profiler_ = Exchange.Profiler_;

    NumSends_ = Exchange.Sends_.Count();

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc");
    core::StartProfile(*Profiler_, MemAllocTime);

    NextBufferEntry_.Resize({NumSends_});

    core::EndProfile(*Profiler_, MemAllocTime);

  }

  request Send(array_view<array_view<const value_type>> DonorValues) {

    const exchange &Exchange = *Exchange_;
    const connectivity_d &Donors = *Donors_;

    long long NumDonors;
    GetConnectivityDonorSideCount(Donors, NumDonors);

    MPI_Datatype MPIDataType = core::GetMPIDataType<mpi_value_type>();

    int MemAllocTime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MemAlloc");
    int PackTime = core::GetProfilerTimerID(*Profiler_, "SendRecv::Pack");
    int MPITime = core::GetProfilerTimerID(*Profiler_, "SendRecv::MPI");

    core::StartProfile(*Profiler_, MemAllocTime);

    // Will be moved into request object
    array<array<mpi_value_type,2>> Buffers({NumSends_});
    for (int iSend = 0; iSend < NumSends_; ++iSend) {
      const exchange::send &Send = Exchange.Sends_(iSend);
      Buffers(iSend).Resize({{Count_,Send.Count}});
    }

    core::EndProfile(*Profiler_, MemAllocTime);
    core::StartProfile(*Profiler_, PackTime);

    for (auto &iBuffer : NextBufferEntry_) {
      iBuffer = 0;
    }

    for (long long iDonorOrder = 0; iDonorOrder < NumDonors; ++iDonorOrder) {
      long long iDonor = Exchange.DonorsSorted_(iDonorOrder);
      int iSend = Exchange.DonorSendIndices_(iDonor);
      if (iSend >= 0) {
        long long iBuffer = NextBufferEntry_(iSend);
        for (int iCount = 0; iCount < Count_; ++iCount) {
          Buffers(iSend)(iCount,iBuffer) = mpi_value_type(DonorValues(iCount)(iDonor));
        }
        ++NextBufferEntry_(iSend);
      }
    }

    core::EndProfile(*Profiler_, PackTime);
    core::StartProfile(*Profiler_, MemAllocTime);

    // Will be moved into request object
    array<MPI_Request> MPIRequests({NumSends_});

    core::EndProfile(*Profiler_, MemAllocTime);
    core::StartProfile(*Profiler_, MPITime);

    for (int iSend = 0; iSend < NumSends_; ++iSend) {
      const exchange::send &Send = Exchange.Sends_(iSend);
      MPI_Isend(Buffers(iSend).Data(), Count_*Send.Count, MPIDataType, Send.Rank, Tag_,
        Exchange.Comm_, &MPIRequests(iSend));
    }

    core::EndProfile(*Profiler_, MPITime);

    return send_request<value_type>(Exchange, Count_, NumSends_, std::move(Buffers),
      std::move(MPIRequests));

  }

private:

  const exchange *Exchange_;
  const connectivity_d *Donors_;
  core::profiler *Profiler_;
  int Count_;
  int Tag_;
  int NumSends_;
  array<long long> NextBufferEntry_;

};

template <typename T> void send_request<T>::Wait() {

  core::StartProfile(*Profiler_, MPITime_);

  MPI_Waitall(NumSends_, MPIRequests_.Data(), MPI_STATUSES_IGNORE);

  core::EndProfile(*Profiler_, MPITime_);

  Buffers_.Clear();
  MPIRequests_.Clear();

}

}

send MakeSend(data_type ValueType) {

  switch (ValueType) {
  case data_type::BOOL: return send_impl<bool>();
  case data_type::BYTE: return send_impl<unsigned char>();
  case data_type::INT: return send_impl<int>();
  case data_type::LONG: return send_impl<long>();
  case data_type::LONG_LONG: return send_impl<long long>();
  case data_type::UNSIGNED_INT: return send_impl<unsigned int>();
  case data_type::UNSIGNED_LONG: return send_impl<unsigned long>();
  case data_type::UNSIGNED_LONG_LONG: return send_impl<unsigned long long>();
  case data_type::FLOAT: return send_impl<float>();
  case data_type::DOUBLE: return send_impl<double>();
  }

  return {};

}

}}
