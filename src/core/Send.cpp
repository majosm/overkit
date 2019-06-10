// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Send.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/SendMap.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

namespace {

template <typename T> class send_impl {

public:

  using value_type = T;

private:

  using mpi_value_type = mpi_compatible_type<value_type>;

  class send_request {
  public:
    send_request(send_impl &Send):
      Send_(Send.FloatingRefGenerator_.Generate())
    {}
    array_view<MPI_Request> MPIRequests() { return Send_->MPIRequests_; }
    void Finish(int) { /* Nothing to finish */ }
    void Wait() {

      send_impl &Send = *Send_;

      profiler &Profiler = Send.Context_->core_Profiler();

      Profiler.Start(WAIT_TIME);
      Profiler.Start(MPI_TIME);

      MPI_Waitall(Send.MPIRequests_.Count(), Send.MPIRequests_.Data(), MPI_STATUSES_IGNORE);

      Profiler.Stop(MPI_TIME);
      Profiler.Stop(WAIT_TIME);

    }
    void StartWaitTime() const {
      profiler &Profiler = Send_->Context_->core_Profiler();
      Profiler.Start(WAIT_TIME);
    }
    void StopWaitTime() const {
      profiler &Profiler = Send_->Context_->core_Profiler();
      Profiler.Stop(WAIT_TIME);
    }
    void StartMPITime() const {
      profiler &Profiler = Send_->Context_->core_Profiler();
      Profiler.Start(MPI_TIME);
    }
    void StopMPITime() const {
      profiler &Profiler = Send_->Context_->core_Profiler();
      Profiler.Stop(MPI_TIME);
    }
  private:
    floating_ref<send_impl> Send_;
    static constexpr int WAIT_TIME = profiler::EXCHANGER_SEND_RECV_TIME;
  };

public:

  send_impl(std::shared_ptr<context> &&Context, comm_view Comm, const send_map &SendMap, int Count,
    int Tag):
    FloatingRefGenerator_(*this),
    Context_(std::move(Context)),
    Comm_(Comm),
    SendMap_(SendMap.GetFloatingRef()),
    Count_(Count),
    Tag_(Tag)
  {

    const array<send_map::send> &Sends = SendMap.Sends();

    Values_.Resize({Count_});

    Buffers_.Resize({Sends.Count()});
    for (int iSend = 0; iSend < Sends.Count(); ++iSend) {
      const send_map::send &Send = Sends(iSend);
      Buffers_(iSend).Resize({{Count_,Send.NumValues}});
    }

    NextBufferEntry_.Resize({Sends.Count()});

    MPIRequests_.Resize({Sends.Count()});

  }

  send_impl(const send_impl &Other) = delete;
  send_impl(send_impl &&Other) = default;

  request Send(const void * const *ValuesVoid) {

    const send_map &SendMap = *SendMap_;

    profiler &Profiler = Context_->core_Profiler();

    long long NumValues = SendMap.Count();

    OVK_DEBUG_ASSERT(ValuesVoid || Count_ == 0, "Invalid values pointer.");
    for (int iCount = 0; iCount < Count_; ++iCount) {
      OVK_DEBUG_ASSERT(ValuesVoid[iCount] || NumValues == 0, "Invalid values pointer.");
      Values_(iCount) = {static_cast<const value_type *>(ValuesVoid[iCount]), {NumValues}};
    }

    MPI_Datatype MPIDataType = GetMPIDataType<mpi_value_type>();

    const array<send_map::send> &Sends = SendMap.Sends();

    Profiler.Start(PACK_TIME);

    const array<int> &SendOrder = SendMap.SendOrder();
    const array<int> &SendIndices = SendMap.SendIndices();

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

    Profiler.Stop(PACK_TIME);
    Profiler.Start(MPI_TIME);

    for (int iSend = 0; iSend < Sends.Count(); ++iSend) {
      const send_map::send &Send = Sends(iSend);
      MPI_Isend(Buffers_(iSend).Data(), Count_*Send.NumValues, MPIDataType, Send.Rank, Tag_, Comm_,
        MPIRequests_.Data(iSend));
    }

    Profiler.Stop(MPI_TIME);

    return send_request(*this);

  }

private:

  floating_ref_generator<send_impl> FloatingRefGenerator_;

  std::shared_ptr<context> Context_;

  comm_view Comm_;

  floating_ref<const send_map> SendMap_;

  int Count_;
  int Tag_;

  array<array_view<const value_type>> Values_;
  array<array<mpi_value_type,2>> Buffers_;
  array<long long> NextBufferEntry_;
  array<MPI_Request> MPIRequests_;

  static constexpr int PACK_TIME = profiler::EXCHANGER_SEND_RECV_PACK_TIME;
  static constexpr int MPI_TIME = profiler::EXCHANGER_SEND_RECV_MPI_TIME;

};

}

send CreateSend(std::shared_ptr<context> Context, comm_view Comm, const send_map &SendMap, data_type
  ValueType, int Count, int Tag) {

  send Send;

  switch (ValueType) {
  case data_type::BOOL:
    Send = send_impl<bool>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  case data_type::BYTE:
    Send = send_impl<unsigned char>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  case data_type::INT:
    Send = send_impl<int>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  case data_type::LONG:
    Send = send_impl<long>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  case data_type::LONG_LONG:
    Send = send_impl<long long>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  case data_type::UNSIGNED_INT:
    Send = send_impl<unsigned int>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  case data_type::UNSIGNED_LONG:
    Send = send_impl<unsigned long>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Send = send_impl<unsigned long long>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  case data_type::FLOAT:
    Send = send_impl<float>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  case data_type::DOUBLE:
    Send = send_impl<double>(std::move(Context), Comm, SendMap, Count, Tag);
    break;
  }

  return Send;

}

}}
