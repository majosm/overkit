// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Recv.hpp"

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
#include "ovk/core/RecvMap.hpp"
#include "ovk/core/Request.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

template <typename T> class recv_impl {

public:

  using value_type = T;

private:

  using mpi_value_type = mpi_compatible_type<value_type>;

  class recv_request {
  public:
    recv_request(recv_impl &Recv):
      Recv_(Recv.FloatingRefGenerator_.Generate(Recv))
    {}
    array_view<MPI_Request> MPIRequests() { return Recv_->MPIRequests_; }
    void OnMPIRequestComplete(int) {}
    void OnComplete() {

      recv_impl &Recv = *Recv_;
      const recv_map &RecvMap = *Recv.RecvMap_;

      profiler &Profiler = Recv.Context_->core_Profiler();

      long long NumValues = RecvMap.Count();

      Profiler.Start(UNPACK_TIME);

      const array<int> &RecvOrder = RecvMap.RecvOrder();
      const array<int> &RecvIndices = RecvMap.RecvIndices();

      Recv.NextBufferEntry_.Fill(0);

      for (long long iOrder = 0; iOrder < NumValues; ++iOrder) {
        long long iValue = RecvOrder(iOrder);
        int iRecv = RecvIndices(iValue);
        if (iRecv >= 0) {
          long long iBuffer = Recv.NextBufferEntry_(iRecv);
          for (int iCount = 0; iCount < Recv.Count_; ++iCount) {
            Recv.Values_(iCount)(iValue) = value_type(Recv.Buffers_(iRecv)(iCount,iBuffer));
          }
          ++Recv.NextBufferEntry_(iRecv);
        }
      }

      Profiler.Stop(UNPACK_TIME);

    }
    void StartWaitTime() const {
      profiler &Profiler = Recv_->Context_->core_Profiler();
      Profiler.Start(WAIT_TIME);
    }
    void StopWaitTime() const {
      profiler &Profiler = Recv_->Context_->core_Profiler();
      Profiler.Stop(WAIT_TIME);
    }
    void StartMPITime() const {
      profiler &Profiler = Recv_->Context_->core_Profiler();
      Profiler.Start(MPI_TIME);
    }
    void StopMPITime() const {
      profiler &Profiler = Recv_->Context_->core_Profiler();
      Profiler.Stop(MPI_TIME);
    }
  private:
    floating_ref<recv_impl> Recv_;
    static constexpr int WAIT_TIME = profiler::EXCHANGER_SEND_RECV_TIME;
  };

public:

  recv_impl(std::shared_ptr<context> &&Context, comm_view Comm, const recv_map &RecvMap, int Count,
    int Tag):
    Context_(std::move(Context)),
    Comm_(Comm),
    RecvMap_(RecvMap.GetFloatingRef()),
    Count_(Count),
    Tag_(Tag)
  {

    const array<recv_map::recv> &Recvs = RecvMap.Recvs();

    Values_.Resize({Count_});

    Buffers_.Resize({Recvs.Count()});
    for (int iRecv = 0; iRecv < Recvs.Count(); ++iRecv) {
      const recv_map::recv &Recv = Recvs(iRecv);
      Buffers_(iRecv).Resize({{Count_,Recv.NumValues}});
    }

    NextBufferEntry_.Resize({Recvs.Count()});

    MPIRequests_.Resize({Recvs.Count()});

  }

  recv_impl(const recv_impl &Other) = delete;
  recv_impl(recv_impl &&Other) = default;

  request Recv(void *ValuesVoid) {

    const recv_map &RecvMap = *RecvMap_;

    profiler &Profiler = Context_->core_Profiler();

    long long NumValues = RecvMap.Count();

    auto ValuesRaw = static_cast<value_type **>(ValuesVoid);

    OVK_DEBUG_ASSERT(ValuesRaw || Count_ == 0, "Invalid values pointer.");

    for (int iCount = 0; iCount < Count_; ++iCount) {
      OVK_DEBUG_ASSERT(ValuesRaw[iCount] || NumValues == 0, "Invalid values pointer.");
      Values_(iCount) = {ValuesRaw[iCount], {NumValues}};
    }

    MPI_Datatype MPIDataType = GetMPIDataType<mpi_value_type>();

    const array<recv_map::recv> &Recvs = RecvMap.Recvs();

    Profiler.Start(MPI_TIME);

    for (int iRecv = 0; iRecv < Recvs.Count(); ++iRecv) {
      const recv_map::recv &Recv = Recvs(iRecv);
      MPI_Irecv(Buffers_(iRecv).Data(), Count_*Recv.NumValues, MPIDataType, Recv.Rank, Tag_, Comm_,
        MPIRequests_.Data(iRecv));
    }

    Profiler.Stop(MPI_TIME);

    return recv_request(*this);

  }

private:

  floating_ref_generator FloatingRefGenerator_;

  std::shared_ptr<context> Context_;

  comm_view Comm_;

  floating_ref<const recv_map> RecvMap_;

  int Count_;
  int Tag_;

  array<array_view<value_type>> Values_;
  array<array<mpi_value_type,2>> Buffers_;
  array<long long> NextBufferEntry_;
  array<MPI_Request> MPIRequests_;

  static constexpr int MPI_TIME = profiler::EXCHANGER_SEND_RECV_MPI_TIME;
  static constexpr int UNPACK_TIME = profiler::EXCHANGER_SEND_RECV_UNPACK_TIME;

};

recv CreateRecv(std::shared_ptr<context> Context, comm_view Comm, const recv_map &RecvMap, data_type
  ValueType, int Count, int Tag) {

  recv Recv;

  switch (ValueType) {
  case data_type::BOOL:
    Recv = recv_impl<bool>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  case data_type::BYTE:
    Recv = recv_impl<unsigned char>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  case data_type::INT:
    Recv = recv_impl<int>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  case data_type::LONG:
    Recv = recv_impl<long>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  case data_type::LONG_LONG:
    Recv = recv_impl<long long>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  case data_type::UNSIGNED_INT:
    Recv = recv_impl<unsigned int>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  case data_type::UNSIGNED_LONG:
    Recv = recv_impl<unsigned long>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  case data_type::UNSIGNED_LONG_LONG:
    Recv = recv_impl<unsigned long long>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  case data_type::FLOAT:
    Recv = recv_impl<float>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  case data_type::DOUBLE:
    Recv = recv_impl<double>(std::move(Context), Comm, RecvMap, Count, Tag);
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Unhandled enum value.");
    break;
  }

  return Recv;

}

}}
