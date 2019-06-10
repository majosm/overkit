// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_HALO_HPP_INCLUDED
#define OVK_CORE_HALO_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/PartitionBase.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Request.hpp>
#include <ovk/core/ScopeGuard.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <map>
#include <memory>
#include <utility>

namespace ovk {
namespace core {

namespace halo_internal {

class halo_map {

public:

  halo_map();
  halo_map(const cart &Cart, const range &LocalRange, const range &ExtendedRange,
    const array<partition_info> &Neighbors);

  floating_ref<const halo_map> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate();
  }

  floating_ref<halo_map> GetFloatingRef() {
    return FloatingRefGenerator_.Generate();
  }

  const array<int> &NeighborRanks() const { return NeighborRanks_; }

  const array<long long> &NeighborSendIndices(int iNeighbor) const {
    return NeighborSendIndices_(iNeighbor);
  }

  const array<long long> &NeighborRecvIndices(int iNeighbor) const {
    return NeighborRecvIndices_(iNeighbor);
  }

  const array<long long> &LocalToLocalSourceIndices() const { return LocalToLocalSourceIndices_; }

  const array<long long> &LocalToLocalDestIndices() const { return LocalToLocalDestIndices_; }

private:

  floating_ref_generator<halo_map> FloatingRefGenerator_;

  array<int> NeighborRanks_;
  array<array<long long>> NeighborSendIndices_;
  array<array<long long>> NeighborRecvIndices_;
  array<long long> LocalToLocalSourceIndices_;
  array<long long> LocalToLocalDestIndices_;

};

class halo_exchanger {

public:

  template <typename T> halo_exchanger(T &&HaloExchanger):
    HaloExchanger_(new model<T>(std::forward<T>(HaloExchanger)))
  {}

  bool Active() const { return HaloExchanger_->Active(); }

  request Exchange(void *ArrayDataVoid) { return HaloExchanger_->Exchange(ArrayDataVoid); }

private:

  struct concept {
    virtual ~concept() noexcept {}
    virtual bool Active() const = 0;
    virtual request Exchange(void *ArrayDataVoid) = 0;
  };

  template <typename T> struct model final : concept {
    explicit model(T HaloExchanger):
      HaloExchanger_(std::move(HaloExchanger))
    {}
    virtual bool Active() const override {
      return HaloExchanger_.Active();
    }
    virtual request Exchange(void *ArrayDataVoid) override {
      return HaloExchanger_.Exchange(static_cast<typename T::value_type *>(ArrayDataVoid));
    }
    T HaloExchanger_;
  };

  std::unique_ptr<concept> HaloExchanger_;

};

}

class halo {

public:

  halo(std::shared_ptr<context> Context, const cart &Cart, comm_view Comm, const range &LocalRange,
    const range &ExtendedRange, const array<partition_info> &Neighbors);

  halo(const halo &Other) = delete;
  halo(halo &&Other) noexcept = default;

  halo &operator=(const halo &Other) = delete;
  halo &operator=(halo &&Other) noexcept = default;

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

// Intel 17 didn't like this for some reason
//   template <typename ArrayType, OVK_FUNCDECL_REQUIRES(IsArray<ArrayType>() && ArrayHasFootprint<
//     ArrayType, MAX_DIMS, array_layout::GRID>())> request Exchange(ArrayType &Array) const;
  template <typename ArrayType, OVK_FUNCDECL_REQUIRES(IsArray<ArrayType>() && ArrayRank<ArrayType>()
    == MAX_DIMS && ArrayLayout<ArrayType>() == array_layout::GRID)> request Exchange(ArrayType
    &Array) const;

private:

  using halo_map = halo_internal::halo_map;
  using halo_exchanger = halo_internal::halo_exchanger;

  std::shared_ptr<context> Context_;

  comm_view Comm_;

  halo_map HaloMap_;

  mutable std::map<type_id_type, array<halo_exchanger>> HaloExchangers_;

  static constexpr int TOTAL_TIME = profiler::HALO_TIME;
  static constexpr int SETUP_TIME = profiler::HALO_SETUP_TIME;
  static constexpr int EXCHANGE_TIME = profiler::HALO_EXCHANGE_TIME;

};

namespace halo_internal {

template <typename T> class halo_exchanger_for_type {

public:

  using value_type = T;

private:

  using mpi_value_type = mpi_compatible_type<value_type>;

public:

  halo_exchanger_for_type(context &Context, comm_view Comm, const halo_map &HaloMap);

  halo_exchanger_for_type(const halo_exchanger_for_type &Other) = delete;
  halo_exchanger_for_type(halo_exchanger_for_type &&Other) noexcept = default;

  halo_exchanger_for_type &operator=(const halo_exchanger_for_type &Other) = delete;
  halo_exchanger_for_type &operator=(halo_exchanger_for_type &&Other) noexcept = default;

  bool Active() const { return Active_; }

  request Exchange(value_type *ArrayData);

private:

  class exchange_request {
  public:
    exchange_request(halo_exchanger_for_type &HaloExchanger, value_type *ArrayData);
    array_view<MPI_Request> MPIRequests() { return HaloExchanger_->MPIRequests_; }
    void Finish(int iMPIRequest);
    void Wait();
    void StartWaitTime() const {
      profiler &Profiler = HaloExchanger_->Context_->core_Profiler();
      Profiler.Start(WAIT_TIME);
    }
    void StopWaitTime() const {
      profiler &Profiler = HaloExchanger_->Context_->core_Profiler();
      Profiler.Stop(WAIT_TIME);
    }
    void StartMPITime() const {
      profiler &Profiler = HaloExchanger_->Context_->core_Profiler();
      Profiler.Start(MPI_TIME);
    }
    void StopMPITime() const {
      profiler &Profiler = HaloExchanger_->Context_->core_Profiler();
      Profiler.Stop(MPI_TIME);
    }
  private:
    floating_ref<halo_exchanger_for_type> HaloExchanger_;
    value_type *ArrayData_;
    static constexpr int WAIT_TIME = profiler::HALO_EXCHANGE_TIME;
  };

  floating_ref_generator<halo_exchanger_for_type> FloatingRefGenerator_;

  floating_ref<context> Context_;

  comm_view Comm_;

  floating_ref<const halo_map> HaloMap_;

  array<array<mpi_value_type>> SendBuffers_;
  array<array<mpi_value_type>> RecvBuffers_;
  array<MPI_Request> MPIRequests_;

  bool Active_ = false;

  static constexpr int PACK_TIME = profiler::HALO_EXCHANGE_PACK_TIME;
  static constexpr int MPI_TIME = profiler::HALO_EXCHANGE_MPI_TIME;
  static constexpr int UNPACK_TIME = profiler::HALO_EXCHANGE_UNPACK_TIME;

};

}

}}

#include <ovk/core/Halo.inl>

#endif
