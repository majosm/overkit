// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_HALO_HPP_INCLUDED
#define OVK_CORE_HALO_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/DataType.hpp>
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

class exchanger {

public:

  template <typename T> exchanger(T &&Exchanger):
    Exchanger_(new model<T>(std::forward<T>(Exchanger)))
  {}

  bool Active() const { return Exchanger_->Active(); }

  request Exchange(void *ArrayDataVoid) { return Exchanger_->Exchange(ArrayDataVoid); }

private:

  struct concept {
    virtual ~concept() {}
    virtual bool Active() const = 0;
    virtual request Exchange(void *ArrayDataVoid) = 0;
  };

  template <typename T> struct model : concept {
    explicit model(T Exchanger):
      Exchanger_(std::move(Exchanger))
    {}
    virtual bool Active() const override {
      return Exchanger_.Active();
    }
    virtual request Exchange(void *ArrayDataVoid) override {
      return Exchanger_.Exchange(static_cast<typename T::value_type *>(ArrayDataVoid));
    }
    T Exchanger_;
  };

  std::unique_ptr<concept> Exchanger_;

};

}

class halo {

public:

  halo(const cart &Cart, comm_view Comm, const range &LocalRange, const range &ExtendedRange, const
    array<partition_info> &Neighbors, profiler &Profiler);

  halo(const halo &Other) = delete;
  halo(halo &&Other) noexcept = default;

  halo &operator=(const halo &Other) = delete;
  halo &operator=(halo &&Other) noexcept = default;

// Intel 17 didn't like this for some reason
//   template <typename ArrayType, OVK_FUNCDECL_REQUIRES(IsArray<ArrayType>() && ArrayHasFootprint<
//     ArrayType, MAX_DIMS, array_layout::GRID>())> request Exchange(ArrayType &Array) const;
  template <typename ArrayType, OVK_FUNCDECL_REQUIRES(IsArray<ArrayType>() && ArrayRank<ArrayType>()
    == MAX_DIMS && ArrayLayout<ArrayType>() == array_layout::GRID)> request Exchange(ArrayType &Array) const;

private:

  using exchanger = halo_internal::exchanger;

  comm_view Comm_;
  array<int> NeighborRanks_;
  array<array<long long>> NeighborSendIndices_;
  array<array<long long>> NeighborRecvIndices_;
  array<long long> LocalToLocalSourceIndices_;
  array<long long> LocalToLocalDestIndices_;
  mutable std::map<type_id_type, array<exchanger>> Exchangers_;
  mutable profiler *Profiler_;
  int TotalTime_;
  int SetupTime_;
  int ExchangeTime_;


};

namespace halo_internal {

template <typename T> class exchanger_request;

template <typename T> class exchanger_for_type {

public:

  using value_type = T;

private:

  using mpi_value_type = mpi_compatible_type<value_type>;

public:

  exchanger_for_type(comm_view Comm, const array<int> &NeighborRanks, const array<array<
    long long>> &NeighborSendIndices, const array<array<long long>> &NeighborRecvIndices, const
    array<long long> &LocalToLocalSourceIndices, const array<long long> &LocalToLocalDestIndices,
    profiler &Profiler);

  exchanger_for_type(const exchanger_for_type &Other) = delete;
  exchanger_for_type(exchanger_for_type &&Other) noexcept = default;

  exchanger_for_type &operator=(const exchanger_for_type &Other) = delete;
  exchanger_for_type &operator=(exchanger_for_type &&Other) noexcept = default;

  bool Active() const { return *Active_; }

  request Exchange(value_type *ArrayData);

private:

  comm_view Comm_;
  array_view<const int> NeighborRanks_;
  array_view<const array<long long>> NeighborSendIndices_;
  array_view<const array<long long>> NeighborRecvIndices_;
  array_view<const long long> LocalToLocalSourceIndices_;
  array_view<const long long> LocalToLocalDestIndices_;
  array<array<mpi_value_type>> SendBuffers_;
  array<array<mpi_value_type>> RecvBuffers_;
  array<MPI_Request> MPIRequests_;
  // static location so pointer in request stays valid if exchanger is moved
  std::unique_ptr<bool> Active_;
  mutable profiler *Profiler_;
  int PackTime_;
  int UnpackTime_;
  int MPITime_;

};

template <typename T> class exchanger_request {

public:

  using value_type = T;

private:

  using mpi_value_type = mpi_compatible_type<value_type>;

public:

  exchanger_request(value_type *ArrayData, array_view<const array<long long>> NeighborRecvIndices,
    array_view<array<mpi_value_type>> RecvBuffers, array_view<MPI_Request> MPIRequests,
    bool &Active, profiler &Profiler);

  array_view<MPI_Request> MPIRequests() { return MPIRequests_; }
  void Finish(int iMPIRequest);

  void Wait();

  void StartProfileMemAlloc() const {}
  void EndProfileMemAlloc() const {}
  void StartProfileMPI() const { StartProfile(*Profiler_, MPITime_); }
  void EndProfileMPI() const { EndProfile(*Profiler_, MPITime_); }

private:

  value_type *ArrayData_;
  array_view<const array<long long>> NeighborRecvIndices_;
  array_view<array<mpi_value_type>> RecvBuffers_;
  array_view<MPI_Request> MPIRequests_;
  bool *Active_;
  mutable profiler *Profiler_;
  int UnpackTime_;
  int MPITime_;

};

}

}}

#include <ovk/core/Halo.inl>

#endif
