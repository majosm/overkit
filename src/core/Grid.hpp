// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GRID_HPP_INCLUDED
#define OVK_CORE_GRID_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Partition.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/StringWrapper.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <string>

namespace ovk {

namespace grid_internal {

// For doing stuff before creation and after destruction
class grid_base {

protected:

  grid_base(std::shared_ptr<context> &&Context, std::string &&Name, MPI_Comm Comm);

  grid_base(const grid_base &Other) = delete;
  grid_base(grid_base &&Other) noexcept = default;

  grid_base &operator=(const grid_base &Other) = delete;
  grid_base &operator=(grid_base &&Other) noexcept = default;

  ~grid_base() noexcept;

  std::shared_ptr<context> Context_;

  core::string_wrapper Name_;

  comm Comm_;

};

}

class grid : private grid_internal::grid_base {

public:

  class params {
  public:
    params() = default;
    const std::string &Name() const { return *Name_; }
    params &SetName(std::string Name);
    int Dimension() const { return NumDims_; }
    params &SetDimension(int NumDims);
    MPI_Comm Comm() const { return Comm_; }
    params &SetComm(MPI_Comm Comm);
    const cart &Cart() const { return Cart_; }
    params &SetCart(const cart &Cart);
    const range &GlobalRange() const { return Cart_.Range(); }
    params &SetGlobalRange(const range &GlobalRange);
    const range &LocalRange() const { return LocalRange_; }
    params &SetLocalRange(const range &LocalRange);
    const tuple<bool> &Periodic() const { return Cart_.Periodic(); }
    params &SetPeriodic(const tuple<bool> &Periodic);
    periodic_storage PeriodicStorage() const { return Cart_.PeriodicStorage(); }
    params &SetPeriodicStorage(periodic_storage PeriodicStorage);
  private:
    core::string_wrapper Name_ = "Grid";
    int NumDims_ = 2;
    MPI_Comm Comm_ = MPI_COMM_NULL;
    cart Cart_ = cart(2, {{1,1,1}}, {false, false, false}, periodic_storage::UNIQUE);
    range LocalRange_ = {{0,0,0}, {1,1,1}};
    friend class grid;
  };

  grid(const grid &Other) = delete;
  grid(grid &&Other) noexcept = default;

  grid &operator=(const grid &Other) = delete;
  grid &operator=(grid &&Other) noexcept = default;

  ~grid() noexcept;

  floating_ref<const grid> GetFloatingRef() const { return FloatingRefGenerator_.Generate(*this); }
  floating_ref<grid> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const std::string &Name() const { return *Name_; }
  int Dimension() const { return NumDims_; }

  const comm &Comm() const { return Comm_; }

  const cart &Cart() const { return Cart_; }

  const range &GlobalRange() const { return Cart_.Range(); }
  const range &LocalRange() const { return Partition_->LocalRange(); }
  const range &ExtendedRange() const { return Partition_->ExtendedRange(); }

  tuple<int> Size() const { return Cart_.Range().Size(); }
  int Size(int iDim) const { return Cart_.Range().Size(iDim); }

  const tuple<bool> &Periodic() const { return Cart_.Periodic(); }
  bool Periodic(int iDim) const { return Cart_.Periodic(iDim); }
  periodic_storage PeriodicStorage() const { return Cart_.PeriodicStorage(); }

  const core::partition_hash &core_PartitionHash() const { return PartitionHash_; }
  const std::shared_ptr<core::partition> &core_PartitionShared() const { return Partition_; }
  const core::partition &core_Partition() const { return *Partition_; }

  static grid internal_Create(std::shared_ptr<context> &&Context, params &&Params);

private:

  floating_ref_generator FloatingRefGenerator_;

  int NumDims_;
  cart Cart_;

  core::partition_hash PartitionHash_;
  std::shared_ptr<core::partition> Partition_;

  grid(std::shared_ptr<context> &&Context, params &&Params);

};

namespace core {

grid CreateGrid(std::shared_ptr<context> Context, grid::params Params);

}

class grid_info {

public:

  grid_info() = default;

  const std::string &Name() const { return *Name_; }
  int RootRank() const { return RootRank_; }
  const cart &Cart() const { return Cart_; }
  const range &GlobalRange() const { return Cart_.Range(); }
  tuple<int> Size() const { return Cart_.Range().Size(); }
  int Size(int iDim) const { return Cart_.Range().Size(iDim); }
  bool IsLocal() const { return IsLocal_; }

  static grid_info internal_Create(grid *MaybeGrid, comm_view Comm);

private:

  core::string_wrapper Name_;
  int RootRank_ = -1;
  cart Cart_ = MakeEmptyCart(2);
  bool IsLocal_ = false;

  grid_info(grid *MaybeGrid, comm_view Comm);

};

namespace core {
grid_info CreateGridInfo(grid *MaybeGrid, comm_view Comm);
}

}

#endif
