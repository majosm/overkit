// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GRID_HPP_INCLUDED
#define OVK_CORE_GRID_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/GridBase.h>
#include <ovk/core/Partition.hpp>
#include <ovk/core/PartitionHash.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/StringWrapper.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <string>

namespace ovk {

enum class geometry_type {
  UNIFORM = OVK_GEOMETRY_UNIFORM,
  ORIENTED_UNIFORM = OVK_GEOMETRY_ORIENTED_UNIFORM,
  RECTILINEAR = OVK_GEOMETRY_RECTILINEAR,
  ORIENTED_RECTILINEAR = OVK_GEOMETRY_ORIENTED_RECTILINEAR,
  CURVILINEAR = OVK_GEOMETRY_CURVILINEAR
};

inline bool ValidGeometryType(geometry_type GeometryType) {
  return ovkValidGeometryType(ovk_geometry_type(GeometryType));
}

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
    const tuple<int> &Size() const { return Size_; }
    params &SetSize(const tuple<int> &Size);
    const tuple<bool> &Periodic() const { return Periodic_; }
    params &SetPeriodic(const tuple<bool> &Periodic);
    periodic_storage PeriodicStorage() const { return PeriodicStorage_; }
    params &SetPeriodicStorage(periodic_storage PeriodicStorage);
    const tuple<double> &PeriodicLength() const { return PeriodicLength_; }
    params &SetPeriodicLength(const tuple<double> &PeriodicLength);
    geometry_type GeometryType() const { return GeometryType_; }
    params &SetGeometryType(geometry_type GeometryType);
    const range &LocalRange() const { return LocalRange_; }
    params &SetLocalRange(const range &LocalRange);
  private:
    core::string_wrapper Name_ = "Grid";
    int NumDims_ = 2;
    MPI_Comm Comm_ = MPI_COMM_NULL;
    tuple<int> Size_ = {1,1,1};
    tuple<bool> Periodic_ = {false, false, false};
    periodic_storage PeriodicStorage_ = periodic_storage::UNIQUE;
    tuple<double> PeriodicLength_ = {0., 0., 0.};
    geometry_type GeometryType_ = geometry_type::CURVILINEAR;
    range LocalRange_ = {{0,0,0}, {1,1,1}};
    friend class grid;
  };

  grid(const grid &Other) = delete;
  grid(grid &&Other) noexcept = default;

  grid &operator=(const grid &Other) = delete;
  grid &operator=(grid &&Other) noexcept = default;

  ~grid() noexcept;

  floating_ref<const grid> GetFloatingRef() const { return FloatingRefGenerator_.Generate(); }
  floating_ref<grid> GetFloatingRef() { return FloatingRefGenerator_.Generate(); }

  const context &Context() const { return *Context_; }
  context &Context() { return *Context_; }
  const std::shared_ptr<context> &SharedContext() const { return Context_; }

  const std::string &Name() const { return *Name_; }
  int Dimension() const { return NumDims_; }

  const comm &Comm() const { return Comm_; }

  const cart &Cart() const { return Cart_; }
  tuple<int> Size() const { return Cart_.Range().Size(); }
  int Size(int iDim) const { return Cart_.Range().Size(iDim); }
  const tuple<bool> &Periodic() const { return Cart_.Periodic(); }
  bool Periodic(int iDim) const { return Cart_.Periodic(iDim); }
  periodic_storage PeriodicStorage() const { return Cart_.PeriodicStorage(); }
  const tuple<double> &PeriodicLength() const { return PeriodicLength_; }
  double PeriodicLength(int iDim) const { return PeriodicLength_[iDim]; }

  geometry_type GeometryType() const { return GeometryType_; }

  const range &GlobalRange() const { return Cart_.Range(); }
  const range &LocalRange() const { return Partition_->LocalRange(); }

  const core::partition_hash &core_PartitionHash() const { return PartitionHash_; }
  const std::shared_ptr<core::partition> &core_PartitionShared() const { return Partition_; }
  const core::partition &core_Partition() const { return *Partition_; }

  static grid internal_Create(std::shared_ptr<context> &&Context, params &&Params);

private:

  floating_ref_generator<grid> FloatingRefGenerator_;

  int NumDims_;
  cart Cart_;
  tuple<double> PeriodicLength_;
  geometry_type GeometryType_;

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
  const tuple<double> &PeriodicLength() const { return PeriodicLength_; }
  geometry_type GeometryType() const { return GeometryType_; }
  const range &GlobalRange() const { return Cart_.Range(); }
  bool IsLocal() const { return IsLocal_; }

  static grid_info internal_Create(grid *MaybeGrid, comm_view Comm);

private:

  core::string_wrapper Name_;
  int RootRank_ = -1;
  cart Cart_ = MakeEmptyCart(2);
  tuple<double> PeriodicLength_ = MakeUniformTuple<double>(2, 0.);
  geometry_type GeometryType_ = geometry_type::CURVILINEAR;
  bool IsLocal_ = false;

  grid_info(grid *MaybeGrid, comm_view Comm);

};

namespace core {
grid_info CreateGridInfo(grid *MaybeGrid, comm_view Comm);
}

}

#endif
