// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GRID_HPP_INCLUDED
#define OVK_CORE_GRID_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/ErrorHandler.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Logger.hpp>
#include <ovk/core/Moveabool.hpp>
#include <ovk/core/Partition.hpp>
#include <ovk/core/PartitionHash.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <string>

namespace ovk {

struct grid_params {
  std::string Name_;
  int NumDims_;
  MPI_Comm Comm_;
  tuple<int> Size_;
  tuple<bool> Periodic_;
  periodic_storage PeriodicStorage_;
  tuple<double> PeriodicLength_;
  geometry_type GeometryType_;
  range LocalRange_;
};

namespace grid_internal {

// For doing stuff before grid creation and after grid destruction
class grid_base {

protected:

  grid_base(core::logger &Logger, core::error_handler &ErrorHandler, core::profiler &Profiler,
    int ID, const std::string &Name, int NumDims, MPI_Comm Comm);

  grid_base(const grid_base &Other) = delete;
  grid_base(grid_base &&Other) noexcept = default;

  grid_base &operator=(const grid_base &Other) = delete;
  grid_base &operator=(grid_base &&Other) noexcept = default;

  ~grid_base();

  core::moveabool Exists_;

  mutable core::logger *Logger_;
  mutable core::error_handler *ErrorHandler_;
  mutable core::profiler *Profiler_;

  int ID_;
  std::string Name_;
  int NumDims_;
  core::comm Comm_;

};

}

class grid : private grid_internal::grid_base {

public:

  grid(int ID, const grid_params &Params, core::logger &Logger, core::error_handler &ErrorHandler,
    core::profiler &Profiler);

  grid(const grid &Other) = delete;
  grid(grid &&Other) noexcept = default;

  grid &operator=(const grid &Other) = delete;
  grid &operator=(grid &&Other) noexcept = default;

  ~grid();

  int ID() const { return ID_; }
  const std::string &Name() const { return Name_; }
  int Dimension() const { return NumDims_; }
  MPI_Comm Comm() const { return Comm_.Get(); }
  int CommSize() const { return Comm_.Size(); }
  int CommRank() const { return Comm_.Rank(); }
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

  const core::comm &core_Comm() const { return Comm_; }
  const core::partition_hash &core_PartitionHash() const { return PartitionHash_; }
  const core::partition &core_Partition() const { return *Partition_; }
  const std::shared_ptr<core::partition> &core_PartitionShared() const { return Partition_; }
  core::logger &core_Logger() const { return *Logger_; }
  core::error_handler &core_ErrorHandler() const { return *ErrorHandler_; }
  core::profiler &core_Profiler() const { return *Profiler_; }

private:

  cart Cart_;
  tuple<double> PeriodicLength_;
  geometry_type GeometryType_;
  core::partition_hash PartitionHash_;
  std::shared_ptr<core::partition> Partition_;

  void PrintSummary_() const;
  void PrintDecomposition_() const;

};

struct grid_info {
  int ID_;
  std::string Name_;
  int RootRank_;
  cart Cart_;
  tuple<double> PeriodicLength_;
  geometry_type GeometryType_;
  bool IsLocal_;
};

void CreateGridParams(grid_params &Params, int NumDims);
void DestroyGridParams(grid_params &Params);

void GetGridParamName(const grid_params &Params, std::string &Name);
void SetGridParamName(grid_params &Params, std::string Name);
void GetGridParamDimension(const grid_params &Params, int &NumDims);
void GetGridParamComm(const grid_params &Params, MPI_Comm &Comm);
void SetGridParamComm(grid_params &Params, MPI_Comm Comm);
void GetGridParamSize(const grid_params &Params, int *Size);
void SetGridParamSize(grid_params &Params, const int *Size);
void GetGridParamPeriodic(const grid_params &Params, bool *Periodic);
void SetGridParamPeriodic(grid_params &Params, const bool *Periodic);
void GetGridParamPeriodicStorage(const grid_params &Params, periodic_storage &PeriodicStorage);
void SetGridParamPeriodicStorage(grid_params &Params, periodic_storage PeriodicStorage);
void GetGridParamPeriodicLength(const grid_params &Params, double *PeriodicLength);
void SetGridParamPeriodicLength(grid_params &Params, const double *PeriodicLength);
void GetGridParamGeometryType(const grid_params &Params, geometry_type &GeometryType);
void SetGridParamGeometryType(grid_params &Params, geometry_type GeometryType);
void SetGridParamLocalRange(grid_params &Params, const range &LocalRange);
void GetGridParamLocalRange(const grid_params &Params, range &LocalRange);

namespace core {
void CreateGridInfo(grid_info &Info, const grid *Grid, comm_view Comm);
void DestroyGridInfo(grid_info &Info);
}

void GetGridInfoID(const grid_info &Info, int &ID);
void GetGridInfoName(const grid_info &Info, std::string &Name);
void GetGridInfoRootRank(const grid_info &Info, int &RootRank);
void GetGridInfoCart(const grid_info &Info, cart &Cart);
void GetGridInfoPeriodicLength(const grid_info &Info, double *PeriodicLength);
void GetGridInfoGeometryType(const grid_info &Info, geometry_type &GeometryType);
void GetGridInfoGlobalRange(const grid_info &Info, range &GlobalRange);
void GetGridInfoIsLocal(const grid_info &Info, bool &IsLocal);

}

#endif
