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
#include <ovk/core/PartitionHash.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

#include <string>

namespace ovk {

struct grid_params {
  std::string Name_;
  int NumDims_;
  MPI_Comm Comm_;
  int Size_[MAX_DIMS];
  bool Periodic_[MAX_DIMS];
  periodic_storage PeriodicStorage_;
  double PeriodicLength_[MAX_DIMS];
  geometry_type GeometryType_;
  range LocalRange_;
};

namespace core {

struct grid_neighbor {
  int Rank;
  range LocalRange;
};

}

struct grid {
  mutable core::logger *Logger_;
  mutable core::error_handler *ErrorHandler_;
  int ID_;
  std::string Name_;
  int NumDims_;
  core::comm Comm_;
  cart Cart_;
  periodic_storage PeriodicStorage_;
  double PeriodicLength_[MAX_DIMS];
  geometry_type GeometryType_;
  range GlobalRange_;
  range LocalRange_;
  array<core::grid_neighbor> Neighbors_;
  core::partition_hash PartitionHash_;
};

struct grid_info {
  int ID_;
  std::string Name_;
  int NumDims_;
  int RootRank_;
  cart Cart_;
  double PeriodicLength_[MAX_DIMS];
  geometry_type GeometryType_;
  range GlobalRange_;
  bool IsLocal_;
};

namespace core {
void CreateGrid(grid &Grid, int ID, const grid_params &Params, core::logger &Logger,
  core::error_handler &ErrorHandler);
void DestroyGrid(grid &Grid);
}

void GetGridID(const grid &Grid, int &ID);
void GetGridName(const grid &Grid, std::string &Name);
void GetGridDimension(const grid &Grid, int &NumDims);
void GetGridComm(const grid &Grid, MPI_Comm &Comm);
void GetGridCommSize(const grid &Grid, int &CommSize);
void GetGridCommRank(const grid &Grid, int &CommRank);
void GetGridCart(const grid &Grid, cart &Cart);
void GetGridSize(const grid &Grid, int *Size);
void GetGridPeriodic(const grid &Grid, bool *Periodic);
void GetGridPeriodicStorage(const grid &Grid, periodic_storage &PeriodicStorage);
void GetGridPeriodicLength(const grid &Grid, double *PeriodicLength);
void GetGridGeometryType(const grid &Grid, geometry_type &GeometryType);
void GetGridGlobalRange(const grid &Grid, range &GlobalRange);
void GetGridLocalRange(const grid &Grid, range &LocalRange);
void GetGridGlobalCount(const grid &Grid, long long &NumGlobal);
void GetGridLocalCount(const grid &Grid, long long &NumLocal);

namespace core {
const comm &GetGridComm(const grid &Grid);
const array<grid_neighbor> &GetGridNeighbors(const grid &Grid);
const partition_hash &GetGridPartitionHash(const grid &Grid);
}

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
void CreateGridInfo(grid_info &Info, const grid *Grid, const comm &Comm);
void DestroyGridInfo(grid_info &Info);
}

void GetGridInfoID(const grid_info &Info, int &ID);
void GetGridInfoName(const grid_info &Info, std::string &Name);
void GetGridInfoDimension(const grid_info &Info, int &NumDims);
void GetGridInfoRootRank(const grid_info &Info, int &RootRank);
void GetGridInfoCart(const grid_info &Info, cart &Cart);
void GetGridInfoPeriodicLength(const grid_info &Info, double *PeriodicLength);
void GetGridInfoGeometryType(const grid_info &Info, geometry_type &GeometryType);
void GetGridInfoGlobalRange(const grid_info &Info, range &GlobalRange);
void GetGridInfoIsLocal(const grid_info &Info, bool &IsLocal);

// namespace core {
// static inline void GetGridLogger(const ovk_grid *Grid, t_logger **Logger) {
//   *Logger = (t_logger *)Grid->logger;
// }
// static inline void GetGridErrorHandler(const ovk_grid *Grid, t_error_handler **ErrorHandler) {
//   *ErrorHandler = (t_error_handler *)Grid->error_handler;
// }
// }

}

#endif
