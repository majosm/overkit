// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_GRID_H_INCLUDED
#define OVK_CORE_C_GRID_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core/Grid.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_grid;
typedef struct ovk_grid ovk_grid;

struct ovk_grid_params;
typedef struct ovk_grid_params ovk_grid_params;

struct ovk_grid_info;
typedef struct ovk_grid_info ovk_grid_info;

void ovkGetGridContextC(const ovk_grid *Grid, const ovk_context **Context);
void ovkGetGridContext(ovk_grid *Grid, ovk_context **Context);
void ovkGetGridSharedContext(ovk_grid *Grid, ovk_shared_context **Context);

void ovkGetGridName(const ovk_grid *Grid, char *Name);
void ovkGetGridDimension(const ovk_grid *Grid, int *NumDims);
void ovkGetGridComm(const ovk_grid *Grid, MPI_Comm *Comm);
void ovkGetGridCommSize(const ovk_grid *Grid, int *CommSize);
void ovkGetGridCommRank(const ovk_grid *Grid, int *CommRank);

void ovkGetGridGlobalRange(const ovk_grid *Grid, int *GlobalBegin, int *GlobalEnd);
void ovkGetGridLocalRange(const ovk_grid *Grid, int *LocalBegin, int *LocalEnd);
void ovkGetGridExtendedRange(const ovk_grid *Grid, int *ExtendedBegin, int *ExtendedEnd);

void ovkGetGridGlobalCount(const ovk_grid *Grid, long long *NumGlobal);
void ovkGetGridLocalCount(const ovk_grid *Grid, long long *NumLocal);
void ovkGetGridExtendedCount(const ovk_grid *Grid, long long *NumExtended);

void ovkGetGridSize(const ovk_grid *Grid, int *Size);

void ovkGetGridPeriodic(const ovk_grid *Grid, bool *Periodic);
void ovkGetGridPeriodicStorage(const ovk_grid *Grid, ovk_periodic_storage *PeriodicStorage);
void ovkGetGridPeriodicLength(const ovk_grid *Grid, double *PeriodicLength);

void ovkGetGridGeometryType(const ovk_grid *Grid, ovk_geometry_type *GeometryType);

void ovkCreateGridParams(ovk_grid_params **Params);
void ovkDestroyGridParams(ovk_grid_params **Params);

void ovkGetGridParamName(const ovk_grid_params *Params, char *Name);
void ovkSetGridParamName(ovk_grid_params *Params, const char *Name);
void ovkGetGridParamDimension(const ovk_grid_params *Params, int *NumDims);
void ovkSetGridParamDimension(ovk_grid_params *Params, int NumDims);
void ovkGetGridParamComm(const ovk_grid_params *Params, MPI_Comm *Comm);
void ovkSetGridParamComm(ovk_grid_params *Params, MPI_Comm Comm);
void ovkGetGridParamGlobalRange(const ovk_grid_params *Params, int *GlobalBegin, int *GlobalEnd);
void ovkSetGridParamGlobalRange(ovk_grid_params *Params, const int *GlobalBegin, const int
  *GlobalEnd);
void ovkGetGridParamLocalRange(const ovk_grid_params *Params, int *LocalBegin, int *LocalEnd);
void ovkSetGridParamLocalRange(ovk_grid_params *Params, const int *LocalBegin, const int *LocalEnd);
void ovkGetGridParamPeriodic(const ovk_grid_params *Params, bool *Periodic);
void ovkSetGridParamPeriodic(ovk_grid_params *Params, const bool *Periodic);
void ovkGetGridParamPeriodicStorage(const ovk_grid_params *Params, ovk_periodic_storage
  *PeriodicStorage);
void ovkSetGridParamPeriodicStorage(ovk_grid_params *Params, ovk_periodic_storage PeriodicStorage);
void ovkGetGridParamPeriodicLength(const ovk_grid_params *Params, double *PeriodicLength);
void ovkSetGridParamPeriodicLength(ovk_grid_params *Params, const double *PeriodicLength);
void ovkGetGridParamGeometryType(const ovk_grid_params *Params, ovk_geometry_type *GeometryType);
void ovkSetGridParamGeometryType(ovk_grid_params *Params, ovk_geometry_type GeometryType);

void ovkGetGridInfoName(const ovk_grid_info *Info, char *Name);
void ovkGetGridInfoRootRank(const ovk_grid_info *Info, int *RootRank);
void ovkGetGridInfoDimension(const ovk_grid_info *Info, int *NumDims);
void ovkGetGridInfoGlobalRange(const ovk_grid_info *Info, int *GlobalBegin, int *GlobalEnd);
void ovkGetGridInfoSize(const ovk_grid_info *Info, int *Size);
void ovkGetGridInfoPeriodic(const ovk_grid_info *Info, bool *Periodic);
void ovkGetGridInfoPeriodicStorage(const ovk_grid_info *Info, ovk_periodic_storage *PeriodicStorage);
void ovkGetGridInfoPeriodicLength(const ovk_grid_info *Info, double *PeriodicLength);
void ovkGetGridInfoGeometryType(const ovk_grid_info *Info, ovk_geometry_type *GeometryType);
void ovkGetGridInfoIsLocal(const ovk_grid_info *Info, bool *IsLocal);

#ifdef __cplusplus
}
#endif

#endif
