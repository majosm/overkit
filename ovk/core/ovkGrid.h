// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_GRID_INCLUDED
#define OVK_CORE_PUBLIC_GRID_INCLUDED

#include <ovk/core/ovkCart.h>
#include <ovk/core/ovkGlobal.h>
#include <ovk/core/ovkRange.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_grid_params;
typedef struct ovk_grid_params ovk_grid_params;

struct ovk_grid;
typedef struct ovk_grid ovk_grid;

struct ovk_grid_info;
typedef struct ovk_grid_info ovk_grid_info;

void ovkGetGridID(const ovk_grid *Grid, int *ID);
void ovkGetGridName(const ovk_grid *Grid, char *Name);
void ovkGetGridDimension(const ovk_grid *Grid, int *NumDims);
void ovkGetGridComm(const ovk_grid *Grid, MPI_Comm *Comm);
void ovkGetGridCommSize(const ovk_grid *Grid, int *CommSize);
void ovkGetGridCommRank(const ovk_grid *Grid, int *CommRank);
void ovkGetGridCart(const ovk_grid *Grid, ovk_cart *Cart);
void ovkGetGridSize(const ovk_grid *Grid, int *Size);
void ovkGetGridPeriodic(const ovk_grid *Grid, bool *Periodic);
void ovkGetGridPeriodicStorage(const ovk_grid *Grid, ovk_periodic_storage *PeriodicStorage);
void ovkGetGridPeriodicLength(const ovk_grid *Grid, double *PeriodicLength);
void ovkGetGridGeometryType(const ovk_grid *Grid, ovk_geometry_type *GeometryType);
void ovkGetGridGlobalRange(const ovk_grid *Grid, ovk_range *GlobalRange);
void ovkGetGridLocalRange(const ovk_grid *Grid, ovk_range *LocalRange);

void ovkCreateGridParams(ovk_grid_params **Params, int NumDims);
void ovkDestroyGridParams(ovk_grid_params **Params);

void ovkGetGridParamName(const ovk_grid_params *Params, char *Name);
void ovkSetGridParamName(ovk_grid_params *Params, const char *Name);
void ovkGetGridParamDimension(const ovk_grid_params *Params, int *NumDims);
void ovkGetGridParamComm(const ovk_grid_params *Params, MPI_Comm *Comm);
void ovkSetGridParamComm(ovk_grid_params *Params, MPI_Comm Comm);
void ovkGetGridParamSize(const ovk_grid_params *Params, int *Size);
void ovkSetGridParamSize(ovk_grid_params *Params, const int *Size);
void ovkGetGridParamPeriodic(const ovk_grid_params *Params, bool *Periodic);
void ovkSetGridParamPeriodic(ovk_grid_params *Params, const bool *Periodic);
void ovkGetGridParamPeriodicStorage(const ovk_grid_params *Params,
  ovk_periodic_storage *PeriodicStorage);
void ovkSetGridParamPeriodicStorage(ovk_grid_params *Params, ovk_periodic_storage PeriodicStorage);
void ovkGetGridParamPeriodicLength(const ovk_grid_params *Params, double *PeriodicLength);
void ovkSetGridParamPeriodicLength(ovk_grid_params *Params, const double *PeriodicLength);
void ovkGetGridParamGeometryType(const ovk_grid_params *Params, ovk_geometry_type *GeometryType);
void ovkSetGridParamGeometryType(ovk_grid_params *Params, ovk_geometry_type GeometryType);
void ovkSetGridParamLocalRange(ovk_grid_params *Params, const ovk_range *LocalRange);
void ovkGetGridParamLocalRange(const ovk_grid_params *Params, ovk_range *LocalRange);

void ovkGetGridInfoID(const ovk_grid_info *Info, int *ID);
void ovkGetGridInfoName(const ovk_grid_info *Info, char *Name);
void ovkGetGridInfoDimension(const ovk_grid_info *Info, int *NumDims);
void ovkGetGridInfoRootRank(const ovk_grid_info *Info, int *RootRank);
void ovkGetGridInfoCart(const ovk_grid_info *Info, ovk_cart *Cart);
void ovkGetGridInfoPeriodicLength(const ovk_grid_info *Info, double *PeriodicLength);
void ovkGetGridInfoGeometryType(const ovk_grid_info *Info, ovk_geometry_type *GeometryType);
void ovkGetGridInfoGlobalRange(const ovk_grid_info *Info, ovk_range *GlobalRange);

#ifdef __cplusplus
}
#endif

#endif
