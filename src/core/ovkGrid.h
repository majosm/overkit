// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_GRID_INCLUDED
#define OVK_CORE_PUBLIC_GRID_INCLUDED

#include <ovkGlobal.h>

struct ovk_grid_params;
typedef struct ovk_grid_params ovk_grid_params;

struct ovk_grid;
typedef struct ovk_grid ovk_grid;

struct ovk_grid_properties;
typedef struct ovk_grid_properties ovk_grid_properties;

void ovkGetGridProperties(const ovk_grid *Grid, const ovk_grid_properties **Properties);
// void ovkEditGridProperties(ovk_grid *Grid, ovk_grid_properties **Properties);
// void ovkReleaseGridProperties(ovk_grid *Grid, ovk_grid_properties **Properties);

void ovkGetGridParamID(const ovk_grid_params *Params, int *ID);
void ovkSetGridParamID(ovk_grid_params *Params, int ID);
void ovkGetGridParamGlobalSize(const ovk_grid_params *Params, int *GlobalSize);
void ovkSetGridParamGlobalSize(ovk_grid_params *Params, const int *GlobalSize);
void ovkGetGridParamLocalStart(const ovk_grid_params *Params, int *LocalStart);
void ovkSetGridParamLocalStart(ovk_grid_params *Params, const int *LocalStart);
void ovkGetGridParamLocalEnd(const ovk_grid_params *Params, int *LocalEnd);
void ovkSetGridParamLocalEnd(ovk_grid_params *Params, const int *LocalEnd);
void ovkGetGridParamPeriodic(const ovk_grid_params *Params, bool *Periodic);
void ovkSetGridParamPeriodic(ovk_grid_params *Params, const bool *Periodic);
void ovkGetGridParamPeriodicStorage(const ovk_grid_params *Params,
  ovk_periodic_storage *PeriodicStorage);
void ovkSetGridParamPeriodicStorage(ovk_grid_params *Params, ovk_periodic_storage PeriodicStorage);
void ovkGetGridParamPeriodicLength(const ovk_grid_params *Params, double *PeriodicLength);
void ovkSetGridParamPeriodicLength(ovk_grid_params *Params, const double *PeriodicLength);
void ovkGetGridParamGeometryType(const ovk_grid_params *Params, ovk_geometry_type *GeometryType);
void ovkSetGridParamGeometryType(ovk_grid_params *Params, ovk_geometry_type GeometryType);
void ovkGetGridParamComm(const ovk_grid_params *Params, MPI_Comm *Comm);
void ovkSetGridParamComm(ovk_grid_params *Params, MPI_Comm Comm);
void ovkGetGridParamNumNeighbors(ovk_grid_params *Params, int *NumNeighbors);
void ovkGetGridParamNeighborRanks(ovk_grid_params *Params, int *NeighborRanks);
void ovkSetGridParamNeighborRanks(ovk_grid_params *Params, int NumNeighbors,
  const int *NeighborRanks);

void ovkGetGridPropertyID(const ovk_grid_properties *Properties, int *ID);
void ovkGetGridPropertyDimension(const ovk_grid_properties *Properties, int *NumDims);
void ovkGetGridPropertyGlobalSize(const ovk_grid_properties *Properties, int *GlobalSize);
void ovkGetGridPropertyLocalStart(const ovk_grid_properties *Properties, int *LocalStart);
void ovkGetGridPropertyLocalEnd(const ovk_grid_properties *Properties, int *LocalEnd);
void ovkGetGridPropertyPeriodic(const ovk_grid_properties *Properties, bool *Periodic);
void ovkGetGridPropertyPeriodicStorage(const ovk_grid_properties *Properties,
  ovk_periodic_storage *PeriodicStorage);
void ovkGetGridPropertyPeriodicLength(const ovk_grid_properties *Properties,
  double *PeriodicLength);
void ovkGetGridPropertyGeometryType(const ovk_grid_properties *Properties,
  ovk_geometry_type *GeometryType);
void ovkGetGridPropertyComm(const ovk_grid_properties *Properties, MPI_Comm *Comm);
void ovkGetGridPropertyNumNeighbors(ovk_grid_properties *Properties, int *NumNeighbors);
void ovkGetGridPropertyNeighborRanks(ovk_grid_properties *Properties, int *NeighborRanks);

#endif
