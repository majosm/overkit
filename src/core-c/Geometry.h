// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_GEOMETRY_H_INCLUDED
#define OVK_CORE_C_GEOMETRY_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core-c/Geometry.h>
#include <ovk/core-c/Grid.h>
#include <ovk/core/Geometry.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_geometry;
typedef struct ovk_geometry ovk_geometry;

struct ovk_geometry_params;
typedef struct ovk_geometry_params ovk_geometry_params;

void ovkGetGeometryContextC(const ovk_geometry *Geometry, const ovk_context **Context);
void ovkGetGeometryContext(ovk_geometry *Geometry, ovk_context **Context);
void ovkGetGeometrySharedContext(ovk_geometry *Geometry, ovk_shared_context **Context);

void ovkGetGeometryGrid(const ovk_geometry *Geometry, const ovk_grid **Grid);

void ovkGetGeometryDimension(const ovk_geometry *Geometry, int *NumDims);
void ovkGetGeometryComm(const ovk_geometry *Geometry, MPI_Comm *Comm);
void ovkGetGeometryCommSize(const ovk_geometry *Geometry, int *CommSize);
void ovkGetGeometryCommRank(const ovk_geometry *Geometry, int *CommRank);

void ovkGetGeometryType(const ovk_geometry *Geometry, ovk_geometry_type *Type);

void ovkGetGeometryPeriodicLength(const ovk_geometry *Geometry, double *PeriodicLength);
void ovkSetGeometryPeriodicLength(ovk_geometry *Geometry, const double *PeriodicLength);

void ovkGetGeometryCoords(const ovk_geometry *Geometry, int Dimension, const double **Coords);
bool ovkEditingGeometryCoords(const ovk_geometry *Geometry);
void ovkEditGeometryCoords(ovk_geometry *Geometry, int Dimension, double **Coords);
void ovkRestoreGeometryCoords(ovk_geometry *Geometry, int Dimension, double **Coords);

void ovkCreateGeometryParams(ovk_geometry_params **Params);
void ovkDestroyGeometryParams(ovk_geometry_params **Params);

void ovkGetGeometryParamType(const ovk_geometry_params *Params, ovk_geometry_type *Type);
void ovkSetGeometryParamType(ovk_geometry_params *Params, ovk_geometry_type Type);
void ovkGetGeometryParamPeriodicLength(const ovk_geometry_params *Params, double *PeriodicLength);
void ovkSetGeometryParamPeriodicLength(ovk_geometry_params *Params, const double *PeriodicLength);

#ifdef __cplusplus
}
#endif

#endif
