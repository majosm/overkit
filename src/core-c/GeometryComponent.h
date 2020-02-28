// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_GEOMETRY_COMPONENT_H_INCLUDED
#define OVK_CORE_C_GEOMETRY_COMPONENT_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Geometry.h>
#include <ovk/core-c/Global.h>
#include <ovk/core/GeometryComponent.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_geometry_component;
typedef struct ovk_geometry_component ovk_geometry_component;

struct ovk_geometry_component_params;
typedef struct ovk_geometry_component_params ovk_geometry_component_params;

int ovkGeometryCount(const ovk_geometry_component *GeometryComponent);

bool ovkGeometryExists(const ovk_geometry_component *GeometryComponent, int GridID);

void ovkCreateGeometry(ovk_geometry_component *GeometryComponent, int GridID, ovk_geometry_params
  **MaybeParams);
void ovkCreateGeometries(ovk_geometry_component *GeometryComponent, int Count, const int *GridIDs,
  ovk_geometry_params **MaybeParams);

void ovkDestroyGeometry(ovk_geometry_component *GeometryComponent, int GridID);
void ovkDestroyGeometries(ovk_geometry_component *GeometryComponent, int Count, const int *GridIDs);

int LocalGeometryCount(const ovk_geometry_component *GeometryComponent);

void ovkGetGeometry(const ovk_geometry_component *GeometryComponent, int GridID, const ovk_geometry
  **Geometry);
bool ovkEditingGeometry(const ovk_geometry_component *GeometryComponent, int GridID);
void ovkEditGeometry(ovk_geometry_component *GeometryComponent, int GridID, ovk_geometry
  **Geometry);
void ovkRestoreGeometry(ovk_geometry_component *GeometryComponent, int GridID, ovk_geometry
  **Geometry);

void ovkCreateGeometryComponentParams(ovk_geometry_component_params **Params);
void ovkDestroyGeometryComponentParams(ovk_geometry_component_params **Params);
void ovkGetGeometryComponentParamName(const ovk_geometry_component_params *Params, char *Name);
void ovkSetGeometryComponentParamName(ovk_geometry_component_params *Params, const char *Name);

#ifdef __cplusplus
}
#endif

#endif
