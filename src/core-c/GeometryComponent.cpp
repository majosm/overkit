// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/GeometryComponent.h"

#include "ovk/core-c/Geometry.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Geometry.hpp"
#include "ovk/core/GeometryComponent.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <utility>

extern "C" {

int ovkGeometryCount(const ovk_geometry_component *GeometryComponent) {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");

  auto &GeometryComponentCPP = *reinterpret_cast<const ovk::geometry_component *>(
    GeometryComponent);
  return GeometryComponentCPP.GeometryCount();

}

bool ovkGeometryExists(const ovk_geometry_component *GeometryComponent, int GridID) {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");

  auto &GeometryComponentCPP = *reinterpret_cast<const ovk::geometry_component *>(
    GeometryComponent);
  return GeometryComponentCPP.GeometryExists(GridID);

}

void ovkCreateGeometry(ovk_geometry_component *GeometryComponent, int GridID, ovk_geometry_params
  **MaybeParams) {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");

  auto &GeometryComponentCPP = *reinterpret_cast<ovk::geometry_component *>(
    GeometryComponent);

  ovk::geometry::params *ParamsCPPPtr = nullptr;
  if (MaybeParams && *MaybeParams) {
    ParamsCPPPtr = reinterpret_cast<ovk::geometry::params *>(*MaybeParams);
  }

  ovk::optional<ovk::geometry::params> MaybeParamsCPP;
  if (ParamsCPPPtr) {
    MaybeParamsCPP = std::move(*ParamsCPPPtr);
  }

  GeometryComponentCPP.CreateGeometry(GridID, std::move(MaybeParamsCPP));

  if (ParamsCPPPtr) {
    delete ParamsCPPPtr;
    *MaybeParams = nullptr;
  }

}

void ovkCreateGeometries(ovk_geometry_component *GeometryComponent, int Count, const int *GridIDs,
  ovk_geometry_params **MaybeParams) {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");

  auto &GeometryComponentCPP = *reinterpret_cast<ovk::geometry_component *>(
    GeometryComponent);

  if (MaybeParams) {

    ovk::array<ovk::geometry::params *> ParamsCPPPtrs({Count}, nullptr);
    for (int iCreate = 0; iCreate < Count; ++iCreate) {
      if (MaybeParams[iCreate]) {
        ParamsCPPPtrs(iCreate) = reinterpret_cast<ovk::geometry::params *>(MaybeParams[iCreate]);
      }
    }

    ovk::array<ovk::optional<ovk::geometry::params>> MaybeParamsCPP({Count});
    for (int iCreate = 0; iCreate < Count; ++iCreate) {
      if (ParamsCPPPtrs(iCreate)) {
        MaybeParamsCPP(iCreate) = std::move(*ParamsCPPPtrs(iCreate));
      }
    }

    GeometryComponentCPP.CreateGeometries({GridIDs, {Count}}, std::move(MaybeParamsCPP));

    for (int iCreate = 0; iCreate < Count; ++iCreate) {
      if (ParamsCPPPtrs(iCreate)) {
        delete ParamsCPPPtrs(iCreate);
        MaybeParams[iCreate] = nullptr;
      }
    }

  } else {

    GeometryComponentCPP.CreateGeometries({GridIDs, {Count}});

  }

}

void ovkDestroyGeometry(ovk_geometry_component *GeometryComponent, int GridID) {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");

  auto &GeometryComponentCPP = *reinterpret_cast<ovk::geometry_component *>(
    GeometryComponent);
  GeometryComponentCPP.DestroyGeometry(GridID);

}

void ovkDestroyGeometries(ovk_geometry_component *GeometryComponent, int Count, const int *GridIDs)
  {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");

  auto &GeometryComponentCPP = *reinterpret_cast<ovk::geometry_component *>(
    GeometryComponent);
  GeometryComponentCPP.DestroyGeometries({GridIDs, {Count}});

}

int ovkLocalGeometryCount(const ovk_geometry_component *GeometryComponent) {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");

  auto &GeometryComponentCPP = *reinterpret_cast<const ovk::geometry_component *>(
    GeometryComponent);
  return GeometryComponentCPP.LocalGeometryCount();

}

void ovkGetGeometry(const ovk_geometry_component *GeometryComponent, int GridID, const ovk_geometry
  **Geometry) {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");
  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");

  auto &GeometryComponentCPP = *reinterpret_cast<const ovk::geometry_component *>(
    GeometryComponent);
  *Geometry = reinterpret_cast<const ovk_geometry *>(&GeometryComponentCPP.Geometry(GridID));

}

bool ovkEditingGeometry(const ovk_geometry_component *GeometryComponent, int GridID) {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");

  auto &GeometryComponentCPP = *reinterpret_cast<const ovk::geometry_component *>(
    GeometryComponent);
  return GeometryComponentCPP.EditingGeometry(GridID);

}

void ovkEditGeometry(ovk_geometry_component *GeometryComponent, int GridID, ovk_geometry **Geometry)
  {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");
  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");

  auto &GeometryComponentCPP = *reinterpret_cast<ovk::geometry_component *>(
    GeometryComponent);

  ovk::edit_handle<ovk::geometry> EditHandle = GeometryComponentCPP.EditGeometry(GridID);
  auto GeometryCPPPtr = EditHandle.Release();

  *Geometry = reinterpret_cast<ovk_geometry *>(GeometryCPPPtr);

}

void ovkRestoreGeometry(ovk_geometry_component *GeometryComponent, int GridID, ovk_geometry
  **Geometry) {

  OVK_DEBUG_ASSERT(GeometryComponent, "Invalid geometry component pointer.");
  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(*Geometry, "Invalid geometry pointer.");

  auto &GeometryComponentCPP = *reinterpret_cast<ovk::geometry_component *>(
    GeometryComponent);
  GeometryComponentCPP.RestoreGeometry(GridID);

  *Geometry = nullptr;

}

void ovkCreateGeometryComponentParams(ovk_geometry_component_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::geometry_component::params();

  *Params = reinterpret_cast<ovk_geometry_component_params *>(ParamsCPPPtr);

}

void ovkDestroyGeometryComponentParams(ovk_geometry_component_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::geometry_component::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetGeometryComponentParamName(const ovk_geometry_component_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::geometry_component::params *>(Params);
  std::strcpy(Name, ParamsCPP.Name().c_str());

}

void ovkSetGeometryComponentParamName(ovk_geometry_component_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::geometry_component::params *>(Params);
  ParamsCPP.SetName(Name);

}

}
