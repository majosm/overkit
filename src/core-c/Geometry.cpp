// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Geometry.h"

#include "ovk/core-c/Context.h"
#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Geometry.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"

#include <mpi.h>

#include <memory>

void ovkGetGeometryContextC(const ovk_geometry *Geometry, const ovk_context **Context) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  *Context = reinterpret_cast<const ovk_context *>(&GeometryCPP.Context());

}

void ovkGetGeometryContext(ovk_geometry *Geometry, ovk_context **Context) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &GeometryCPP = *reinterpret_cast<ovk::geometry *>(Geometry);
  *Context = reinterpret_cast<ovk_context *>(&GeometryCPP.Context());

}

void ovkGetGeometrySharedContext(ovk_geometry *Geometry, ovk_shared_context **Context) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &GeometryCPP = *reinterpret_cast<ovk::geometry *>(Geometry);
  auto &ContextCPP = GeometryCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

void ovkGetGeometryGrid(const ovk_geometry *Geometry, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  *Grid = reinterpret_cast<const ovk_grid *>(&GeometryCPP.Grid());

}

void ovkGetGeometryDimension(const ovk_geometry *Geometry, int *NumDims) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  *NumDims = GeometryCPP.Dimension();

}

void ovkGetGeometryComm(const ovk_geometry *Geometry, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  *Comm = GeometryCPP.Comm();

}

void ovkGetGeometryCommSize(const ovk_geometry *Geometry, int *CommSize) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  *CommSize = GeometryCPP.Comm().Size();

}

void ovkGetGeometryCommRank(const ovk_geometry *Geometry, int *CommRank) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  *CommRank = GeometryCPP.Comm().Rank();

}

void ovkGetGeometryType(const ovk_geometry *Geometry, ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  *GeometryType = ovk_geometry_type(GeometryCPP.GeometryType());

}

void ovkGetGeometryPeriodicLength(const ovk_geometry *Geometry, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = GeometryCPP.PeriodicLength(iDim);
  }

}

void ovkSetGeometryPeriodicLength(ovk_geometry *Geometry, const double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  auto &GeometryCPP = *reinterpret_cast<ovk::geometry *>(Geometry);
  GeometryCPP.SetPeriodicLength(PeriodicLength);

}

void ovkGetGeometryCoords(const ovk_geometry *Geometry, int Dimension, const double **Coords) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  *Coords = GeometryCPP.Coords()(Dimension).Data();

}

bool ovkEditingGeometryCoords(const ovk_geometry *Geometry) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");

  auto &GeometryCPP = *reinterpret_cast<const ovk::geometry *>(Geometry);
  return GeometryCPP.EditingCoords();

}

void ovkEditGeometryCoords(ovk_geometry *Geometry, int Dimension, double **Coords) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &GeometryCPP = *reinterpret_cast<ovk::geometry *>(Geometry);

  ovk::edit_handle<ovk::array<ovk::field<double>>> EditHandle = GeometryCPP.EditCoords();
  auto &CoordsCPP = *EditHandle.Release();

  *Coords = CoordsCPP(Dimension).Data();

}

void ovkRestoreGeometryCoords(ovk_geometry *Geometry, int Dimension, double **Coords) {

  OVK_DEBUG_ASSERT(Geometry, "Invalid geometry pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &GeometryCPP = *reinterpret_cast<ovk::geometry *>(Geometry);
  GeometryCPP.RestoreCoords();

  *Coords = nullptr;

}

void ovkGetGeometryParamGeometryType(const ovk_geometry_params *Params, ovk_geometry_type
  *GeometryType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::geometry::params *>(Params);
  *GeometryType = ovk_geometry_type(ParamsCPP.GeometryType());

}

void ovkSetGeometryParamGeometryType(ovk_geometry_params *Params, ovk_geometry_type GeometryType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::geometry::params *>(Params);
  ParamsCPP.SetGeometryType(ovk::geometry_type(GeometryType));

}

void ovkGetGeometryParamPeriodicLength(const ovk_geometry_params *Params, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::geometry::params *>(Params);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = ParamsCPP.PeriodicLength()[iDim];
  }

}

void ovkSetGeometryParamPeriodicLength(ovk_geometry_params *Params, const double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::geometry::params *>(Params);
  ParamsCPP.SetPeriodicLength(PeriodicLength);

}
