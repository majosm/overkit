// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Grid.h"

#include "ovk/core-c/Constants.h"
#include "ovk/core-c/Context.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <memory>
#include <cstring>
#include <string>

void ovkGetGridContextC(const ovk_grid *Grid, const ovk_context **Context) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *Context = reinterpret_cast<const ovk_context *>(&GridCPP.Context());

}

void ovkGetGridContext(ovk_grid *Grid, ovk_context **Context) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &GridCPP = *reinterpret_cast<ovk::grid *>(Grid);
  *Context = reinterpret_cast<ovk_context *>(&GridCPP.Context());

}

void ovkGetGridSharedContext(ovk_grid *Grid, ovk_shared_context **Context) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &GridCPP = *reinterpret_cast<ovk::grid *>(Grid);
  auto &ContextCPP = GridCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

void ovkGetGridName(const ovk_grid *Grid, char *Name) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  std::strcpy(Name, GridCPP.Name().c_str());

}

void ovkGetGridDimension(const ovk_grid *Grid, int *NumDims) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *NumDims = GridCPP.Dimension();

}

void ovkGetGridComm(const ovk_grid *Grid, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *Comm = GridCPP.Comm();

}

void ovkGetGridCommSize(const ovk_grid *Grid, int *CommSize) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *CommSize = GridCPP.CommSize();

}

void ovkGetGridCommRank(const ovk_grid *Grid, int *CommRank) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *CommRank = GridCPP.CommRank();

}

void ovkGetGridSize(const ovk_grid *Grid, int *Size) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  for (int iDim = 0; iDim < GridCPP.Dimension(); ++iDim) {
    Size[iDim] = GridCPP.Size(iDim);
  }

}

void ovkGetGridPeriodic(const ovk_grid *Grid, bool *Periodic) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  for (int iDim = 0; iDim < GridCPP.Dimension(); ++iDim) {
    Periodic[iDim] = GridCPP.Periodic(iDim);
  }

}

void ovkGetGridPeriodicStorage(const ovk_grid *Grid, ovk_periodic_storage *PeriodicStorage) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(PeriodicStorage, "Invalid periodic storage pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *PeriodicStorage = ovk_periodic_storage(GridCPP.PeriodicStorage());

}

void ovkGetGridPeriodicLength(const ovk_grid *Grid, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  for (int iDim = 0; iDim < GridCPP.Dimension(); ++iDim) {
    PeriodicLength[iDim] = GridCPP.PeriodicLength(iDim);
  }

}

void ovkGetGridGeometryType(const ovk_grid *Grid, ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *GeometryType = ovk_geometry_type(GridCPP.GeometryType());

}

void ovkGetGridGlobalRange(const ovk_grid *Grid, int *GlobalBegin, int *GlobalEnd) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(GlobalBegin, "Invalid global begin pointer.");
  OVK_DEBUG_ASSERT(GlobalEnd, "Invalid global end pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  for (int iDim = 0; iDim < GridCPP.Dimension(); ++iDim) {
    GlobalBegin[iDim] = GridCPP.GlobalRange().Begin(iDim);
    GlobalEnd[iDim] = GridCPP.GlobalRange().End(iDim);
  }

}

void ovkGetGridLocalRange(const ovk_grid *Grid, int *LocalBegin, int *LocalEnd) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(LocalBegin, "Invalid local begin pointer.");
  OVK_DEBUG_ASSERT(LocalEnd, "Invalid local end pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  for (int iDim = 0; iDim < GridCPP.Dimension(); ++iDim) {
    LocalBegin[iDim] = GridCPP.LocalRange().Begin(iDim);
    LocalEnd[iDim] = GridCPP.LocalRange().End(iDim);
  }

}

void ovkGetGridGlobalCount(const ovk_grid *Grid, long long *NumGlobal) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(NumGlobal, "Invalid num global pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *NumGlobal = GridCPP.GlobalRange().Count();

}

void ovkGetGridLocalCount(const ovk_grid *Grid, long long *NumLocal) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(NumLocal, "Invalid num local pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *NumLocal = GridCPP.LocalRange().Count();

}

void ovkCreateGridParams(ovk_grid_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::grid::params();

  *Params = reinterpret_cast<ovk_grid_params *>(ParamsCPPPtr);

}

void ovkDestroyGridParams(ovk_grid_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::grid::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetGridParamName(const ovk_grid_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  std::strcpy(Name, ParamsCPP.Name().c_str());
  
}

void ovkSetGridParamName(ovk_grid_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);
  ParamsCPP.SetName(Name);

}

void ovkGetGridParamDimension(const ovk_grid_params *Params, int *NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  *NumDims = ParamsCPP.Dimension();

}

void ovkSetGridParamDimension(ovk_grid_params *Params, int NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);
  ParamsCPP.SetDimension(NumDims);

}

void ovkGetGridParamComm(const ovk_grid_params *Params, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  *Comm = ParamsCPP.Comm();

}

void ovkSetGridParamComm(ovk_grid_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);
  ParamsCPP.SetComm(Comm);

}

void ovkGetGridParamSize(const ovk_grid_params *Params, int *Size) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid size pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  for (int iDim = 0; iDim < ParamsCPP.Dimension(); ++iDim) {
    Size[iDim] = ParamsCPP.Size()[iDim];
  }

}

void ovkSetGridParamSize(ovk_grid_params *Params, const int *Size) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);

  ovk::tuple<int> SizeCPP = ovk::MakeUniformTuple<int>(ParamsCPP.Dimension(), 1, 1);
  for (int iDim = 0; iDim < ParamsCPP.Dimension(); ++iDim) {
    SizeCPP[iDim] = Size[iDim]; 
  }

  ParamsCPP.SetSize(SizeCPP);

}

void ovkGetGridParamPeriodic(const ovk_grid_params *Params, bool *Periodic) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  for (int iDim = 0; iDim < ParamsCPP.Dimension(); ++iDim) {
    Periodic[iDim] = ParamsCPP.Periodic()[iDim];
  }

}

void ovkSetGridParamPeriodic(ovk_grid_params *Params, const bool *Periodic) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);

  ovk::tuple<bool> PeriodicCPP = ovk::MakeUniformTuple<bool>(ParamsCPP.Dimension(), false);
  for (int iDim = 0; iDim < ParamsCPP.Dimension(); ++iDim) {
    PeriodicCPP[iDim] = Periodic[iDim]; 
  }

  ParamsCPP.SetPeriodic(PeriodicCPP);

}

void ovkGetGridParamPeriodicStorage(const ovk_grid_params *Params, ovk_periodic_storage
  *PeriodicStorage) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(PeriodicStorage, "Invalid periodic storage pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  *PeriodicStorage = ovk_periodic_storage(ParamsCPP.PeriodicStorage());

}

void ovkSetGridParamPeriodicStorage(ovk_grid_params *Params, ovk_periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);
  ParamsCPP.SetPeriodicStorage(ovk::periodic_storage(PeriodicStorage));

}

void ovkGetGridParamPeriodicLength(const ovk_grid_params *Params, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  for (int iDim = 0; iDim < ParamsCPP.Dimension(); ++iDim) {
    PeriodicLength[iDim] = ParamsCPP.PeriodicLength()[iDim];
  }

}

void ovkSetGridParamPeriodicLength(ovk_grid_params *Params, const double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);

  ovk::tuple<double> PeriodicLengthCPP = ovk::MakeUniformTuple<double>(ParamsCPP.Dimension(), 0.);
  for (int iDim = 0; iDim < ParamsCPP.Dimension(); ++iDim) {
    PeriodicLengthCPP[iDim] = PeriodicLength[iDim]; 
  }

  ParamsCPP.SetPeriodicLength(PeriodicLengthCPP);

}

void ovkGetGridParamGeometryType(const ovk_grid_params *Params, ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  *GeometryType = ovk_geometry_type(ParamsCPP.GeometryType());

}

void ovkSetGridParamGeometryType(ovk_grid_params *Params, ovk_geometry_type GeometryType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);
  ParamsCPP.SetGeometryType(ovk::geometry_type(GeometryType));

}

void ovkGetGridParamLocalRange(const ovk_grid_params *Params, int *LocalBegin, int *LocalEnd) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LocalBegin, "Invalid local begin pointer.");
  OVK_DEBUG_ASSERT(LocalEnd, "Invalid local end pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  for (int iDim = 0; iDim < ParamsCPP.Dimension(); ++iDim) {
    LocalBegin[iDim] = ParamsCPP.LocalRange().Begin(iDim);
    LocalEnd[iDim] = ParamsCPP.LocalRange().End(iDim);
  }

}

void ovkSetGridParamLocalRange(ovk_grid_params *Params, const int *LocalBegin, const int *LocalEnd) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LocalBegin, "Invalid local begin pointer.");
  OVK_DEBUG_ASSERT(LocalEnd, "Invalid local end pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);

  ovk::range LocalRange = ovk::MakeEmptyRange(ParamsCPP.Dimension());
  for (int iDim = 0; iDim < ParamsCPP.Dimension(); ++iDim) {
    LocalRange.Begin(iDim) = LocalBegin[iDim];
    LocalRange.End(iDim) = LocalEnd[iDim];
  }

  ParamsCPP.SetLocalRange(LocalRange);

}

void ovkGetGridInfoName(const ovk_grid_info *Info, char *Name) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  std::strcpy(Name, InfoCPP.Name().c_str());
  
}

void ovkGetGridInfoRootRank(const ovk_grid_info *Info, int *RootRank) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(RootRank, "Invalid root rank pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  *RootRank = InfoCPP.RootRank();

}

void ovkGetGridInfoDimension(const ovk_grid_info *Info, int *NumDims) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  *NumDims = InfoCPP.Cart().Dimension();

}

void ovkGetGridInfoSize(const ovk_grid_info *Info, int *Size) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  for (int iDim = 0; iDim < InfoCPP.Cart().Dimension(); ++iDim) {
    Size[iDim] = InfoCPP.Cart().Range().Size(iDim);
  }

}

void ovkGetGridInfoPeriodic(const ovk_grid_info *Info, bool *Periodic) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  for (int iDim = 0; iDim < InfoCPP.Cart().Dimension(); ++iDim) {
    Periodic[iDim] = InfoCPP.Cart().Periodic(iDim);
  }

}

void ovkGetGridInfoPeriodicStorage(const ovk_grid_info *Info, ovk_periodic_storage *PeriodicStorage)
  {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(PeriodicStorage, "Invalid periodic storage pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  *PeriodicStorage = ovk_periodic_storage(InfoCPP.Cart().PeriodicStorage());

}

void ovkGetGridInfoPeriodicLength(const ovk_grid_info *Info, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  for (int iDim = 0; iDim < InfoCPP.Cart().Dimension(); ++iDim) {
    PeriodicLength[iDim] = InfoCPP.PeriodicLength()[iDim];
  }

}

void ovkGetGridInfoGeometryType(const ovk_grid_info *Info, ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  *GeometryType = ovk_geometry_type(InfoCPP.GeometryType());

}

void ovkGetGridInfoIsLocal(const ovk_grid_info *Info, bool *IsLocal) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(IsLocal, "Invalid is local pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  *IsLocal = InfoCPP.IsLocal();

}
