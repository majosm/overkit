// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Grid.h"

#include "ovk/core-c/Cart.h"
#include "ovk/core-c/Constants.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

#include <cstring>
#include <string>

void ovkGetGridID(const ovk_grid *Grid, int *ID) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(ID, "Invalid ID pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridID(GridCPP, *ID);

}

void ovkGetGridName(const ovk_grid *Grid, char *Name) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);

  std::string NameCPP;
  ovk::GetGridName(GridCPP, NameCPP);

  std::strcpy(Name, NameCPP.c_str());

}

void ovkGetGridDimension(const ovk_grid *Grid, int *NumDims) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridDimension(GridCPP, *NumDims);

}

void ovkGetGridComm(const ovk_grid *Grid, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridComm(GridCPP, *Comm);

}

void ovkGetGridCommSize(const ovk_grid *Grid, int *CommSize) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridCommSize(GridCPP, *CommSize);

}

void ovkGetGridCommRank(const ovk_grid *Grid, int *CommRank) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridCommRank(GridCPP, *CommRank);

}

void ovkGetGridCart(const ovk_grid *Grid, ovk_cart *Cart) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Cart, "Invalid cart pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridCart(GridCPP, *Cart);

}

void ovkGetGridSize(const ovk_grid *Grid, int *Size) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridSize(GridCPP, Size);

}

void ovkGetGridPeriodic(const ovk_grid *Grid, bool *Periodic) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridPeriodic(GridCPP, Periodic);

}

void ovkGetGridPeriodicStorage(const ovk_grid *Grid, ovk_periodic_storage *PeriodicStorage) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(PeriodicStorage, "Invalid periodic storage pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);

  ovk::periodic_storage PeriodicStorageCPP;
  ovk::GetGridPeriodicStorage(GridCPP, PeriodicStorageCPP);

  *PeriodicStorage = ovk_periodic_storage(PeriodicStorageCPP);

}

void ovkGetGridPeriodicLength(const ovk_grid *Grid, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridPeriodicLength(GridCPP, PeriodicLength);

}

void ovkGetGridGeometryType(const ovk_grid *Grid, ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);

  ovk::geometry_type GeometryTypeCPP;
  ovk::GetGridGeometryType(GridCPP, GeometryTypeCPP);

  *GeometryType = ovk_geometry_type(GeometryTypeCPP);

}

void ovkGetGridGlobalRange(const ovk_grid *Grid, int *GlobalBegin, int *GlobalEnd) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(GlobalBegin, "Invalid global begin pointer.");
  OVK_DEBUG_ASSERT(GlobalEnd, "Invalid global end pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);

  int NumDims;
  ovk::GetGridDimension(GridCPP, NumDims);

  ovk::range GlobalRange;
  ovk::GetGridGlobalRange(GridCPP, GlobalRange);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    GlobalBegin[iDim] = GlobalRange.Begin(iDim);
    GlobalEnd[iDim] = GlobalRange.End(iDim);
  }

}

void ovkGetGridLocalRange(const ovk_grid *Grid, int *LocalBegin, int *LocalEnd) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(LocalBegin, "Invalid local begin pointer.");
  OVK_DEBUG_ASSERT(LocalEnd, "Invalid local end pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);

  int NumDims;
  ovk::GetGridDimension(GridCPP, NumDims);

  ovk::range LocalRange;
  ovk::GetGridLocalRange(GridCPP, LocalRange);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    LocalBegin[iDim] = LocalRange.Begin(iDim);
    LocalEnd[iDim] = LocalRange.End(iDim);
  }

}

void ovkGetGridGlobalCount(const ovk_grid *Grid, long long *NumGlobal) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(NumGlobal, "Invalid num global pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridGlobalCount(GridCPP, *NumGlobal);

}

void ovkGetGridLocalCount(const ovk_grid *Grid, long long *NumLocal) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(NumLocal, "Invalid num local pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  ovk::GetGridLocalCount(GridCPP, *NumLocal);

}

void ovkCreateGridParams(ovk_grid_params **Params, int NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::grid_params();

  ovk::CreateGridParams(*ParamsCPPPtr, NumDims);

  *Params = reinterpret_cast<ovk_grid_params *>(ParamsCPPPtr);

}

void ovkDestroyGridParams(ovk_grid_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::grid_params *>(*Params);

  ovk::DestroyGridParams(*ParamsCPPPtr);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetGridParamName(const ovk_grid_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);

  std::string NameCPP;
  ovk::GetGridParamName(ParamsCPP, NameCPP);

  std::strcpy(Name, NameCPP.c_str());

}

void ovkSetGridParamName(ovk_grid_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid_params *>(Params);
  ovk::SetGridParamName(ParamsCPP, Name);

}

void ovkGetGridParamDimension(const ovk_grid_params *Params, int *NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);
  ovk::GetGridParamDimension(ParamsCPP, *NumDims);

}

void ovkGetGridParamComm(const ovk_grid_params *Params, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);
  ovk::GetGridParamComm(ParamsCPP, *Comm);

}

void ovkSetGridParamComm(ovk_grid_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid_params *>(Params);
  ovk::SetGridParamComm(ParamsCPP, Comm);

}

void ovkGetGridParamSize(const ovk_grid_params *Params, int *Size) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);
  ovk::GetGridParamSize(ParamsCPP, Size);

}

void ovkSetGridParamSize(ovk_grid_params *Params, const int *Size) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid_params *>(Params);
  ovk::SetGridParamSize(ParamsCPP, Size);

}

void ovkGetGridParamPeriodic(const ovk_grid_params *Params, bool *Periodic) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);
  ovk::GetGridParamPeriodic(ParamsCPP, Periodic);

}

void ovkSetGridParamPeriodic(ovk_grid_params *Params, const bool *Periodic) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid_params *>(Params);
  ovk::SetGridParamPeriodic(ParamsCPP, Periodic);

}

void ovkGetGridParamPeriodicStorage(const ovk_grid_params *Params,
  ovk_periodic_storage *PeriodicStorage) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(PeriodicStorage, "Invalid periodic storage pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);

  ovk::periodic_storage PeriodicStorageCPP;
  ovk::GetGridParamPeriodicStorage(ParamsCPP, PeriodicStorageCPP);

  *PeriodicStorage = ovk_periodic_storage(PeriodicStorageCPP);

}

void ovkSetGridParamPeriodicStorage(ovk_grid_params *Params, ovk_periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid_params *>(Params);
  ovk::SetGridParamPeriodicStorage(ParamsCPP, ovk::periodic_storage(PeriodicStorage));

}

void ovkGetGridParamPeriodicLength(const ovk_grid_params *Params, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);
  ovk::GetGridParamPeriodicLength(ParamsCPP, PeriodicLength);

}

void ovkSetGridParamPeriodicLength(ovk_grid_params *Params, const double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid_params *>(Params);
  ovk::SetGridParamPeriodicLength(ParamsCPP, PeriodicLength);

}

void ovkGetGridParamGeometryType(const ovk_grid_params *Params, ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);

  ovk::geometry_type GeometryTypeCPP;
  ovk::GetGridParamGeometryType(ParamsCPP, GeometryTypeCPP);

  *GeometryType = ovk_geometry_type(GeometryTypeCPP);

}

void ovkSetGridParamGeometryType(ovk_grid_params *Params, ovk_geometry_type GeometryType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid_params *>(Params);
  ovk::SetGridParamGeometryType(ParamsCPP, ovk::geometry_type(GeometryType));

}

void ovkGetGridParamLocalRange(const ovk_grid_params *Params, int *LocalBegin, int *LocalEnd) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LocalBegin, "Invalid local begin pointer.");
  OVK_DEBUG_ASSERT(LocalEnd, "Invalid local end pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);

  int NumDims;
  ovk::GetGridParamDimension(ParamsCPP, NumDims);

  ovk::range LocalRange;
  ovk::GetGridParamLocalRange(ParamsCPP, LocalRange);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    LocalBegin[iDim] = LocalRange.Begin(iDim);
    LocalEnd[iDim] = LocalRange.End(iDim);
  }

}

void ovkSetGridParamLocalRange(ovk_grid_params *Params, const int *LocalBegin, const int *LocalEnd) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LocalBegin, "Invalid local begin pointer.");
  OVK_DEBUG_ASSERT(LocalEnd, "Invalid local end pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid_params *>(Params);

  int NumDims;
  ovk::GetGridParamDimension(ParamsCPP, NumDims);

  ovk::range LocalRange = ovk::MakeEmptyRange(NumDims);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    LocalRange.Begin(iDim) = LocalBegin[iDim];
    LocalRange.End(iDim) = LocalEnd[iDim];
  }

  ovk::SetGridParamLocalRange(ParamsCPP, LocalRange);

}

void ovkGetGridInfoID(const ovk_grid_info *Info, int *ID) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(ID, "Invalid ID pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  ovk::GetGridInfoID(InfoCPP, *ID);
}

void ovkGetGridInfoName(const ovk_grid_info *Info, char *Name) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);

  std::string NameCPP;
  ovk::GetGridInfoName(InfoCPP, NameCPP);

  std::strcpy(Name, NameCPP.c_str());

}

void ovkGetGridInfoDimension(const ovk_grid_info *Info, int *NumDims) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  ovk::GetGridInfoDimension(InfoCPP, *NumDims);

}

void ovkGetGridInfoRootRank(const ovk_grid_info *Info, int *RootRank) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(RootRank, "Invalid root rank pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  ovk::GetGridInfoRootRank(InfoCPP, *RootRank);

}

void ovkGetGridInfoCart(const ovk_grid_info *Info, ovk_cart *Cart) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Cart, "Invalid cart pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  ovk::GetGridInfoCart(InfoCPP, *Cart);

}


void ovkGetGridInfoPeriodicLength(const ovk_grid_info *Info, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  ovk::GetGridInfoPeriodicLength(InfoCPP, PeriodicLength);

}

void ovkGetGridInfoGeometryType(const ovk_grid_info *Info, ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);

  ovk::geometry_type GeometryTypeCPP;
  ovk::GetGridInfoGeometryType(InfoCPP, GeometryTypeCPP);

  *GeometryType = ovk_geometry_type(GeometryTypeCPP);

}

void ovkGetGridInfoGlobalRange(const ovk_grid_info *Info, int *GlobalBegin, int *GlobalEnd) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(GlobalBegin, "Invalid global begin pointer.");
  OVK_DEBUG_ASSERT(GlobalEnd, "Invalid global end pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);

  int NumDims;
  ovk::GetGridInfoDimension(InfoCPP, NumDims);

  ovk::range GlobalRange;
  ovk::GetGridInfoGlobalRange(InfoCPP, GlobalRange);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    GlobalBegin[iDim] = GlobalRange.Begin(iDim);
    GlobalEnd[iDim] = GlobalRange.End(iDim);
  }

}

void ovkGetGridInfoIsLocal(const ovk_grid_info *Info, bool *IsLocal) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(IsLocal, "Invalid is local pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  ovk::GetGridInfoIsLocal(InfoCPP, *IsLocal);

}
