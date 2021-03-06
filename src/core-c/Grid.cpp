// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Grid.h"

#include "ovk/core-c/Context.h"
#include "ovk/core-c/Global.h"
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
  *CommSize = GridCPP.Comm().Size();

}

void ovkGetGridCommRank(const ovk_grid *Grid, int *CommRank) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *CommRank = GridCPP.Comm().Rank();

}

void ovkGetGridGlobalRange(const ovk_grid *Grid, int *GlobalBegin, int *GlobalEnd) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(GlobalBegin, "Invalid global begin pointer.");
  OVK_DEBUG_ASSERT(GlobalEnd, "Invalid global end pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    GlobalBegin[iDim] = GridCPP.GlobalRange().Begin(iDim);
    GlobalEnd[iDim] = GridCPP.GlobalRange().End(iDim);
  }

}

void ovkGetGridLocalRange(const ovk_grid *Grid, int *LocalBegin, int *LocalEnd) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(LocalBegin, "Invalid local begin pointer.");
  OVK_DEBUG_ASSERT(LocalEnd, "Invalid local end pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    LocalBegin[iDim] = GridCPP.LocalRange().Begin(iDim);
    LocalEnd[iDim] = GridCPP.LocalRange().End(iDim);
  }

}

void ovkGetGridExtendedRange(const ovk_grid *Grid, int *ExtendedBegin, int *ExtendedEnd) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(ExtendedBegin, "Invalid extended begin pointer.");
  OVK_DEBUG_ASSERT(ExtendedEnd, "Invalid extended end pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    ExtendedBegin[iDim] = GridCPP.ExtendedRange().Begin(iDim);
    ExtendedEnd[iDim] = GridCPP.ExtendedRange().End(iDim);
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

void ovkGetGridExtendedCount(const ovk_grid *Grid, long long *NumExtended) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(NumExtended, "Invalid num extended pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *NumExtended = GridCPP.ExtendedRange().Count();

}

void ovkGetGridPeriodic(const ovk_grid *Grid, bool *Periodic) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    Periodic[iDim] = GridCPP.Periodic(iDim);
  }

}

void ovkGetGridPeriodicStorage(const ovk_grid *Grid, ovk_periodic_storage *PeriodicStorage) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(PeriodicStorage, "Invalid periodic storage pointer.");

  auto &GridCPP = *reinterpret_cast<const ovk::grid *>(Grid);
  *PeriodicStorage = ovk_periodic_storage(GridCPP.PeriodicStorage());

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

void ovkGetGridParamGlobalRange(const ovk_grid_params *Params, int *GlobalBegin, int *GlobalEnd) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(GlobalBegin, "Invalid global begin pointer.");
  OVK_DEBUG_ASSERT(GlobalEnd, "Invalid global end pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    GlobalBegin[iDim] = ParamsCPP.GlobalRange().Begin(iDim);
    GlobalEnd[iDim] = ParamsCPP.GlobalRange().End(iDim);
  }

}

void ovkSetGridParamGlobalRange(ovk_grid_params *Params, const int *GlobalBegin, const int
  *GlobalEnd) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(GlobalBegin, "Invalid global begin pointer.");
  OVK_DEBUG_ASSERT(GlobalEnd, "Invalid global end pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);
  ParamsCPP.SetGlobalRange({GlobalBegin, GlobalEnd});

}

void ovkGetGridParamLocalRange(const ovk_grid_params *Params, int *LocalBegin, int *LocalEnd) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LocalBegin, "Invalid local begin pointer.");
  OVK_DEBUG_ASSERT(LocalEnd, "Invalid local end pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    LocalBegin[iDim] = ParamsCPP.LocalRange().Begin(iDim);
    LocalEnd[iDim] = ParamsCPP.LocalRange().End(iDim);
  }

}

void ovkSetGridParamLocalRange(ovk_grid_params *Params, const int *LocalBegin, const int *LocalEnd) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LocalBegin, "Invalid local begin pointer.");
  OVK_DEBUG_ASSERT(LocalEnd, "Invalid local end pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);
  ParamsCPP.SetLocalRange({LocalBegin, LocalEnd});

}

void ovkGetGridParamPeriodic(const ovk_grid_params *Params, bool *Periodic) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::grid::params *>(Params);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    Periodic[iDim] = ParamsCPP.Periodic()[iDim];
  }

}

void ovkSetGridParamPeriodic(ovk_grid_params *Params, const bool *Periodic) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::grid::params *>(Params);
  ParamsCPP.SetPeriodic(Periodic);

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

void ovkGetGridInfoGlobalRange(const ovk_grid_info *Info, int *GlobalBegin, int *GlobalEnd) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(GlobalBegin, "Invalid global begin pointer.");
  OVK_DEBUG_ASSERT(GlobalEnd, "Invalid global end pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    GlobalBegin[iDim] = InfoCPP.Cart().Range().Begin(iDim);
    GlobalEnd[iDim] = InfoCPP.Cart().Range().End(iDim);
  }

}

void ovkGetGridInfoSize(const ovk_grid_info *Info, int *Size) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
    Size[iDim] = InfoCPP.Cart().Range().Size(iDim);
  }

}

void ovkGetGridInfoPeriodic(const ovk_grid_info *Info, bool *Periodic) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
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

void ovkGetGridInfoIsLocal(const ovk_grid_info *Info, bool *IsLocal) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(IsLocal, "Invalid is local pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::grid_info *>(Info);
  *IsLocal = InfoCPP.IsLocal();

}
