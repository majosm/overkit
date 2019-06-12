// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/ConnectivityN.h"

#include "ovk/core-c/Context.h"
#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core/ConnectivityN.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"

#include <mpi.h>

#include <memory>

void ovkGetConnectivityNContextC(const ovk_connectivity_n *ConnectivityN, const ovk_context
  **Context) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *Context = reinterpret_cast<const ovk_context *>(&ConnectivityNCPP.Context());

}

void ovkGetConnectivityNContext(ovk_connectivity_n *ConnectivityN, ovk_context **Context) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<ovk::connectivity_n *>(ConnectivityN);
  *Context = reinterpret_cast<ovk_context *>(&ConnectivityNCPP.Context());

}

void ovkGetConnectivityNSharedContext(ovk_connectivity_n *ConnectivityN, ovk_shared_context
  **Context) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<ovk::connectivity_n *>(ConnectivityN);
  auto &ContextCPP = ConnectivityNCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

void ovkGetConnectivityNGridID(const ovk_connectivity_n *ConnectivityN, int *GridID) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *GridID = ConnectivityNCPP.GridID();

}

void ovkGetConnectivityNGrid(const ovk_connectivity_n *ConnectivityN, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *Grid = reinterpret_cast<const ovk_grid *>(&ConnectivityNCPP.Grid());

}

void ovkGetConnectivityNSourceGridID(const ovk_connectivity_n *ConnectivityN, int *SourceGridID) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(SourceGridID, "Invalid source grid ID pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *SourceGridID = ConnectivityNCPP.SourceGridID();

}

void ovkGetConnectivityNSourceGridInfo(const ovk_connectivity_n *ConnectivityN, const ovk_grid_info
  **SourceGridInfo) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(SourceGridInfo, "Invalid source grid info pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *SourceGridInfo = reinterpret_cast<const ovk_grid_info *>(&ConnectivityNCPP.SourceGridInfo());

}

void ovkGetConnectivityNDimension(const ovk_connectivity_n *ConnectivityN, int *NumDims) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *NumDims = ConnectivityNCPP.Dimension();

}

void ovkGetConnectivityNComm(const ovk_connectivity_n *ConnectivityN, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *Comm = ConnectivityNCPP.Comm();

}

void ovkGetConnectivityNCommSize(const ovk_connectivity_n *ConnectivityN, int *CommSize) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *CommSize = ConnectivityNCPP.CommSize();

}

void ovkGetConnectivityNCommRank(const ovk_connectivity_n *ConnectivityN, int *CommRank) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *CommRank = ConnectivityNCPP.CommRank();

}

void ovkGetConnectivityNCount(const ovk_connectivity_n *ConnectivityN, long long *Count) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Count, "Invalid count pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *Count = ConnectivityNCPP.Count();

}

void ovkResizeConnectivityN(ovk_connectivity_n *ConnectivityN, long long Count) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<ovk::connectivity_n *>(ConnectivityN);
  ConnectivityNCPP.Resize(Count);

}

void ovkGetConnectivityNPoints(const ovk_connectivity_n *ConnectivityN, int Dimension, const int
  **Points) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *Points = ConnectivityNCPP.Points().Data(Dimension,0);

}

bool ovkEditingConnectivityNPoints(const ovk_connectivity_n *ConnectivityN) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  return ConnectivityNCPP.EditingPoints();

}

void ovkEditConnectivityNPoints(ovk_connectivity_n *ConnectivityN, int Dimension, int **Points) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<ovk::connectivity_n *>(ConnectivityN);

  ovk::edit_handle<ovk::array<int,2>> EditHandle = ConnectivityNCPP.EditPoints();
  auto &PointsCPP = *EditHandle.Release();

  *Points = PointsCPP.Data(Dimension,0);

}

void ovkRestoreConnectivityNPoints(ovk_connectivity_n *ConnectivityN, int Dimension, int **Points) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<ovk::connectivity_n *>(ConnectivityN);
  ConnectivityNCPP.RestorePoints();

  *Points = nullptr;

}

void ovkGetConnectivityNSources(const ovk_connectivity_n *ConnectivityN, int Dimension, const int
  **Sources) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *Sources = ConnectivityNCPP.Sources().Data(Dimension,0);

}

bool ovkEditingConnectivityNSources(const ovk_connectivity_n *ConnectivityN) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  return ConnectivityNCPP.EditingSources();

}

void ovkEditConnectivityNSources(ovk_connectivity_n *ConnectivityN, int Dimension, int **Sources) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<ovk::connectivity_n *>(ConnectivityN);

  ovk::edit_handle<ovk::array<int,2>> EditHandle = ConnectivityNCPP.EditSources();
  auto &SourcesCPP = *EditHandle.Release();

  *Sources = SourcesCPP.Data(Dimension,0);

}

void ovkRestoreConnectivityNSources(ovk_connectivity_n *ConnectivityN, int Dimension, int **Sources)
  {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<ovk::connectivity_n *>(ConnectivityN);
  ConnectivityNCPP.RestoreSources();

  *Sources = nullptr;

}

void ovkGetConnectivityNSourceRanks(const ovk_connectivity_n *ConnectivityN, const int
  **SourceRanks) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  *SourceRanks = ConnectivityNCPP.SourceRanks().Data();

}

bool ovkEditingConnectivityNSourceRanks(const ovk_connectivity_n *ConnectivityN) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<const ovk::connectivity_n *>(ConnectivityN);
  return ConnectivityNCPP.EditingSourceRanks();

}

void ovkEditConnectivityNSourceRanks(ovk_connectivity_n *ConnectivityN, int **SourceRanks) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<ovk::connectivity_n *>(ConnectivityN);

  ovk::edit_handle<ovk::array<int>> EditHandle = ConnectivityNCPP.EditSourceRanks();
  auto &SourceRanksCPP = *EditHandle.Release();

  *SourceRanks = SourceRanksCPP.Data();

}

void ovkReleaseConnectivityNSourceRanks(ovk_connectivity_n *ConnectivityN, int **SourceRanks) {

  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  auto &ConnectivityNCPP = *reinterpret_cast<ovk::connectivity_n *>(ConnectivityN);
  ConnectivityNCPP.RestoreSourceRanks();

  *SourceRanks = nullptr;

}
