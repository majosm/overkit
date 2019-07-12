// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/OverlapN.h"

#include "ovk/core-c/Context.h"
#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/OverlapN.hpp"

#include <mpi.h>

#include <memory>

void ovkGetOverlapNContextC(const ovk_overlap_n *OverlapN, const ovk_context **Context) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *Context = reinterpret_cast<const ovk_context *>(&OverlapNCPP.Context());

}

void ovkGetOverlapNContext(ovk_overlap_n *OverlapN, ovk_context **Context) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &OverlapNCPP = *reinterpret_cast<ovk::overlap_n *>(OverlapN);
  *Context = reinterpret_cast<ovk_context *>(&OverlapNCPP.Context());

}

void ovkGetOverlapNSharedContext(ovk_overlap_n *OverlapN, ovk_shared_context **Context) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &OverlapNCPP = *reinterpret_cast<ovk::overlap_n *>(OverlapN);
  auto &ContextCPP = OverlapNCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

void ovkGetOverlapNGrid(const ovk_overlap_n *OverlapN, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *Grid = reinterpret_cast<const ovk_grid *>(&OverlapNCPP.Grid());

}

void ovkGetOverlapNSourceGridInfo(const ovk_overlap_n *OverlapN, const ovk_grid_info
  **SourceGridInfo) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(SourceGridInfo, "Invalid source grid info pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *SourceGridInfo = reinterpret_cast<const ovk_grid_info *>(&OverlapNCPP.SourceGridInfo());

}

void ovkGetOverlapNDimension(const ovk_overlap_n *OverlapN, int *NumDims) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *NumDims = OverlapNCPP.Dimension();

}

void ovkGetOverlapNComm(const ovk_overlap_n *OverlapN, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *Comm = OverlapNCPP.Comm();

}

void ovkGetOverlapNCommSize(const ovk_overlap_n *OverlapN, int *CommSize) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *CommSize = OverlapNCPP.Comm().Size();

}

void ovkGetOverlapNCommRank(const ovk_overlap_n *OverlapN, int *CommRank) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *CommRank = OverlapNCPP.Comm().Rank();

}

void ovkGetOverlapNCount(const ovk_overlap_n *OverlapN, long long *Count) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Count, "Invalid count pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *Count = OverlapNCPP.Count();

}

void ovkResizeOverlapN(ovk_overlap_n *OverlapN, long long Count) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");

  auto &OverlapNCPP = *reinterpret_cast<ovk::overlap_n *>(OverlapN);
  OverlapNCPP.Resize(Count);

}

void ovkGetOverlapNMask(const ovk_overlap_n *OverlapN, const bool **OverlapMask) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(OverlapMask, "Invalid overlap mask pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *OverlapMask = OverlapNCPP.Mask().Data();

}

void ovkGetOverlapNPoints(const ovk_overlap_n *OverlapN, int Dimension, const int **Points) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *Points = OverlapNCPP.Points().Data(Dimension,0);

}

bool ovkEditingOverlapNPoints(const ovk_overlap_n *OverlapN) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  return OverlapNCPP.EditingPoints();

}

void ovkEditOverlapNPoints(ovk_overlap_n *OverlapN, int Dimension, int **Points) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  auto &OverlapNCPP = *reinterpret_cast<ovk::overlap_n *>(OverlapN);

  ovk::edit_handle<ovk::array<int,2>> EditHandle = OverlapNCPP.EditPoints();
  auto &PointsCPP = *EditHandle.Release();

  *Points = PointsCPP.Data(Dimension,0);

}

void ovkRestoreOverlapNPoints(ovk_overlap_n *OverlapN, int Dimension, int **Points) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  auto &OverlapNCPP = *reinterpret_cast<ovk::overlap_n *>(OverlapN);
  OverlapNCPP.RestorePoints();

  *Points = nullptr;

}

void ovkGetOverlapNSources(const ovk_overlap_n *OverlapN, int Dimension, const int **Sources) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *Sources = OverlapNCPP.Sources().Data(Dimension,0);

}

bool ovkEditingOverlapNSources(const ovk_overlap_n *OverlapN) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  return OverlapNCPP.EditingSources();

}

void ovkEditOverlapNSources(ovk_overlap_n *OverlapN, int Dimension, int **Sources) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  auto &OverlapNCPP = *reinterpret_cast<ovk::overlap_n *>(OverlapN);

  ovk::edit_handle<ovk::array<int,2>> EditHandle = OverlapNCPP.EditSources();
  auto &SourcesCPP = *EditHandle.Release();

  *Sources = SourcesCPP.Data(Dimension,0);

}

void ovkRestoreOverlapNSources(ovk_overlap_n *OverlapN, int Dimension, int **Sources) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  auto &OverlapNCPP = *reinterpret_cast<ovk::overlap_n *>(OverlapN);
  OverlapNCPP.RestoreSources();

  *Sources = nullptr;

}

void ovkGetOverlapNSourceRanks(const ovk_overlap_n *OverlapN, const int **SourceRanks) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  *SourceRanks = OverlapNCPP.SourceRanks().Data();

}

bool ovkEditingOverlapNSourceRanks(const ovk_overlap_n *OverlapN) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");

  auto &OverlapNCPP = *reinterpret_cast<const ovk::overlap_n *>(OverlapN);
  return OverlapNCPP.EditingSourceRanks();

}

void ovkEditOverlapNSourceRanks(ovk_overlap_n *OverlapN, int **SourceRanks) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  auto &OverlapNCPP = *reinterpret_cast<ovk::overlap_n *>(OverlapN);

  ovk::edit_handle<ovk::array<int>> EditHandle = OverlapNCPP.EditSourceRanks();
  auto &SourceRanksCPP = *EditHandle.Release();

  *SourceRanks = SourceRanksCPP.Data();

}

void ovkReleaseOverlapNSourceRanks(ovk_overlap_n *OverlapN, int **SourceRanks) {

  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  auto &OverlapNCPP = *reinterpret_cast<ovk::overlap_n *>(OverlapN);
  OverlapNCPP.RestoreSourceRanks();

  *SourceRanks = nullptr;

}
