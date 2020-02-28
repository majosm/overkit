// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/OverlapM.h"

#include "ovk/core-c/Context.h"
#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/OverlapM.hpp"

#include <mpi.h>

#include <memory>

void ovkGetOverlapMContextC(const ovk_overlap_m *OverlapM, const ovk_context **Context) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *Context = reinterpret_cast<const ovk_context *>(&OverlapMCPP.Context());

}

void ovkGetOverlapMContext(ovk_overlap_m *OverlapM, ovk_context **Context) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);
  *Context = reinterpret_cast<ovk_context *>(&OverlapMCPP.Context());

}

void ovkGetOverlapMSharedContext(ovk_overlap_m *OverlapM, ovk_shared_context **Context) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);
  auto &ContextCPP = OverlapMCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

void ovkGetOverlapMGrid(const ovk_overlap_m *OverlapM, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *Grid = reinterpret_cast<const ovk_grid*>(&OverlapMCPP.Grid());

}

void ovkGetOverlapMDestinationGridInfo(const ovk_overlap_m *OverlapM, const ovk_grid_info
  **DestinationGridInfo) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(DestinationGridInfo, "Invalid destination grid info pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *DestinationGridInfo = reinterpret_cast<const ovk_grid_info *>(
    &OverlapMCPP.DestinationGridInfo());

}

void ovkGetOverlapMDimension(const ovk_overlap_m *OverlapM, int *NumDims) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *NumDims = OverlapMCPP.Dimension();

}

void ovkGetOverlapMComm(const ovk_overlap_m *OverlapM, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *Comm = OverlapMCPP.Comm();

}

void ovkGetOverlapMCommSize(const ovk_overlap_m *OverlapM, int *CommSize) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *CommSize = OverlapMCPP.Comm().Size();

}

void ovkGetOverlapMCommRank(const ovk_overlap_m *OverlapM, int *CommRank) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *CommRank = OverlapMCPP.Comm().Rank();

}

long long ovkGetOverlapMSize(const ovk_overlap_m *OverlapM) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  return OverlapMCPP.Size();

}

void ovkResizeOverlapM(ovk_overlap_m *OverlapM, long long NumCells) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);
  OverlapMCPP.Resize(NumCells);

}

void ovkGetOverlapMCells(const ovk_overlap_m *OverlapM, int Dimension, const int **Cells) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Cells, "Invalid cells pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *Cells = OverlapMCPP.Cells().Data(Dimension,0);

}

bool ovkEditingOverlapMCells(const ovk_overlap_m *OverlapM) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  return OverlapMCPP.EditingCells();

}

void ovkEditOverlapMCells(ovk_overlap_m *OverlapM, int Dimension, int **Cells) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Cells, "Invalid cells pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);

  ovk::edit_handle<ovk::array<int,2>> EditHandle = OverlapMCPP.EditCells();
  auto &CellsCPP = *EditHandle.Release();

  *Cells = CellsCPP.Data(Dimension,0);

}

void ovkRestoreOverlapMCells(ovk_overlap_m *OverlapM, int Dimension, int **Cells) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Cells, "Invalid cells pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);
  OverlapMCPP.RestoreCells();

  *Cells = nullptr;

}

void ovkGetOverlapMCoords(const ovk_overlap_m *OverlapM, int Dimension, const double **Coords) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *Coords = OverlapMCPP.Coords().Data(Dimension,0);

}

bool ovkEditingOverlapMCoords(const ovk_overlap_m *OverlapM) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  return OverlapMCPP.EditingCoords();

}

void ovkEditOverlapMCoords(ovk_overlap_m *OverlapM, int Dimension, double **Coords) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);

  ovk::edit_handle<ovk::array<double,2>> EditHandle = OverlapMCPP.EditCoords();
  auto &CoordsCPP = *EditHandle.Release();

  *Coords = CoordsCPP.Data(Dimension,0);

}

void ovkRestoreOverlapMCoords(ovk_overlap_m *OverlapM, int Dimension, double **Coords) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);
  OverlapMCPP.RestoreCoords();

  *Coords = nullptr;

}

void ovkGetOverlapMDestinations(const ovk_overlap_m *OverlapM, int Dimension, const int
  **Destinations) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *Destinations = OverlapMCPP.Destinations().Data(Dimension,0);

}

bool ovkEditingOverlapMDestinations(const ovk_overlap_m *OverlapM) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  return OverlapMCPP.EditingDestinations();

}

void ovkEditOverlapMDestinations(ovk_overlap_m *OverlapM, int Dimension, int **Destinations) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);

  ovk::edit_handle<ovk::array<int,2>> EditHandle = OverlapMCPP.EditDestinations();
  auto &DestinationsCPP = *EditHandle.Release();

  *Destinations = DestinationsCPP.Data(Dimension,0);

}

void ovkRestoreOverlapMDestinations(ovk_overlap_m *OverlapM, int Dimension, int **Destinations) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < ovk::MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);
  OverlapMCPP.RestoreDestinations();

  *Destinations = nullptr;

}

void ovkGetOverlapMDestinationRanks(const ovk_overlap_m *OverlapM, const int **DestinationRanks) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  *DestinationRanks = OverlapMCPP.DestinationRanks().Data();

}

bool ovkEditingOverlapMDestinationRanks(const ovk_overlap_m *OverlapM) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");

  auto &OverlapMCPP = *reinterpret_cast<const ovk::overlap_m *>(OverlapM);
  return OverlapMCPP.EditingDestinationRanks();

}

void ovkEditOverlapMDestinationRanks(ovk_overlap_m *OverlapM, int **DestinationRanks) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);

  ovk::edit_handle<ovk::array<int>> EditHandle = OverlapMCPP.EditDestinationRanks();
  auto &DestinationRanksCPP = *EditHandle.Release();

  *DestinationRanks = DestinationRanksCPP.Data();

}

void ovkRestoreOverlapMDestinationRanks(ovk_overlap_m *OverlapM, int **DestinationRanks) {

  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  auto &OverlapMCPP = *reinterpret_cast<ovk::overlap_m *>(OverlapM);
  OverlapMCPP.RestoreDestinationRanks();

  *DestinationRanks = nullptr;

}
