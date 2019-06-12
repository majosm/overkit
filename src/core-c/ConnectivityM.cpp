// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/ConnectivityM.h"

#include "ovk/core-c/Context.h"
#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core/ConnectivityM.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"

#include <mpi.h>

#include <memory>

void ovkGetConnectivityMContextC(const ovk_connectivity_m *ConnectivityM, const ovk_context
  **Context) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *Context = reinterpret_cast<const ovk_context *>(&ConnectivityMCPP.Context());

}

void ovkGetConnectivityMContext(ovk_connectivity_m *ConnectivityM, ovk_context **Context) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);
  *Context = reinterpret_cast<ovk_context *>(&ConnectivityMCPP.Context());

}

void ovkGetConnectivityMSharedContext(ovk_connectivity_m *ConnectivityM, ovk_shared_context
  **Context) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);
  auto &ContextCPP = ConnectivityMCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

void ovkGetConnectivityMGridID(const ovk_connectivity_m *ConnectivityM, int *GridID) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *GridID = ConnectivityMCPP.GridID();

}

void ovkGetConnectivityMGrid(const ovk_connectivity_m *ConnectivityM, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *Grid = reinterpret_cast<const ovk_grid*>(&ConnectivityMCPP.Grid());

}

void ovkGetConnectivityMDestinationGridID(const ovk_connectivity_m *ConnectivityM, int
  *DestinationGridID) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(DestinationGridID, "Invalid destination grid ID pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *DestinationGridID = ConnectivityMCPP.DestinationGridID();

}

void ovkGetConnectivityMDestinationGridInfo(const ovk_connectivity_m *ConnectivityM, const
  ovk_grid_info **DestinationGridInfo) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(DestinationGridInfo, "Invalid destination grid info pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *DestinationGridInfo = reinterpret_cast<const ovk_grid_info *>(
    &ConnectivityMCPP.DestinationGridInfo());

}

void ovkGetConnectivityMDimension(const ovk_connectivity_m *ConnectivityM, int *NumDims) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *NumDims = ConnectivityMCPP.Dimension();

}

void ovkGetConnectivityMComm(const ovk_connectivity_m *ConnectivityM, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *Comm = ConnectivityMCPP.Comm();

}

void ovkGetConnectivityMCommSize(const ovk_connectivity_m *ConnectivityM, int *CommSize) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *CommSize = ConnectivityMCPP.CommSize();

}

void ovkGetConnectivityMCommRank(const ovk_connectivity_m *ConnectivityM, int *CommRank) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *CommRank = ConnectivityMCPP.CommRank();

}

void ovkGetConnectivityMCount(const ovk_connectivity_m *ConnectivityM, long long *Count) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Count, "Invalid count pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *Count = ConnectivityMCPP.Count();

}

void ovkGetConnectivityMMaxSize(const ovk_connectivity_m *ConnectivityM, int *MaxSize) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(MaxSize, "Invalid max size pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *MaxSize = ConnectivityMCPP.MaxSize();

}

void ovkResizeConnectivityM(ovk_connectivity_m *ConnectivityM, long long Count, int MaxSize) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);
  ConnectivityMCPP.Resize(Count, MaxSize);

}

void ovkGetConnectivityMExtents(const ovk_connectivity_m *ConnectivityM, int Dimension, const int
  **Begins, const int **Ends) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Begins, "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(Ends, "Invalid ends pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *Begins = ConnectivityMCPP.Extents().Data(0,Dimension,0);
  *Ends = ConnectivityMCPP.Extents().Data(1,Dimension,0);

}

bool ovkEditingConnectivityMExtents(const ovk_connectivity_m *ConnectivityM) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  return ConnectivityMCPP.EditingExtents();

}

void ovkEditConnectivityMExtents(ovk_connectivity_m *ConnectivityM, int Dimension, int **Begins,
  int **Ends) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Begins, "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(Ends, "Invalid ends pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);

  ovk::edit_handle<ovk::array<int,3>> EditHandle = ConnectivityMCPP.EditExtents();
  auto &ExtentsCPP = *EditHandle.Release();

  *Begins = ExtentsCPP.Data(0,Dimension,0);
  *Ends = ExtentsCPP.Data(1,Dimension,0);

}

void ovkRestoreConnectivityMExtents(ovk_connectivity_m *ConnectivityM, int Dimension, int **Begins,
  int **Ends) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Begins, "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(Ends, "Invalid ends pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);
  ConnectivityMCPP.RestoreExtents();

  *Begins = nullptr;
  *Ends = nullptr;

}

void ovkGetConnectivityMCoords(const ovk_connectivity_m *ConnectivityM, int Dimension, const
  double **Coords) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *Coords = ConnectivityMCPP.Coords().Data(Dimension,0);

}

bool ovkEditingConnectivityMCoords(const ovk_connectivity_m *ConnectivityM) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  return ConnectivityMCPP.EditingCoords();

}

void ovkEditConnectivityMCoords(ovk_connectivity_m *ConnectivityM, int Dimension, double **Coords) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);

  ovk::edit_handle<ovk::array<double,2>> EditHandle = ConnectivityMCPP.EditCoords();
  auto &CoordsCPP = *EditHandle.Release();

  *Coords = CoordsCPP.Data(Dimension,0);

}

void ovkRestoreConnectivityMCoords(ovk_connectivity_m *ConnectivityM, int Dimension, double
  **Coords) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);
  ConnectivityMCPP.RestoreCoords();

  *Coords = nullptr;

}

void ovkGetConnectivityMInterpCoefs(const ovk_connectivity_m *ConnectivityM, int Dimension, int
  Point, const double **InterpCoefs) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(InterpCoefs, "Invalid interp coefs pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);

  OVK_DEBUG_ASSERT(Point >= 0 && Point < ConnectivityMCPP.MaxSize(), "Invalid point.");

  *InterpCoefs = ConnectivityMCPP.InterpCoefs().Data(Dimension,Point,0);

}

bool ovkEditingConnectivityMInterpCoefs(const ovk_connectivity_m *ConnectivityM) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  return ConnectivityMCPP.EditingInterpCoefs();

}

void ovkEditConnectivityMInterpCoefs(ovk_connectivity_m *ConnectivityM, int Dimension, int Point,
  double **InterpCoefs) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(InterpCoefs, "Invalid interp coefs pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);

  OVK_DEBUG_ASSERT(Point >= 0 && Point < ConnectivityMCPP.MaxSize(), "Invalid point.");

  ovk::edit_handle<ovk::array<double,3>> EditHandle = ConnectivityMCPP.EditInterpCoefs();
  auto &InterpCoefsCPP = *EditHandle.Release();

  *InterpCoefs = InterpCoefsCPP.Data(Dimension,Point,0);

}

void ovkRestoreConnectivityMInterpCoefs(ovk_connectivity_m *ConnectivityM, int Dimension, int Point,
  double **InterpCoefs) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(InterpCoefs, "Invalid interp coefs pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);

  OVK_DEBUG_ASSERT(Point >= 0 && Point < ConnectivityMCPP.MaxSize(), "Invalid point.");

  ConnectivityMCPP.RestoreInterpCoefs();

  *InterpCoefs = nullptr;

}

void ovkGetConnectivityMDestinations(const ovk_connectivity_m *ConnectivityM, int Dimension, const
  int **Destinations) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *Destinations = ConnectivityMCPP.Destinations().Data(Dimension,0);

}

bool ovkEditingConnectivityMDestinations(const ovk_connectivity_m *ConnectivityM) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  return ConnectivityMCPP.EditingDestinations();

}

void ovkEditConnectivityMDestinations(ovk_connectivity_m *ConnectivityM, int Dimension, int
  **Destinations) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);

  ovk::edit_handle<ovk::array<int,2>> EditHandle = ConnectivityMCPP.EditDestinations();
  auto &DestinationsCPP = *EditHandle.Release();

  *Destinations = DestinationsCPP.Data(Dimension,0);

}

void ovkRestoreConnectivityMDestinations(ovk_connectivity_m *ConnectivityM, int Dimension, int
  **Destinations) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < OVK_MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);
  ConnectivityMCPP.RestoreDestinations();

  *Destinations = nullptr;

}

void ovkGetConnectivityMDestinationRanks(const ovk_connectivity_m *ConnectivityM, const int
  **DestinationRanks) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  *DestinationRanks = ConnectivityMCPP.DestinationRanks().Data();

}

bool ovkEditingConnectivityMDestinationRanks(const ovk_connectivity_m *ConnectivityM) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<const ovk::connectivity_m *>(ConnectivityM);
  return ConnectivityMCPP.EditingDestinationRanks();

}

void ovkEditConnectivityMDestinationRanks(ovk_connectivity_m *ConnectivityM, int
  **DestinationRanks) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);

  ovk::edit_handle<ovk::array<int>> EditHandle = ConnectivityMCPP.EditDestinationRanks();
  auto &DestinationRanksCPP = *EditHandle.Release();

  *DestinationRanks = DestinationRanksCPP.Data();

}

void ovkRestoreConnectivityMDestinationRanks(ovk_connectivity_m *ConnectivityM, int
  **DestinationRanks) {

  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  auto &ConnectivityMCPP = *reinterpret_cast<ovk::connectivity_m *>(ConnectivityM);
  ConnectivityMCPP.RestoreDestinationRanks();

  *DestinationRanks = nullptr;

}
