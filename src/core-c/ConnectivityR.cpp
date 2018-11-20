// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/ConnectivityR.h"

#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core/ConnectivityR.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

void ovkGetConnectivityReceiverSideGridID(const ovk_connectivity_r *Receivers, int *GridID) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetConnectivityReceiverSideGridID(ReceiversCPP, *GridID);

}

void ovkGetConnectivityReceiverSideSourceGridID(const ovk_connectivity_r *Receivers, int
  *SourceGridID) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(SourceGridID, "Invalid source grid ID pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetConnectivityReceiverSideSourceGridID(ReceiversCPP, *SourceGridID);

}

void ovkGetConnectivityReceiverSideDimension(const ovk_connectivity_r *Receivers, int *NumDims) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetConnectivityReceiverSideDimension(ReceiversCPP, *NumDims);

}

void ovkGetConnectivityReceiverSideComm(const ovk_connectivity_r *Receivers, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetConnectivityReceiverSideComm(ReceiversCPP, *Comm);

}

void ovkGetConnectivityReceiverSideCommSize(const ovk_connectivity_r *Receivers, int *CommSize) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetConnectivityReceiverSideCommSize(ReceiversCPP, *CommSize);

}

void ovkGetConnectivityReceiverSideCommRank(const ovk_connectivity_r *Receivers, int *CommRank) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetConnectivityReceiverSideCommRank(ReceiversCPP, *CommRank);

}

void ovkGetConnectivityReceiverSideCount(const ovk_connectivity_r *Receivers, long long
  *NumReceivers) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(NumReceivers, "Invalid num receivers pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetConnectivityReceiverSideCount(ReceiversCPP, *NumReceivers);

}

void ovkGetConnectivityReceiverSideGrid(const ovk_connectivity_r *Receivers, const ovk_grid
  **Grid) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid receiver grid pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);

  const ovk::grid *GridCPPPtr;
  ovk::GetConnectivityReceiverSideGrid(ReceiversCPP, GridCPPPtr);

  *Grid = reinterpret_cast<const ovk_grid *>(GridCPPPtr);

}

void ovkResizeReceivers(ovk_connectivity_r *Receivers, long long NumReceivers) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");

  auto &ReceiversCPP = *reinterpret_cast<ovk::connectivity_r *>(Receivers);
  ovk::ResizeReceivers(ReceiversCPP, NumReceivers);

}

void ovkGetReceiverPoints(const ovk_connectivity_r *Receivers, int Dimension, const int **Points) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetReceiverPoints(ReceiversCPP, Dimension, *Points);

}

void ovkEditReceiverPoints(ovk_connectivity_r *Receivers, int Dimension, int **Points) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  auto &ReceiversCPP = *reinterpret_cast<ovk::connectivity_r *>(Receivers);
  ovk::EditReceiverPoints(ReceiversCPP, Dimension, *Points);

}

void ovkReleaseReceiverPoints(ovk_connectivity_r *Receivers, int Dimension, int **Points) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  auto &ReceiversCPP = *reinterpret_cast<ovk::connectivity_r *>(Receivers);
  ovk::ReleaseReceiverPoints(ReceiversCPP, Dimension, *Points);

}

void ovkGetReceiverSources(const ovk_connectivity_r *Receivers, int Dimension, const int **Sources)
  {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetReceiverSources(ReceiversCPP, Dimension, *Sources);

}

void ovkEditReceiverSources(ovk_connectivity_r *Receivers, int Dimension, int **Sources) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  auto &ReceiversCPP = *reinterpret_cast<ovk::connectivity_r *>(Receivers);
  ovk::EditReceiverSources(ReceiversCPP, Dimension, *Sources);

}

void ovkReleaseReceiverSources(ovk_connectivity_r *Receivers, int Dimension, int **Sources) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  auto &ReceiversCPP = *reinterpret_cast<ovk::connectivity_r *>(Receivers);
  ovk::ReleaseReceiverSources(ReceiversCPP, Dimension, *Sources);

}

void ovkGetReceiverSourceRanks(const ovk_connectivity_r *Receivers, const int **SourceRanks) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  auto &ReceiversCPP = *reinterpret_cast<const ovk::connectivity_r *>(Receivers);
  ovk::GetReceiverSourceRanks(ReceiversCPP, *SourceRanks);

}

void ovkEditReceiverSourceRanks(ovk_connectivity_r *Receivers, int **SourceRanks) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  auto &ReceiversCPP = *reinterpret_cast<ovk::connectivity_r *>(Receivers);
  ovk::EditReceiverSourceRanks(ReceiversCPP, *SourceRanks);

}

void ovkReleaseReceiverSourceRanks(ovk_connectivity_r *Receivers, int **SourceRanks) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  auto &ReceiversCPP = *reinterpret_cast<ovk::connectivity_r *>(Receivers);
  ovk::ReleaseReceiverSourceRanks(ReceiversCPP, *SourceRanks);

}
