// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_CONNECTIVITY_R_H_INCLUDED
#define OVK_CORE_C_CONNECTIVITY_R_H_INCLUDED

#include <ovk/core-c/Global.h>
#include <ovk/core-c/Grid.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity_r;
typedef struct ovk_connectivity_r ovk_connectivity_r;

void ovkGetConnectivityReceiverSideGridID(const ovk_connectivity_r *Receivers, int *GridID);
void ovkGetConnectivityReceiverSideSourceGridID(const ovk_connectivity_r *Receivers, int
  *SourceGridID);
void ovkGetConnectivityReceiverSideDimension(const ovk_connectivity_r *Receivers, int *NumDims);
void ovkGetConnectivityReceiverSideComm(const ovk_connectivity_r *Receivers, MPI_Comm *Comm);
void ovkGetConnectivityReceiverSideCommSize(const ovk_connectivity_r *Receivers, int *CommSize);
void ovkGetConnectivityReceiverSideCommRank(const ovk_connectivity_r *Receivers, int *CommRank);
void ovkGetConnectivityReceiverSideCount(const ovk_connectivity_r *Receivers, long long
  *NumReceivers);
void ovkGetConnectivityReceiverSideGrid(const ovk_connectivity_r *Receivers, const ovk_grid
  **ReceiverGrid);

void ovkResizeReceivers(ovk_connectivity_r *Receivers, long long NumReceivers);

void ovkGetReceiverPoints(const ovk_connectivity_r *Receivers, int Dimension, const int **Points);
void ovkEditReceiverPoints(ovk_connectivity_r *Receivers, int Dimension, int **Points);
void ovkReleaseReceiverPoints(ovk_connectivity_r *Receivers, int Dimension, int **Points);

void ovkGetReceiverSources(const ovk_connectivity_r *Receivers, int iDim, const int **Sources);
void ovkEditReceiverSources(ovk_connectivity_r *Receivers, int iDim, int **Sources);
void ovkReleaseReceiverSources(ovk_connectivity_r *Receivers, int iDim, int **Sources);

void ovkGetReceiverSourceRanks(const ovk_connectivity_r *Receivers, const int **SourceRanks);
void ovkEditReceiverSourceRanks(ovk_connectivity_r *Receivers, int **SourceRanks);
void ovkReleaseReceiverSourceRanks(ovk_connectivity_r *Receivers, int **SourceRanks);

#ifdef __cplusplus
}
#endif

#endif