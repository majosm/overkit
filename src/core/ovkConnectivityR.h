// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CONNECTIVITY_R_INCLUDED
#define OVK_CORE_PUBLIC_CONNECTIVITY_R_INCLUDED

#include <ovkGlobal.h>
#include <ovkGrid.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity_r_properties;
typedef struct ovk_connectivity_r_properties ovk_connectivity_r_properties;

struct ovk_connectivity_r;
typedef struct ovk_connectivity_r ovk_connectivity_r;

void ovkGetConnectivityReceiverSideProperties(const ovk_connectivity_r *Receivers,
  const ovk_connectivity_r_properties **Properties);

void ovkResizeReceivers(ovk_connectivity_r *Receivers, size_t NumReceivers);
void ovkEditReceiverPoints(ovk_connectivity_r *Receivers, int Dimension, int **Points);
void ovkReleaseReceiverPoints(ovk_connectivity_r *Receivers, int Dimension, int **Points);
void ovkEditReceiverSources(ovk_connectivity_r *Receivers, int iDim, int **Sources);
void ovkReleaseReceiverSources(ovk_connectivity_r *Receivers, int iDim, int **Sources);
void ovkEditReceiverSourceRanks(ovk_connectivity_r *Receivers, int **SourceRanks);
void ovkReleaseReceiverSourceRanks(ovk_connectivity_r *Receivers, int **SourceRanks);

void ovkGetConnectivityReceiverSideGrid(const ovk_connectivity_r *Receivers,
  const ovk_grid **ReceiverGrid);

void ovkGetConnectivityReceiverSidePropertyGridID(const ovk_connectivity_r_properties *Properties,
  int *GridID);
void ovkGetConnectivityReceiverSidePropertyDimension(const ovk_connectivity_r_properties *Properties,
  int *NumDims);
void ovkGetConnectivityReceiverSidePropertyComm(const ovk_connectivity_r_properties *Properties,
  MPI_Comm *Comm);
void ovkGetConnectivityReceiverSidePropertyCommSize(const ovk_connectivity_r_properties *Properties,
  int *CommSize);
void ovkGetConnectivityReceiverSidePropertyCommRank(const ovk_connectivity_r_properties *Properties,
  int *CommRank);
void ovkGetConnectivityReceiverSidePropertyReceiverCount(const ovk_connectivity_r_properties *Properties,
  size_t *NumReceivers);

#ifdef __cplusplus
}
#endif

#endif
