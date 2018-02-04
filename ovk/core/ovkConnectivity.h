// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CONNECTIVITY_INCLUDED
#define OVK_CORE_PUBLIC_CONNECTIVITY_INCLUDED

#include <ovk/core/ovkConnectivityD.h>
#include <ovk/core/ovkConnectivityR.h>
#include <ovk/core/ovkGlobal.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity;
typedef struct ovk_connectivity ovk_connectivity;

struct ovk_connectivity_properties;
typedef struct ovk_connectivity_properties ovk_connectivity_properties;

struct ovk_connectivity_info;
typedef struct ovk_connectivity_info ovk_connectivity_info;

void ovkGetConnectivityProperties(const ovk_connectivity *Connectivity,
  const ovk_connectivity_properties **Properties);

void ovkGetConnectivityDonorGridInfo(const ovk_connectivity *Connectivity,
  const ovk_grid_info **DonorGridInfo);
void ovkGetConnectivityReceiverGridInfo(const ovk_connectivity *Connectivity,
  const ovk_grid_info **ReceiverGridInfo);

bool ovkRankHasConnectivityDonorSide(const ovk_connectivity *Connectivity);
void ovkGetConnectivityDonorSide(const ovk_connectivity *Connectivity,
  const ovk_connectivity_d **Donors);
void ovkEditConnectivityDonorSideLocal(ovk_connectivity *Connectivity,
  ovk_connectivity_d **Donors);
void ovkEditConnectivityDonorSideRemote(ovk_connectivity *Connectivity);
void ovkReleaseConnectivityDonorSideLocal(ovk_connectivity *Connectivity,
  ovk_connectivity_d **Donors);
void ovkReleaseConnectivityDonorSideRemote(ovk_connectivity *Connectivity);

bool ovkRankHasConnectivityReceiverSide(const ovk_connectivity *Connectivity);
void ovkGetConnectivityReceiverSide(const ovk_connectivity *Connectivity,
  const ovk_connectivity_r **Receivers);
void ovkEditConnectivityReceiverSideLocal(ovk_connectivity *Connectivity,
  ovk_connectivity_r **Receivers);
void ovkEditConnectivityReceiverSideRemote(ovk_connectivity *Connectivity);
void ovkReleaseConnectivityReceiverSideLocal(ovk_connectivity *Connectivity,
  ovk_connectivity_r **Receivers);
void ovkReleaseConnectivityReceiverSideRemote(ovk_connectivity *Connectivity);

void ovkGetConnectivityPropertyDonorGridID(const ovk_connectivity_properties *Properties,
  int *DonorGridID);
void ovkGetConnectivityPropertyReceiverGridID(const ovk_connectivity_properties *Properties,
  int *ReceiverGridID);
void ovkGetConnectivityPropertyName(const ovk_connectivity_properties *Properties, char *Name);
void ovkGetConnectivityPropertyDimension(const ovk_connectivity_properties *Properties,
  int *NumDims);
void ovkGetConnectivityPropertyComm(const ovk_connectivity_properties *Properties, MPI_Comm *Comm);
void ovkGetConnectivityPropertyCommSize(const ovk_connectivity_properties *Properties,
  int *CommSize);
void ovkGetConnectivityPropertyCommRank(const ovk_connectivity_properties *Properties,
  int *CommRank);

void ovkGetConnectivityInfoDonorGridID(const ovk_connectivity_info *Info, int *DonorGridID);
void ovkGetConnectivityInfoReceiverGridID(const ovk_connectivity_info *Info, int *ReceiverGridID);
void ovkGetConnectivityInfoName(const ovk_connectivity_info *Info, char *Name);
void ovkGetConnectivityInfoDimension(const ovk_connectivity_info *Info, int *NumDims);
void ovkGetConnectivityInfoRootRank(const ovk_connectivity_info *Info, int *RootRank);

#ifdef __cplusplus
}
#endif

#endif
