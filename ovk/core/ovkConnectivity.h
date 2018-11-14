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

struct ovk_connectivity_info;
typedef struct ovk_connectivity_info ovk_connectivity_info;

void ovkGetConnectivityDonorGridID(const ovk_connectivity *Connectivity, int *DonorGridID);
void ovkGetConnectivityReceiverGridID(const ovk_connectivity *Connectivity, int *ReceiverGridID);
void ovkGetConnectivityName(const ovk_connectivity *Connectivity, char *Name);
void ovkGetConnectivityDimension(const ovk_connectivity *Connectivity, int *NumDims);
void ovkGetConnectivityComm(const ovk_connectivity *Connectivity, MPI_Comm *Comm);
void ovkGetConnectivityCommSize(const ovk_connectivity *Connectivity, int *CommSize);
void ovkGetConnectivityCommRank(const ovk_connectivity *Connectivity, int *CommRank);

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

void ovkGetConnectivityInfoDonorGridID(const ovk_connectivity_info *Info, int *DonorGridID);
void ovkGetConnectivityInfoReceiverGridID(const ovk_connectivity_info *Info, int *ReceiverGridID);
void ovkGetConnectivityInfoName(const ovk_connectivity_info *Info, char *Name);
void ovkGetConnectivityInfoDimension(const ovk_connectivity_info *Info, int *NumDims);
void ovkGetConnectivityInfoRootRank(const ovk_connectivity_info *Info, int *RootRank);

#ifdef __cplusplus
}
#endif

#endif
