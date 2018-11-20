// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Connectivity.h"

#include "ovk/core-c/ConnectivityD.h"
#include "ovk/core-c/ConnectivityR.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/ConnectivityD.hpp"
#include "ovk/core/ConnectivityR.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <cstring>
#include <string>

void ovkGetConnectivityDonorGridID(const ovk_connectivity *Connectivity, int *DonorGridID) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(DonorGridID, "Invalid donor grid ID pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  ovk::GetConnectivityDonorGridID(ConnectivityCPP, *DonorGridID);

}

void ovkGetConnectivityReceiverGridID(const ovk_connectivity *Connectivity, int *ReceiverGridID) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(ReceiverGridID, "Invalid receiver grid ID pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  ovk::GetConnectivityReceiverGridID(ConnectivityCPP, *ReceiverGridID);

}

void ovkGetConnectivityName(const ovk_connectivity *Connectivity, char *Name) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);

  std::string NameCPP;
  ovk::GetConnectivityName(ConnectivityCPP, NameCPP);

  std::strcpy(Name, NameCPP.c_str());

}

void ovkGetConnectivityDimension(const ovk_connectivity *Connectivity, int *NumDims) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  ovk::GetConnectivityDimension(ConnectivityCPP, *NumDims);

}

void ovkGetConnectivityComm(const ovk_connectivity *Connectivity, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  ovk::GetConnectivityComm(ConnectivityCPP, *Comm);

}

void ovkGetConnectivityCommSize(const ovk_connectivity *Connectivity, int *CommSize) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  ovk::GetConnectivityCommSize(ConnectivityCPP, *CommSize);

}

void ovkGetConnectivityCommRank(const ovk_connectivity *Connectivity, int *CommRank) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  ovk::GetConnectivityCommRank(ConnectivityCPP, *CommRank);

}

void ovkGetConnectivityDonorGridInfo(const ovk_connectivity *Connectivity,
  const ovk_grid_info **DonorGridInfo) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(DonorGridInfo, "Invalid donor grid info pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  
  const ovk::grid_info *DonorGridInfoCPPPtr;
  ovk::GetConnectivityDonorGridInfo(ConnectivityCPP, DonorGridInfoCPPPtr);

  *DonorGridInfo = reinterpret_cast<const ovk_grid_info *>(DonorGridInfoCPPPtr);

}

void ovkGetConnectivityReceiverGridInfo(const ovk_connectivity *Connectivity,
  const ovk_grid_info **ReceiverGridInfo) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(ReceiverGridInfo, "Invalid receiver grid info pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  
  const ovk::grid_info *ReceiverGridInfoCPPPtr;
  ovk::GetConnectivityReceiverGridInfo(ConnectivityCPP, ReceiverGridInfoCPPPtr);

  *ReceiverGridInfo = reinterpret_cast<const ovk_grid_info *>(ReceiverGridInfoCPPPtr);

}

bool ovkRankHasConnectivityDonorSide(const ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  return ovk::RankHasConnectivityDonorSide(ConnectivityCPP);

}

void ovkGetConnectivityDonorSide(const ovk_connectivity *Connectivity, const ovk_connectivity_d
  **Donors) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);

  const ovk::connectivity_d *DonorsCPPPtr;
  ovk::GetConnectivityDonorSide(ConnectivityCPP, DonorsCPPPtr);

  *Donors = reinterpret_cast<const ovk_connectivity_d *>(DonorsCPPPtr);

}

void ovkEditConnectivityDonorSideLocal(ovk_connectivity *Connectivity, ovk_connectivity_d **Donors)
  {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<ovk::connectivity *>(Connectivity);

  ovk::connectivity_d *DonorsCPPPtr;
  ovk::EditConnectivityDonorSideLocal(ConnectivityCPP, DonorsCPPPtr);

  *Donors = reinterpret_cast<ovk_connectivity_d *>(DonorsCPPPtr);

}

void ovkEditConnectivityDonorSideRemote(ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<ovk::connectivity *>(Connectivity);
  ovk::EditConnectivityDonorSideRemote(ConnectivityCPP);
  
}

void ovkReleaseConnectivityDonorSideLocal(ovk_connectivity *Connectivity, ovk_connectivity_d
  **Donors) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(*Donors, "Invalid donors pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<ovk::connectivity *>(Connectivity);

  auto DonorsCPPPtr = reinterpret_cast<ovk::connectivity_d *>(*Donors);
  ovk::ReleaseConnectivityDonorSideLocal(ConnectivityCPP, DonorsCPPPtr);

  *Donors = nullptr;

}

void ovkReleaseConnectivityDonorSideRemote(ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<ovk::connectivity *>(Connectivity);
  ovk::ReleaseConnectivityDonorSideRemote(ConnectivityCPP);
  
}

bool ovkRankHasConnectivityReceiverSide(const ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);
  return ovk::RankHasConnectivityReceiverSide(ConnectivityCPP);

}

void ovkGetConnectivityReceiverSide(const ovk_connectivity *Connectivity, const ovk_connectivity_r
  **Receivers) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<const ovk::connectivity *>(Connectivity);

  const ovk::connectivity_r *ReceiversCPPPtr;
  ovk::GetConnectivityReceiverSide(ConnectivityCPP, ReceiversCPPPtr);

  *Receivers = reinterpret_cast<const ovk_connectivity_r *>(ReceiversCPPPtr);

}

void ovkEditConnectivityReceiverSideLocal(ovk_connectivity *Connectivity,
  ovk_connectivity_r **Receivers) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<ovk::connectivity *>(Connectivity);

  ovk::connectivity_r *ReceiversCPPPtr;
  ovk::EditConnectivityReceiverSideLocal(ConnectivityCPP, ReceiversCPPPtr);

  *Receivers = reinterpret_cast<ovk_connectivity_r *>(ReceiversCPPPtr);

}

void ovkEditConnectivityReceiverSideRemote(ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<ovk::connectivity *>(Connectivity);
  ovk::EditConnectivityReceiverSideRemote(ConnectivityCPP);
  
}

void ovkReleaseConnectivityReceiverSideLocal(ovk_connectivity *Connectivity, ovk_connectivity_r
  **Receivers) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(*Receivers, "Invalid receivers pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<ovk::connectivity *>(Connectivity);

  auto ReceiversCPPPtr = reinterpret_cast<ovk::connectivity_r *>(*Receivers);
  ovk::ReleaseConnectivityReceiverSideLocal(ConnectivityCPP, ReceiversCPPPtr);

  *Receivers = nullptr;

}

void ovkReleaseConnectivityReceiverSideRemote(ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  auto &ConnectivityCPP = *reinterpret_cast<ovk::connectivity *>(Connectivity);
  ovk::ReleaseConnectivityReceiverSideRemote(ConnectivityCPP);
  
}

void ovkGetConnectivityInfoDonorGridID(const ovk_connectivity_info *Info, int *DonorGridID) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(DonorGridID, "Invalid donor grid ID pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::connectivity_info *>(Info);
  ovk::GetConnectivityInfoDonorGridID(InfoCPP, *DonorGridID);

}

void ovkGetConnectivityInfoReceiverGridID(const ovk_connectivity_info *Info, int *ReceiverGridID) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(ReceiverGridID, "Invalid receiver grid ID pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::connectivity_info *>(Info);
  ovk::GetConnectivityInfoReceiverGridID(InfoCPP, *ReceiverGridID);

}

void ovkGetConnectivityInfoName(const ovk_connectivity_info *Info, char *Name) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::connectivity_info *>(Info);

  std::string NameCPP;
  ovk::GetConnectivityInfoName(InfoCPP, NameCPP);

  std::strcpy(Name, NameCPP.c_str());

}

void ovkGetConnectivityInfoDimension(const ovk_connectivity_info *Info, int *NumDims) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::connectivity_info *>(Info);
  ovk::GetConnectivityInfoDimension(InfoCPP, *NumDims);

}

void ovkGetConnectivityInfoRootRank(const ovk_connectivity_info *Info, int *RootRank) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(RootRank, "Invalid root rank pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::connectivity_info *>(Info);
  ovk::GetConnectivityInfoRootRank(InfoCPP, *RootRank);

}

void ovkGetConnectivityInfoIsLocal(const ovk_connectivity_info *Info, bool *IsLocal) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(IsLocal, "Invalid is local pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::connectivity_info *>(Info);
  ovk::GetConnectivityInfoIsLocal(InfoCPP, *IsLocal);

}
