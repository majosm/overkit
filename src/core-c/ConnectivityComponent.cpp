// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/ConnectivityComponent.h"

#include "ovk/core-c/ConnectivityM.h"
#include "ovk/core-c/ConnectivityN.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/ConnectivityComponent.hpp"
#include "ovk/core/ConnectivityM.hpp"
#include "ovk/core/ConnectivityN.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <utility>

extern "C" {

int ovkConnectivityCount(const ovk_connectivity_component *ConnectivityComponent) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  return ConnectivityComponentCPP.ConnectivityCount();

}

bool ovkConnectivityExists(const ovk_connectivity_component *ConnectivityComponent, int MGridID,
  int NGridID) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  return ConnectivityComponentCPP.ConnectivityExists({MGridID,NGridID});

}

void ovkCreateConnectivity(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<ovk::connectivity_component *>(
    ConnectivityComponent);
  ConnectivityComponentCPP.CreateConnectivity({MGridID,NGridID});

}

void ovkCreateConnectivities(ovk_connectivity_component *ConnectivityComponent, int Count, const
  int *MGridIDs, const int *NGridIDs) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<ovk::connectivity_component *>(
    ConnectivityComponent);

  ovk::array<ovk::elem<int,2>> ConnectivityIDs({Count});
  for (int iCreate = 0; iCreate < Count; ++iCreate) {
    ConnectivityIDs(iCreate) = {MGridIDs[iCreate],NGridIDs[iCreate]};
  }

  ConnectivityComponentCPP.CreateConnectivities(ConnectivityIDs);

}

void ovkDestroyConnectivity(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<ovk::connectivity_component *>(
    ConnectivityComponent);
  ConnectivityComponentCPP.DestroyConnectivity({MGridID,NGridID});

}

void ovkDestroyConnectivities(ovk_connectivity_component *ConnectivityComponent, int Count, const
  int *MGridIDs, const int *NGridIDs) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<ovk::connectivity_component *>(
    ConnectivityComponent);

  ovk::array<ovk::elem<int,2>> ConnectivityIDs({Count});
  for (int iDestroy = 0; iDestroy < Count; ++iDestroy) {
    ConnectivityIDs(iDestroy) = {MGridIDs[iDestroy],NGridIDs[iDestroy]};
  }

  ConnectivityComponentCPP.DestroyConnectivities(ConnectivityIDs);

}

int ovkLocalConnectivityMCount(const ovk_connectivity_component *ConnectivityComponent) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  return ConnectivityComponentCPP.LocalConnectivityMCount();

}

int ovkLocalConnectivityMCountForGrid(const ovk_connectivity_component *ConnectivityComponent, int
  MGridID) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  return ConnectivityComponentCPP.LocalConnectivityMCountForGrid(MGridID);

}

void ovkGetConnectivityM(const ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, const ovk_connectivity_m **ConnectivityM) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");
  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  *ConnectivityM = reinterpret_cast<const ovk_connectivity_m *>(&ConnectivityComponentCPP
    .ConnectivityM({MGridID,NGridID}));

}

bool ovkEditingConnectivityM(const ovk_connectivity_component *ConnectivityComponent, int MGridID,
  int NGridID) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  return ConnectivityComponentCPP.EditingConnectivityM({MGridID,NGridID});

}

void ovkEditConnectivityM(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, ovk_connectivity_m **ConnectivityM) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");
  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<ovk::connectivity_component *>(
    ConnectivityComponent);

  ovk::edit_handle<ovk::connectivity_m> EditHandle = ConnectivityComponentCPP.EditConnectivityM(
    {MGridID,NGridID});
  auto ConnectivityMCPPPtr = EditHandle.Release();

  *ConnectivityM = reinterpret_cast<ovk_connectivity_m *>(ConnectivityMCPPPtr);

}

void ovkRestoreConnectivityM(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, ovk_connectivity_m **ConnectivityM) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");
  OVK_DEBUG_ASSERT(ConnectivityM, "Invalid connectivity M pointer.");
  OVK_DEBUG_ASSERT(*ConnectivityM, "Invalid connectivity M pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<ovk::connectivity_component *>(
    ConnectivityComponent);
  ConnectivityComponentCPP.RestoreConnectivityM({MGridID,NGridID});

  *ConnectivityM = nullptr;

}

int ovkLocalConnectivityNCount(const ovk_connectivity_component *ConnectivityComponent) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  return ConnectivityComponentCPP.LocalConnectivityNCount();

}

int ovkLocalConnectivityNCountForGrid(const ovk_connectivity_component *ConnectivityComponent, int
  NGridID) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  return ConnectivityComponentCPP.LocalConnectivityNCountForGrid(NGridID);

}

void ovkGetConnectivityN(const ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, const ovk_connectivity_n **ConnectivityN) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");
  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  *ConnectivityN = reinterpret_cast<const ovk_connectivity_n *>(&ConnectivityComponentCPP
    .ConnectivityN({MGridID,NGridID}));

}

bool ovkEditingConnectivityN(const ovk_connectivity_component *ConnectivityComponent, int MGridID,
  int NGridID) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<const ovk::connectivity_component *>(
    ConnectivityComponent);
  return ConnectivityComponentCPP.EditingConnectivityN({MGridID,NGridID});

}

void ovkEditConnectivityN(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, ovk_connectivity_n **ConnectivityN) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");
  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<ovk::connectivity_component *>(
    ConnectivityComponent);

  ovk::edit_handle<ovk::connectivity_n> EditHandle = ConnectivityComponentCPP.EditConnectivityN(
    {MGridID,NGridID});
  auto ConnectivityNCPPPtr = EditHandle.Release();

  *ConnectivityN = reinterpret_cast<ovk_connectivity_n *>(ConnectivityNCPPPtr);

}

void ovkRestoreConnectivityN(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, ovk_connectivity_n **ConnectivityN) {

  OVK_DEBUG_ASSERT(ConnectivityComponent, "Invalid connectivity component pointer.");
  OVK_DEBUG_ASSERT(ConnectivityN, "Invalid connectivity N pointer.");
  OVK_DEBUG_ASSERT(*ConnectivityN, "Invalid connectivity N pointer.");

  auto &ConnectivityComponentCPP = *reinterpret_cast<ovk::connectivity_component *>(
    ConnectivityComponent);
  ConnectivityComponentCPP.RestoreConnectivityN({MGridID,NGridID});

  *ConnectivityN = nullptr;

}

void ovkCreateConnectivityComponentParams(ovk_connectivity_component_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::connectivity_component::params();

  *Params = reinterpret_cast<ovk_connectivity_component_params *>(ParamsCPPPtr);

}

void ovkDestroyConnectivityComponentParams(ovk_connectivity_component_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::connectivity_component::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetConnectivityComponentParamName(const ovk_connectivity_component_params *Params, char
  *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::connectivity_component::params *>(Params);
  std::strcpy(Name, ParamsCPP.Name().c_str());

}

void ovkSetConnectivityComponentParamName(ovk_connectivity_component_params *Params, const char
  *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::connectivity_component::params *>(Params);
  ParamsCPP.SetName(Name);

}

}
