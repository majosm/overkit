// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_CONNECTIVITY_COMPONENT_H_INCLUDED
#define OVK_CORE_C_CONNECTIVITY_COMPONENT_H_INCLUDED

#include <ovk/core-c/ConnectivityM.h>
#include <ovk/core-c/ConnectivityN.h>
#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core/ConnectivityComponent.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity_component;
typedef struct ovk_connectivity_component ovk_connectivity_component;

struct ovk_connectivity_component_params;
typedef struct ovk_connectivity_component_params ovk_connectivity_component_params;

int ovkConnectivityCount(const ovk_connectivity_component *ConnectivityComponent);

bool ovkConnectivityExists(const ovk_connectivity_component *ConnectivityComponent, int MGridID,
  int NGridID);

void ovkCreateConnectivity(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID);
void ovkCreateConnectivities(ovk_connectivity_component *ConnectivityComponent, int Count, const int
  *MGridIDs, const int *NGridIDs);

void ovkDestroyConnectivity(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID);
void ovkDestroyConnectivities(ovk_connectivity_component *ConnectivityComponent, int Count, const
  int *MGridIDs, const int *NGridIDs);

int LocalConnectivityMCount(const ovk_connectivity_component *ConnectivityComponent);
int LocalConnectivityMCountForGrid(const ovk_connectivity_component *ConnectivityComponent, int
  MGridID);

void ovkGetConnectivityM(const ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, const ovk_connectivity_m **ConnectivityM);
bool ovkEditingConnectivityM(const ovk_connectivity_component *ConnectivityComponent, int MGridID,
  int NGridID);
void ovkEditConnectivityM(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, ovk_connectivity_m **ConnectivityM);
void ovkRestoreConnectivityM(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, ovk_connectivity_m **ConnectivityM);

int LocalConnectivityNCount(const ovk_connectivity_component *ConnectivityComponent);
int LocalConnectivityNCountForGrid(const ovk_connectivity_component *ConnectivityComponent, int
  NGridID);

void ovkGetConnectivityN(const ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, const ovk_connectivity_n **ConnectivityN);
bool ovkEditingConnectivityN(const ovk_connectivity_component *ConnectivityComponent, int MGridID,
  int NGridID);
void ovkEditConnectivityN(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, ovk_connectivity_n **ConnectivityN);
void ovkRestoreConnectivityN(ovk_connectivity_component *ConnectivityComponent, int MGridID, int
  NGridID, ovk_connectivity_n **ConnectivityN);

void ovkCreateConnectivityComponentParams(ovk_connectivity_component_params **Params);
void ovkDestroyConnectivityComponentParams(ovk_connectivity_component_params **Params);
void ovkGetConnectivityComponentParamName(const ovk_connectivity_component_params *Params, char
  *Name);
void ovkSetConnectivityComponentParamName(ovk_connectivity_component_params *Params, const char
  *Name);

#ifdef __cplusplus
}
#endif

#endif
