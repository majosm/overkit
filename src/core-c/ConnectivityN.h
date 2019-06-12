// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_CONNECTIVITY_N_H_INCLUDED
#define OVK_CORE_C_CONNECTIVITY_N_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core-c/Grid.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity_n;
typedef struct ovk_connectivity_n ovk_connectivity_n;

void ovkGetConnectivityNContextC(const ovk_connectivity_n *ConnectivityN, const ovk_context
  **Context);
void ovkGetConnectivityNContext(ovk_connectivity_n *ConnectivityN, ovk_context **Context);
void ovkGetConnectivityNSharedContext(ovk_connectivity_n *ConnectivityN, ovk_shared_context
  **Context);

void ovkGetConnectivityNGridID(const ovk_connectivity_n *ConnectivityN, int *GridID);
void ovkGetConnectivityNGrid(const ovk_connectivity_n *ConnectivityN, const ovk_grid **Grid);

void ovkGetConnectivityNSourceGridID(const ovk_connectivity_n *ConnectivityN, int *SourceGridID);
void ovkGetConnectivityNSourceGridInfo(const ovk_connectivity_n *ConnectivityN, const ovk_grid_info
  **SourceGridInfo);

void ovkGetConnectivityNDimension(const ovk_connectivity_n *ConnectivityN, int *NumDims);
void ovkGetConnectivityNComm(const ovk_connectivity_n *ConnectivityN, MPI_Comm *Comm);
void ovkGetConnectivityNCommSize(const ovk_connectivity_n *ConnectivityN, int *CommSize);
void ovkGetConnectivityNCommRank(const ovk_connectivity_n *ConnectivityN, int *CommRank);
void ovkGetConnectivityNCount(const ovk_connectivity_n *ConnectivityN, long long *Count);

void ovkResizeConnectivityN(ovk_connectivity_n *ConnectivityN, long long Count);

void ovkGetConnectivityNPoints(const ovk_connectivity_n *ConnectivityN, int Dimension, const int
  **Points);
void ovkEditConnectivityNPoints(ovk_connectivity_n *ConnectivityN, int Dimension, int **Points);
void ovkRestoreConnectivityNPoints(ovk_connectivity_n *ConnectivityN, int Dimension, int **Points);

void ovkGetConnectivityNSources(const ovk_connectivity_n *ConnectivityN, int iDim, const int
  **Sources);
void ovkEditConnectivityNSources(ovk_connectivity_n *ConnectivityN, int iDim, int **Sources);
void ovkRestoreConnectivityNSources(ovk_connectivity_n *ConnectivityN, int iDim, int **Sources);

void ovkGetConnectivityNSourceRanks(const ovk_connectivity_n *ConnectivityN, const int
  **SourceRanks);
void ovkEditConnectivityNSourceRanks(ovk_connectivity_n *ConnectivityN, int **SourceRanks);
void ovkRestoreConnectivityNSourceRanks(ovk_connectivity_n *ConnectivityN, int **SourceRanks);

#ifdef __cplusplus
}
#endif

#endif
