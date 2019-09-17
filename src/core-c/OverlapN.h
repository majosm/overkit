// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_OVERLAP_N_H_INCLUDED
#define OVK_CORE_C_OVERLAP_N_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core-c/Grid.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_overlap_n;
typedef struct ovk_overlap_n ovk_overlap_n;

void ovkGetOverlapNContextC(const ovk_overlap_n *OverlapN, const ovk_context **Context);
void ovkGetOverlapNContext(ovk_overlap_n *OverlapN, ovk_context **Context);
void ovkGetOverlapNSharedContext(ovk_overlap_n *OverlapN, ovk_shared_context **Context);

void ovkGetOverlapNGrid(const ovk_overlap_n *OverlapN, const ovk_grid **Grid);

void ovkGetOverlapNSourceGridInfo(const ovk_overlap_n *OverlapN, const ovk_grid_info
  **SourceGridInfo);

void ovkGetOverlapNDimension(const ovk_overlap_n *OverlapN, int *NumDims);
void ovkGetOverlapNComm(const ovk_overlap_n *OverlapN, MPI_Comm *Comm);
void ovkGetOverlapNCommSize(const ovk_overlap_n *OverlapN, int *CommSize);
void ovkGetOverlapNCommRank(const ovk_overlap_n *OverlapN, int *CommRank);

long long ovkGetOverlapNSize(const ovk_overlap_n *OverlapN);

void ovkResizeOverlapN(ovk_overlap_n *OverlapN, long long NumPoints);

void ovkGetOverlapNMask(const ovk_overlap_n *OverlapN, const bool **OverlapMask);

void ovkGetOverlapNPoints(const ovk_overlap_n *OverlapN, int Dimension, const int **Points);
bool ovkEditingOverlapNPoints(const ovk_overlap_n *OverlapN);
void ovkEditOverlapNPoints(ovk_overlap_n *OverlapN, int Dimension, int **Points);
void ovkRestoreOverlapNPoints(ovk_overlap_n *OverlapN, int Dimension, int **Points);

void ovkGetOverlapNSources(const ovk_overlap_n *OverlapN, int Dimension, const int **Sources);
bool ovkEditingOverlapNSources(const ovk_overlap_n *OverlapN);
void ovkEditOverlapNSources(ovk_overlap_n *OverlapN, int Dimension, int **Sources);
void ovkRestoreOverlapNSources(ovk_overlap_n *OverlapN, int Dimension, int **Sources);

void ovkGetOverlapNSourceRanks(const ovk_overlap_n *OverlapN, const int **SourceRanks);
bool ovkEditingOverlapNSourceRanks(const ovk_overlap_n *OverlapN);
void ovkEditOverlapNSourceRanks(ovk_overlap_n *OverlapN, int **SourceRanks);
void ovkRestoreOverlapNSourceRanks(ovk_overlap_n *OverlapN, int **SourceRanks);

#ifdef __cplusplus
}
#endif

#endif
