// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_OVERLAP_M_H_INCLUDED
#define OVK_CORE_C_OVERLAP_M_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core-c/Grid.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_overlap_m;
typedef struct ovk_overlap_m ovk_overlap_m;

void ovkGetOverlapMContextC(const ovk_overlap_m *OverlapM, const ovk_context **Context);
void ovkGetOverlapMContext(ovk_overlap_m *OverlapM, ovk_context **Context);
void ovkGetOverlapMSharedContext(ovk_overlap_m *OverlapM, ovk_shared_context **Context);

void ovkGetOverlapMGrid(const ovk_overlap_m *OverlapM, const ovk_grid **Grid);

void ovkGetOverlapMDestinationGridInfo(const ovk_overlap_m *OverlapM, ovk_grid_info
  **DestinationGridInfo);

void ovkGetOverlapMDimension(const ovk_overlap_m *OverlapM, int *NumDims);
void ovkGetOverlapMComm(const ovk_overlap_m *OverlapM, MPI_Comm *Comm);
void ovkGetOverlapMCommSize(const ovk_overlap_m *OverlapM, int *CommSize);
void ovkGetOverlapMCommRank(const ovk_overlap_m *OverlapM, int *CommRank);

long long ovkGetOverlapMSize(const ovk_overlap_m *OverlapM);

void ovkResizeOverlapM(ovk_overlap_m *OverlapM, long long NumCells);

void ovkGetOverlapMCells(const ovk_overlap_m *OverlapM, int Dimension, const int **Cells);
bool ovkEditingOverlapMCells(const ovk_overlap_m *OverlapM);
void ovkEditOverlapMCells(ovk_overlap_m *OverlapM, int Dimension, int **Cells);
void ovkRestoreOverlapMCells(ovk_overlap_m *OverlapM, int Dimension, int **Cells);

void ovkGetOverlapMCoords(const ovk_overlap_m *OverlapM, int Dimension, const double **Coords);
bool ovkEditingOverlapMCoords(const ovk_overlap_m *OverlapM);
void ovkEditOverlapMCoords(ovk_overlap_m *OverlapM, int Dimension, double **Coords);
void ovkRestoreOverlapMCoords(ovk_overlap_m *OverlapM, int Dimension, double **Coords);

void ovkGetOverlapMDestinations(const ovk_overlap_m *OverlapM, int Dimension, const int
  **Destinations);
bool ovkEditingOverlapMDestinations(const ovk_overlap_m *OverlapM);
void ovkEditOverlapMDestinations(ovk_overlap_m *OverlapM, int Dimension, int **Destinations);
void ovkRestoreOverlapMDestinations(ovk_overlap_m *OverlapM, int Dimension, int **Destinations);

void ovkGetOverlapMDestinationRanks(const ovk_overlap_m *OverlapM, const int **DestinationRanks);
bool ovkEditingOverlapMDestinationRanks(const ovk_overlap_m *OverlapM);
void ovkEditOverlapMDestinationRanks(ovk_overlap_m *OverlapM, int **DestinationRanks);
void ovkRestoreOverlapMDestinationRanks(ovk_overlap_m *OverlapM, int **DestinationRanks);

#ifdef __cplusplus
}
#endif

#endif
