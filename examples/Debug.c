// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <overkit.h>
#include "Global.h"
// #include "Grid.h"
// #include "Connectivity.h"
// #include "Exchange.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define Printf(...) printf(__VA_ARGS__); fflush(stdout)

typedef struct {
  int id;
  char name[64];
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  int comm_dims[2];
  int comm_coords[2];
  int global_size[2];
  int local_size[2];
  int local_count;
  int is[2], ie[2];
  double *xyz;
} input_grid;

typedef struct {
  double *values;
} input_state;

void CreateInputs(int N, int *NumLocalGrids, input_grid **Grids, input_state **States);
void DestroyInputs(int NumLocalGrids, input_grid **Grids, input_state **States);

input_grid *FindLocalGrid(int NumLocalGrids, input_grid *Grids, int GridID);

int Debug(int argc, char **argv);

int main(int argc, char **argv) {

  return Debug(argc, argv);

}

int Debug(int argc, char **argv) {

  int OtherRank;
//   int iGrid, iLocalGrid;

  MPI_Init(&argc, &argv);

  int NumProcs, Rank;
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

  if (Rank == 0) {
    Printf("Running with %i processes.\n", NumProcs);
  }

  for (OtherRank = 0; OtherRank < NumProcs; ++OtherRank) {
    if (OtherRank == Rank) {
      Printf("Rank %i reporting for duty.\n", Rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

//   int N = 10000;

//   int NumLocalGrids;
//   input_grid *InputGrids;
//   input_state *InputStates;
//   CreateInputs(N, &NumLocalGrids, &InputGrids, &InputStates);

//   t_logger *Logger;
//   CreateLogger(&Logger, Rank, OVK_LOG_ALL);

//   t_error_handler *ErrorHandler;
//   CreateErrorHandler(&ErrorHandler, OVK_ERROR_HANDLER_ABORT);

//   ovk_grid *Grids[2] = {NULL, NULL};
//   for (iGrid = 0; iGrid < 2; ++iGrid) {
//     input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, iGrid+1);
//     if (InputGrid) {
//       ovk_grid_params *GridParams;
//       CreateGridParams(&GridParams, 2, MPI_COMM_WORLD);
//       ovkSetGridParamName(GridParams, InputGrid->name);
//       ovkSetGridParamComm(GridParams, InputGrid->comm);
//       ovkSetGridParamSize(GridParams, InputGrid->global_size);
//       ovkSetGridParamLocalBegin(GridParams, InputGrid->is);
//       ovkSetGridParamLocalEnd(GridParams, InputGrid->ie);
//       CreateGrid(Grids+iGrid, InputGrid->id, GridParams, Logger, ErrorHandler);
//       DestroyGridParams(&GridParams);
//     }
//   }

//   ovk_connectivity *Connectivity1;
//   CreateConnectivity(&Connectivity1, 2, MPI_COMM_WORLD, Grids[0], Grids[1], Logger, ErrorHandler);

//   ovk_connectivity *Connectivity2;
//   CreateConnectivity(&Connectivity2, 2, MPI_COMM_WORLD, Grids[1], Grids[0], Logger, ErrorHandler);

//   ovk_exchange *Exchange1;
//   CreateExchange(&Exchange1, Connectivity1, Logger, ErrorHandler);

//   ovk_exchange *Exchange2;
//   CreateExchange(&Exchange2, Connectivity2, Logger, ErrorHandler);

//   DestroyExchange(&Exchange1);
//   DestroyExchange(&Exchange2);

//   DestroyConnectivity(&Connectivity1);
//   DestroyConnectivity(&Connectivity2);

//   for (iGrid = 0; iGrid < 2; ++iGrid) {
//     input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, iGrid+1);
//     if (InputGrid) {
//       DestroyGrid(Grids+iGrid);
//     }
//   }

//   DestroyInputs(NumLocalGrids, &InputGrids, &InputStates);

  MPI_Finalize();

  return 0;

}

void CreateInputs(int N, int *NumLocalGrids, input_grid **Grids, input_state **States) {

  int d, i, j, l, m, n, p;

  int NumProcs, Rank;
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

  int NumProcsPerGrid[2];
  int GridToProcMap[4];

  if (NumProcs >= 2) {
    NumProcsPerGrid[0] = (NumProcs+1)/2;
    NumProcsPerGrid[1] = NumProcs - (NumProcs+1)/2;
    GridToProcMap[0] = 0;
    GridToProcMap[1] = (NumProcs+1)/2;
    GridToProcMap[2] = (NumProcs+1)/2;
    GridToProcMap[3] = NumProcs;
  } else {
    NumProcsPerGrid[0] = 1;
    NumProcsPerGrid[1] = 1;
    GridToProcMap[0] = 0;
    GridToProcMap[1] = 1;
    GridToProcMap[2] = 0;
    GridToProcMap[3] = 1;
  }

  *NumLocalGrids = 0;
  for (m = 0; m < 2; ++m) {
    if (Rank >= GridToProcMap[2*m] && Rank < GridToProcMap[2*m+1]) {
      *NumLocalGrids += 1;
    }
  }

  *Grids = malloc((*NumLocalGrids)*sizeof(input_grid));

  input_grid *Grids_ = *Grids;

  for (m = 0, n = 0; m < 2; ++m) {
    if (Rank >= GridToProcMap[2*m] && Rank < GridToProcMap[2*m+1]) {
      Grids_[n].id = m+1;
      ++n;
    }
  }

  for (n = 0; n < *NumLocalGrids; ++n) {
    if (Grids_[n].id == 1) {
      strcpy(Grids_[n].name, "Left");
    } else {
      strcpy(Grids_[n].name, "Right");
    }
  }

  MPI_Group WorldGroup;
  MPI_Comm_group(MPI_COMM_WORLD, &WorldGroup);

  for (m = 0; m < 2; ++m) {
    int *Ranks = malloc(NumProcsPerGrid[m]*sizeof(int));
    for (p = 0; p < NumProcsPerGrid[m]; ++p) {
      Ranks[p] = GridToProcMap[2*m] + p;
    }
    MPI_Group GridGroup;
    MPI_Comm GridComm, GridCartComm;
    MPI_Group_incl(WorldGroup, NumProcsPerGrid[m], Ranks, &GridGroup);
    MPI_Comm_create(MPI_COMM_WORLD, GridGroup, &GridComm);
    if (GridComm != MPI_COMM_NULL) {
      for (n = 0; n < *NumLocalGrids; ++n) {
        input_grid *Grid = &Grids_[n];
        if (Grid->id == m+1) {
//           int Dims[2] = {0, 0};
          int Dims[2] = {1, 0};
          MPI_Dims_create(NumProcsPerGrid[m], 2, Dims);
          int Periods[2] = {0, 0};
          MPI_Cart_create(GridComm, 2, Dims, Periods, 0, &GridCartComm);
          Grid->comm = GridCartComm;
          MPI_Comm_size(Grid->comm, &Grid->comm_size);
          MPI_Comm_rank(Grid->comm, &Grid->comm_rank);
          MPI_Cart_get(Grid->comm, 2, Grid->comm_dims, Periods, Grid->comm_coords);
        }
      }
      MPI_Comm_free(&GridComm);
    }
    MPI_Group_free(&GridGroup);
    free(Ranks);
  }

  MPI_Group_free(&WorldGroup);

  for (n = 0; n < *NumLocalGrids; ++n) {
    input_grid *Grid = &Grids_[n];
    Grid->global_size[0] = N;
    Grid->global_size[1] = N;
    for (d = 0; d < 2; ++d) {
      int NumPerRank = N/Grid->comm_dims[d];
      int Remainder = N - Grid->comm_dims[d] * NumPerRank;
      Grid->is[d] = NumPerRank*Grid->comm_coords[d] + min(Remainder, Grid->comm_coords[d]);
      Grid->ie[d] = NumPerRank*(Grid->comm_coords[d]+1) + min(Remainder, (Grid->comm_coords[d]+1));
      Grid->local_size[d] = Grid->ie[d] - Grid->is[d];
    }
    Grid->local_count = Grid->local_size[0] * Grid->local_size[1];
  }

  for (n = 0; n < *NumLocalGrids; ++n) {
    input_grid *Grid = &Grids_[n];
    Grid->xyz = malloc(Grid->local_count*2*sizeof(double));
    double xs, xe, ys, ye;
    switch (Grid->id) {
    case 1:
      xs = -1.; xe = 0.;
      ys = -1.; ye = 1.;
      break;
    case 2:
      xs = 0.; xe = 1.;
      ys = -1.; ye = 1.;
      break;
    }
    for (j = Grid->is[1]; j < Grid->ie[1]; ++j) {
      for (i = Grid->is[0]; i < Grid->ie[0]; ++i) {
        l = (i-Grid->is[0]) + (j-Grid->is[1]) * Grid->local_size[0];
        double u = (double)i/(double)(N-1);
        double v = (double)j/(double)(N-1);
        Grid->xyz[l] = (1.-u)*xs + u*xe;
        Grid->xyz[l+Grid->local_count] = (1.-v)*ys + v*ye;
      }
    }
  }

  *States = malloc((*NumLocalGrids)*sizeof(input_state));

  input_state *States_ = *States;

  for (n = 0; n < *NumLocalGrids; ++n) {
    input_grid *Grid = &Grids_[n];
    input_state *State = &States_[n];
    State->values = malloc(Grid->local_count*sizeof(double));
    double fs, fe;
    switch (Grid->id) {
    case 1:
      fs = -1.; fe = 1.;
      break;
    case 2:
      fs = 1.; fe = -1.;
      break;
    }
    for (j = Grid->is[1]; j < Grid->ie[1]; ++j) {
      for (i = Grid->is[0]; i < Grid->ie[0]; ++i) {
        l = (i-Grid->is[0]) + (j-Grid->is[1]) * Grid->local_size[0];
        double v = (double)j/(double)(N-1);
        State->values[l] = (1.-v)*fs + v*fe;
      }
    }
  }

}

void DestroyInputs(int NumLocalGrids, input_grid **Grids, input_state **States) {

  int n;

  for (n = 0; n < NumLocalGrids; ++n) {
    input_grid *Grid = (*Grids) + n;
    free(Grid->xyz);
    Grid->xyz = NULL;
    MPI_Comm_free(&Grid->comm);
  }
  free(*Grids);
  *Grids = NULL;

  for (n = 0; n < NumLocalGrids; ++n) {
    input_state *State = (*States) + n;
    free(State->values);
    State->values = NULL;
  }
  free(*States);
  *States = NULL;

}

input_grid *FindLocalGrid(int NumLocalGrids, input_grid *Grids, int GridID) {

  int n;

  for (n = 0; n < NumLocalGrids; ++n) {
    if (Grids[n].id == GridID) return Grids + n;
  }

  return NULL;

}
