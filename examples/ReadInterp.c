// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <overkit.h>

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
// #include <sys/types.h>
// #include <unistd.h>

#define min(a, b) ovk_min(a, b)
#define max(a, b) ovk_max(a, b)

#define Printf(...) printf(__VA_ARGS__); fflush(stdout)

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
  int id;
  char name[64];
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  int comm_dims[2];
  int comm_coords[2];
  int global_size[2];
  bool periodic[2];
  double periodic_length[2];
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

double StateFuncOne(double u, double v);
double StateFuncTwo(double u, double v);

int main(int argc, char **argv) {

  int OtherRank;
  int iGrid, iLocalGrid;
  int iSend, iReceive;
  int i, j, l;

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
//       pid_t ProcessID = getpid();
//       Printf("Rank %i (pid %i) reporting for duty.\n", Rank, ProcessID);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

//   sleep(10);

  int N = 1000;
  const char *HOPath = "../../data/small/XINTOUT.HO.2D";
  const char *XPath = "../../data/small/XINTOUT.X.2D";

  int NumLocalGrids;
  input_grid *InputGrids;
  input_state *InputStates;
  CreateInputs(N, &NumLocalGrids, &InputGrids, &InputStates);

  ovk_context_params *ContextParams;
  ovkCreateContextParams(&ContextParams);
  ovkSetContextParamComm(ContextParams, MPI_COMM_WORLD);
  ovkSetContextParamLogLevel(ContextParams, OVK_LOG_ALL);
  ovkSetContextParamProfiling(ContextParams, true);

  ovk_context *Context;
  ovk_error CreateError;
  ovkCreateContext(&Context, &ContextParams, &CreateError);
  if (CreateError != OVK_ERROR_NONE) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  ovk_shared_context *SharedContext;
  ovkShareContext(&Context, &SharedContext);

  ovk_domain_params *DomainParams;
  ovkCreateDomainParams(&DomainParams);
  ovkSetDomainParamName(DomainParams, "Domain");
  ovkSetDomainParamDimension(DomainParams, 2);
  ovkSetDomainParamComm(DomainParams, MPI_COMM_WORLD);

  ovk_domain *Domain;
  ovkCreateDomain(&Domain, SharedContext, &DomainParams);

  const int CONNECTIVITY_ID = 1;
  ovkCreateComponent(Domain, CONNECTIVITY_ID, OVK_COMPONENT_TYPE_CONNECTIVITY, NULL);

  int GridIDs[] = {1, 2};
  ovk_grid_params *MaybeGridParams[] = {NULL, NULL};

  for (iGrid = 0; iGrid < 2; ++iGrid) {
    int GridID = iGrid+1;
    input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, GridID);
    if (InputGrid) {
      ovkCreateGridParams(MaybeGridParams+iGrid);
      ovk_grid_params *GridParams = MaybeGridParams[iGrid];
      ovkSetGridParamName(GridParams, InputGrid->name);
      ovkSetGridParamDimension(GridParams, 2);
      ovkSetGridParamComm(GridParams, InputGrid->comm);
      ovkSetGridParamSize(GridParams, InputGrid->global_size);
      ovkSetGridParamPeriodic(GridParams, InputGrid->periodic);
      ovkSetGridParamPeriodicStorage(GridParams, OVK_PERIODIC_STORAGE_UNIQUE);
      ovkSetGridParamPeriodicLength(GridParams, InputGrid->periodic_length);
      ovkSetGridParamLocalRange(GridParams, InputGrid->is, InputGrid->ie);
    }
  }

  ovkCreateGrids(Domain, 2, GridIDs, MaybeGridParams);

  ovk_error ImportError;
  ovkImportXINTOUT(Domain, CONNECTIVITY_ID, HOPath, XPath, 0, MPI_INFO_NULL, &ImportError);
  if (ImportError != OVK_ERROR_NONE) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  ovk_exchanger *Exchanger;
  ovkCreateExchanger(&Exchanger, SharedContext, NULL);

  ovk_exchanger_bindings *ExchangerBindings;
  ovkCreateExchangerBindings(&ExchangerBindings);
  ovkSetExchangerBindingsConnectivityComponentID(ExchangerBindings, 1);
  ovkBindExchanger(Exchanger, Domain, &ExchangerBindings);

  ovk_connectivity_component *ConnectivityComponent;
  ovkGetComponent(Domain, CONNECTIVITY_ID, OVK_COMPONENT_TYPE_CONNECTIVITY, &ConnectivityComponent);

  int NumSends = 0;
  int NumReceives = 0;
  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    int LocalGridID = InputGrids[iLocalGrid].id;
    for (iGrid = 0; iGrid < 2; ++iGrid) {
      int OtherGridID = iGrid+1;
      if (ovkConnectivityExists(ConnectivityComponent, LocalGridID, OtherGridID)) ++NumSends;
      if (ovkConnectivityExists(ConnectivityComponent, OtherGridID, LocalGridID)) ++NumReceives;
    }
  }

  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    input_grid *InputGrid = InputGrids+iLocalGrid;
    int LocalGridID = InputGrid->id;
    for (iGrid = 0; iGrid < 2; ++iGrid) {
      int OtherGridID = iGrid+1;
      if (ovkConnectivityExists(ConnectivityComponent, LocalGridID, OtherGridID)) {
        ovkCreateExchangerCollect(Exchanger, LocalGridID, OtherGridID, 1, OVK_COLLECT_INTERPOLATE,
          OVK_DOUBLE, 1, InputGrid->is, InputGrid->ie, OVK_COLUMN_MAJOR);
        ovkCreateExchangerSend(Exchanger, LocalGridID, OtherGridID, 1, OVK_DOUBLE, 1, 1);
      }
    }
  }

  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    input_grid *InputGrid = InputGrids+iLocalGrid;
    int LocalGridID = InputGrid->id;
    for (iGrid = 0; iGrid < 2; ++iGrid) {
      int OtherGridID = iGrid+1;
      if (ovkConnectivityExists(ConnectivityComponent, OtherGridID, LocalGridID)) {
        ovkCreateExchangerReceive(Exchanger, OtherGridID, LocalGridID, 1, OVK_DOUBLE, 1, 1);
        ovkCreateExchangerDisperse(Exchanger, OtherGridID, LocalGridID, 1, OVK_DISPERSE_OVERWRITE,
          OVK_DOUBLE, 1, InputGrid->is, InputGrid->ie, OVK_COLUMN_MAJOR);
      }
    }
  }

  double **SendBuffers = malloc(NumSends*sizeof(double *));
  double **ReceiveBuffers = malloc(NumReceives*sizeof(double *));
  iSend = 0;
  iReceive = 0;
  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    int LocalGridID = InputGrids[iLocalGrid].id;
    for (iGrid = 0; iGrid < 2; ++iGrid) {
      int OtherGridID = iGrid+1;
      if (ovkConnectivityExists(ConnectivityComponent, LocalGridID, OtherGridID)) {
        const ovk_connectivity_m *ConnectivityM;
        ovkGetConnectivityM(ConnectivityComponent, LocalGridID, OtherGridID, &ConnectivityM);
        long long NumDonors;
        ovkGetConnectivityMCount(ConnectivityM, &NumDonors);
        SendBuffers[iSend] = malloc(NumDonors*sizeof(double));
        ++iSend;
      }
      if (ovkConnectivityExists(ConnectivityComponent, OtherGridID, LocalGridID)) {
        const ovk_connectivity_n *ConnectivityN;
        ovkGetConnectivityN(ConnectivityComponent, OtherGridID, LocalGridID, &ConnectivityN);
        long long NumReceivers;
        ovkGetConnectivityNCount(ConnectivityN, &NumReceivers);
        ReceiveBuffers[iReceive] = malloc(NumReceivers*sizeof(double));
        ++iReceive;
      }
    }
  }

  int iExchange = 0;

  int NumExchanges = 1;

  for (iExchange = 0; iExchange < NumExchanges; ++iExchange) {

    iSend = 0;
    for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
      input_grid *InputGrid = InputGrids+iLocalGrid;
      input_state *InputState = InputStates+iLocalGrid;
      int LocalGridID = InputGrid->id;
      for (iGrid = 0; iGrid < 2; ++iGrid) {
        int OtherGridID = iGrid+1;
        if (ovkConnectivityExists(ConnectivityComponent, LocalGridID, OtherGridID)) {
          ovkExchangerCollect(Exchanger, LocalGridID, OtherGridID, 1, &InputState->values,
            &SendBuffers[iSend]);
          ++iSend;
        }
      }
    }

    ovk_request **Requests = malloc((NumSends+NumReceives)*sizeof(ovk_request *));

    iReceive = 0;
    for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
      int LocalGridID = InputGrids[iLocalGrid].id;
      for (iGrid = 0; iGrid < 2; ++iGrid) {
        int OtherGridID = iGrid+1;
        if (ovkConnectivityExists(ConnectivityComponent, OtherGridID, LocalGridID)) {
          ovkExchangerReceive(Exchanger, OtherGridID, LocalGridID, 1, &ReceiveBuffers[iReceive],
            &Requests[iReceive]);
          ++iReceive;
        }
      }
    }

    iSend = 0;
    for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
      int LocalGridID = InputGrids[iLocalGrid].id;
      for (iGrid = 0; iGrid < 2; ++iGrid) {
        int OtherGridID = iGrid+1;
        if (ovkConnectivityExists(ConnectivityComponent, LocalGridID, OtherGridID)) {
          ovkExchangerSend(Exchanger, LocalGridID, OtherGridID, 1, &SendBuffers[iSend],
            &Requests[NumReceives+iSend]);
          ++iSend;
        }
      }
    }

    ovkWaitAll(NumSends+NumReceives, Requests);

    free(Requests);

    iReceive = 0;
    for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
      input_grid *InputGrid = InputGrids+iLocalGrid;
      input_state *InputState = InputStates+iLocalGrid;
      int LocalGridID = InputGrid->id;
      for (iGrid = 0; iGrid < 2; ++iGrid) {
        int OtherGridID = iGrid+1;
        if (ovkConnectivityExists(ConnectivityComponent, OtherGridID, LocalGridID)) {
          ovkExchangerDisperse(Exchanger, OtherGridID, LocalGridID, 1, &ReceiveBuffers[iReceive],
            &InputState->values);
          ++iReceive;
        }
      }
    }

  }

  for (iSend = 0; iSend < NumSends; ++iSend) {
    free(SendBuffers[iSend]);
  }
  free(SendBuffers);

  for (iReceive = 0; iReceive < NumReceives; ++iReceive) {
    free(ReceiveBuffers[iReceive]);
  }
  free(ReceiveBuffers);

  for (OtherRank = 0; OtherRank < NumProcs; ++OtherRank) {
    if (OtherRank == Rank) {
      for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
        input_grid *InputGrid = InputGrids+iLocalGrid;
        input_state *InputState = InputStates+iLocalGrid;
        for (j = InputGrid->is[1]; j != InputGrid->ie[1]; ++j) {
          for (i = InputGrid->is[0]; i != InputGrid->ie[0]; ++i) {
            double u = (double)i/(double)N;
            double v = (double)j/(double)N;
            l = (i-InputGrid->is[0]) + (j-InputGrid->is[1]) * InputGrid->local_size[0];
            double Expected = 1e10;
            switch (InputGrid->id) {
            case 1:
              Expected = StateFuncTwo(u, v);
              break;
            case 2:
              Expected = StateFuncOne(u, v);
              break;
            }
            if (fabs(InputState->values[l] - Expected) > 1.e-12) {
              Printf("Incorrect value at (%i,%i) on grid %i: Got %f, expected %f\n",
                i, j, InputGrid->id, InputState->values[l], Expected);
            }
          }
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  ovkDestroyExchanger(&Exchanger);
  ovkDestroyDomain(&Domain);
  ovkResetSharedContext(&SharedContext);

  DestroyInputs(NumLocalGrids, &InputGrids, &InputStates);

  MPI_Finalize();

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
      strcpy(Grids_[n].name, "One");
    } else {
      strcpy(Grids_[n].name, "Two");
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
          int Dims[2] = {0, 0};
          MPI_Dims_create(NumProcsPerGrid[m], 2, Dims);
          int Periods[2] = {1, 1};
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
    Grid->periodic[0] = true;
    Grid->periodic[1] = true;
    Grid->periodic_length[0] = 1.;
    Grid->periodic_length[1] = 1.;
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
      xs = -0.5; xe = 0.5;
      ys = -0.5; ye = 0.5;
      break;
    case 2:
      xs = -0.5; xe = 0.5;
      ys = -0.5; ye = 0.5;
      break;
    }
    for (j = Grid->is[1]; j < Grid->ie[1]; ++j) {
      for (i = Grid->is[0]; i < Grid->ie[0]; ++i) {
        l = (i-Grid->is[0]) + (j-Grid->is[1]) * Grid->local_size[0];
        double u = (double)i/(double)N;
        double v = (double)j/(double)N;
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
    for (j = Grid->is[1]; j < Grid->ie[1]; ++j) {
      for (i = Grid->is[0]; i < Grid->ie[0]; ++i) {
        double u = (double)i/(double)N;
        double v = (double)j/(double)N;
        l = (i-Grid->is[0]) + (j-Grid->is[1]) * Grid->local_size[0];
        switch (Grid->id) {
        case 1:
          State->values[l] = StateFuncOne(u, v);
          break;
        case 2:
          State->values[l] = StateFuncTwo(u, v);
          break;
        }
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

double StateFuncOne(double u, double v) {

//   return cos(2.*M_PI*v);
//   return (1.-v)*(-1.) + v*(1.);
  return v < 0.5 ? -1. : 0;

}

double StateFuncTwo(double u, double v) {

//   return sin(2.*M_PI*v);
  return v < 0.5 ? 0. : 1;

}
