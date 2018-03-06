// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <overkit.h>
#include <overkit-extras.h>

#include "../ovk/core/ProfileUtils.h"

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

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

  int N = 100;
  const char *HOPath = "../../data/small/XINTOUT.HO.2D";
  const char *XPath = "../../data/small/XINTOUT.X.2D";

  t_profiler *Profiler;
  CreateProfiler(&Profiler, MPI_COMM_WORLD);
  if (OVK_PROFILE) EnableProfiler(Profiler);
  int OverallTime = AddProfilerTimer(Profiler, "ReadInterp");
  int CreateTime = AddProfilerTimer(Profiler, "ReadInterp::Create");
  int DestroyTime = AddProfilerTimer(Profiler, "ReadInterp::Destroy");
  int ImportTime = AddProfilerTimer(Profiler, "ReadInterp::Import");
  int AssembleTime = AddProfilerTimer(Profiler, "ReadInterp::Assemble");
  int ExchangeTime = AddProfilerTimer(Profiler, "ReadInterp::Exchange");
  int CollectTime = AddProfilerTimer(Profiler, "ReadInterp::Exchange::Collect");
  int SendRecvTime = AddProfilerTimer(Profiler, "ReadInterp::Exchange::SendRecv");
  int DisperseTime = AddProfilerTimer(Profiler, "ReadInterp::Exchange::Disperse");

  StartProfileSync(Profiler, OverallTime, MPI_COMM_WORLD);
  StartProfileSync(Profiler, CreateTime, MPI_COMM_WORLD);

  int NumLocalGrids;
  input_grid *InputGrids;
  input_state *InputStates;
  CreateInputs(N, &NumLocalGrids, &InputGrids, &InputStates);

  ovk_context_params *ContextParams;
  ovkCreateContextParams(&ContextParams);
  ovkSetContextParamComm(ContextParams, MPI_COMM_WORLD);
  ovkSetContextParamLogLevel(ContextParams, OVK_LOG_ALL);
  ovkSetContextParamErrorHandlerType(ContextParams, OVK_ERROR_HANDLER_ABORT);

  ovk_context *Context;
  ovkCreateContext(&Context, ContextParams);

  ovkDestroyContextParams(&ContextParams);

  ovk_domain_params *DomainParams;
  ovkCreateDomainParams(Context, &DomainParams);
  ovkSetDomainParamDimension(DomainParams, 2);
  ovkSetDomainParamComm(DomainParams, MPI_COMM_WORLD);

  ovk_domain *Domain;
  ovkCreateDomain(Context, &Domain, DomainParams);

  ovkDestroyDomainParams(Context, &DomainParams);

  ovkConfigureDomain(Domain, OVK_DOMAIN_CONFIG_CONNECTIVITY | OVK_DOMAIN_CONFIG_EXCHANGE);

  for (iGrid = 0; iGrid < 2; ++iGrid) {
    int GridID = iGrid+1;
    input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, GridID);
    if (InputGrid) {
      ovk_grid_params *GridParams;
      ovkCreateGridParams(Domain, &GridParams);
      ovkSetGridParamName(GridParams, InputGrid->name);
      ovkSetGridParamComm(GridParams, InputGrid->comm);
      ovkSetGridParamSize(GridParams, InputGrid->global_size);
      ovkSetGridParamPeriodic(GridParams, InputGrid->periodic);
      ovkSetGridParamPeriodicStorage(GridParams, OVK_NO_OVERLAP_PERIODIC);
      ovkSetGridParamPeriodicLength(GridParams, InputGrid->periodic_length);
      ovkSetGridParamLocalBegin(GridParams, InputGrid->is);
      ovkSetGridParamLocalEnd(GridParams, InputGrid->ie);
      ovkCreateGridLocal(Domain, GridID, GridParams);
      ovkDestroyGridParams(Domain, &GridParams);
    } else {
      ovkCreateGridRemote(Domain, GridID);
    }
  }

  EndProfileSync(Profiler, CreateTime, MPI_COMM_WORLD);
  StartProfileSync(Profiler, ImportTime, MPI_COMM_WORLD);

  int ReadGranularityAdjust = 0;
  int Blah = NumProcs;
  while (Blah > 1024) {
    Blah >>= 1;
    ReadGranularityAdjust -= 1;
  }

  ovkEXTImportXINTOUT(Domain, HOPath, XPath, ReadGranularityAdjust, MPI_INFO_NULL);

  EndProfileSync(Profiler, ImportTime, MPI_COMM_WORLD);
  StartProfileSync(Profiler, AssembleTime, MPI_COMM_WORLD);

  ovkAssemble(Domain);

  EndProfileSync(Profiler, AssembleTime, MPI_COMM_WORLD);
  StartProfileSync(Profiler, ExchangeTime, MPI_COMM_WORLD);

  int NumSends = 0;
  int NumReceives = 0;
  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    int LocalGridID = InputGrids[iLocalGrid].id;
    for (iGrid = 0; iGrid < 2; ++iGrid) {
      int OtherGridID = iGrid+1;
      if (ovkConnectivityExists(Domain, LocalGridID, OtherGridID)) ++NumSends;
      if (ovkConnectivityExists(Domain, OtherGridID, LocalGridID)) ++NumReceives;
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
      if (ovkConnectivityExists(Domain, LocalGridID, OtherGridID)) {
        size_t NumDonors;
        ovkGetLocalDonorCount(Domain, LocalGridID, OtherGridID, &NumDonors);
        SendBuffers[iSend] = malloc(NumDonors*sizeof(double));
        ++iSend;
      }
      if (ovkConnectivityExists(Domain, OtherGridID, LocalGridID)) {
        size_t NumReceivers;
        ovkGetLocalReceiverCount(Domain, OtherGridID, LocalGridID, &NumReceivers);
        ReceiveBuffers[iReceive] = malloc(NumReceivers*sizeof(double));
        ++iReceive;
      }
    }
  }

  int iExchange = 0;

  int NumExchanges = 1;

  for (iExchange = 0; iExchange < NumExchanges; ++iExchange) {

    StartProfileSync(Profiler, CollectTime, MPI_COMM_WORLD);

    iSend = 0;
    for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
      input_grid *InputGrid = InputGrids+iLocalGrid;
      input_state *InputState = InputStates+iLocalGrid;
      int LocalGridID = InputGrid->id;
      ovk_range GridDataRange;
      ovkSetRange(&GridDataRange, 2, InputGrid->is, InputGrid->ie);
      for (iGrid = 0; iGrid < 2; ++iGrid) {
        int OtherGridID = iGrid+1;
        if (ovkConnectivityExists(Domain, LocalGridID, OtherGridID)) {
          const void *GridData = InputState->values;
          void *DonorData = SendBuffers[iSend];
          ovkCollect(Domain, LocalGridID, OtherGridID, OVK_DOUBLE, 1, OVK_COLLECT_INTERPOLATE,
            &GridDataRange, OVK_COLUMN_MAJOR, &GridData, &DonorData);
          ++iSend;
        }
      }
    }

    EndProfileSync(Profiler, CollectTime, MPI_COMM_WORLD);
    StartProfileSync(Profiler, SendRecvTime, MPI_COMM_WORLD);

    ovk_request **Requests = malloc((NumSends+NumReceives)*sizeof(ovk_request *));

    iReceive = 0;
    for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
      int LocalGridID = InputGrids[iLocalGrid].id;
      for (iGrid = 0; iGrid < 2; ++iGrid) {
        int OtherGridID = iGrid+1;
        if (ovkConnectivityExists(Domain, OtherGridID, LocalGridID)) {
          void *ReceiverData = ReceiveBuffers[iReceive];
          ovkReceive(Domain, OtherGridID, LocalGridID, OVK_DOUBLE, 1, &ReceiverData, 1,
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
        if (ovkConnectivityExists(Domain, LocalGridID, OtherGridID)) {
          const void *DonorData = SendBuffers[iSend];
          ovkSend(Domain, LocalGridID, OtherGridID, OVK_DOUBLE, 1, &DonorData, 1,
            &Requests[NumReceives+iSend]);
          ++iSend;
        }
      }
    }

    ovkWaitAll(NumSends+NumReceives, Requests);

    free(Requests);

    EndProfileSync(Profiler, SendRecvTime, MPI_COMM_WORLD);
    StartProfileSync(Profiler, DisperseTime, MPI_COMM_WORLD);

    iReceive = 0;
    for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
      input_grid *InputGrid = InputGrids+iLocalGrid;
      input_state *InputState = InputStates+iLocalGrid;
      int LocalGridID = InputGrid->id;
      ovk_range GridDataRange;
      ovkSetRange(&GridDataRange, 2, InputGrid->is, InputGrid->ie);
      for (iGrid = 0; iGrid < 2; ++iGrid) {
        int OtherGridID = iGrid+1;
        if (ovkConnectivityExists(Domain, OtherGridID, LocalGridID)) {
          const void *ReceiverData = ReceiveBuffers[iReceive];
          void *GridData = InputState->values;
          ovkDisperse(Domain, OtherGridID, LocalGridID, OVK_DOUBLE, 1, OVK_DISPERSE_OVERWRITE,
            &ReceiverData, &GridDataRange, OVK_COLUMN_MAJOR, &GridData);
          ++iReceive;
        }
      }
    }

    EndProfileSync(Profiler, DisperseTime, MPI_COMM_WORLD);

  }

  for (iSend = 0; iSend < NumSends; ++iSend) {
    free(SendBuffers[iSend]);
  }
  free(SendBuffers);

  for (iReceive = 0; iReceive < NumReceives; ++iReceive) {
    free(ReceiveBuffers[iReceive]);
  }
  free(ReceiveBuffers);

  EndProfileSync(Profiler, ExchangeTime, MPI_COMM_WORLD);

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

  StartProfileSync(Profiler, DestroyTime, MPI_COMM_WORLD);

  ovkDestroyContext(&Context);

  DestroyInputs(NumLocalGrids, &InputGrids, &InputStates);

  EndProfileSync(Profiler, DestroyTime, MPI_COMM_WORLD);
  EndProfileSync(Profiler, OverallTime, MPI_COMM_WORLD);

  WriteProfileTimes(Profiler, stdout);
  DestroyProfiler(&Profiler);

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

  return cos(2.*M_PI*v);

}

double StateFuncTwo(double u, double v) {

  return sin(2.*M_PI*v);

}
