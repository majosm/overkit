// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <overkit.h>

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define min(a, b) ovk_min(a, b)
#define max(a, b) ovk_max(a, b)

#define Printf(...) printf(__VA_ARGS__); fflush(stdout)

enum {
  MAX_DIMS = OVK_MAX_DIMS
};

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

void ExchangeTest(int argc, char **argv);
void AssembleTest(int argc, char **argv);

int main(int argc, char **argv) {

  ExchangeTest(argc, argv);
//   AssembleTest(argc, argv);

  return 0;

}

void ExchangeTest(int argc, char **argv) {

  int OtherRank;
  int iDim, iCoef;
  int iGrid, jGrid, iLocalGrid;
  int iConnectivity, iSend, iReceive;
  int i, j;
  size_t iDonor, iReceiver;

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

  int N = 100;

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
  ovkCreateDomainParams(&DomainParams, 2);
  ovkSetDomainParamName(DomainParams, "Domain");
  ovkSetDomainParamComm(DomainParams, MPI_COMM_WORLD);

  ovk_domain *Domain;
  ovkCreateDomain(Context, &Domain, DomainParams);

  ovkDestroyDomainParams(&DomainParams);

  ovkConfigureDomain(Domain, OVK_DOMAIN_CONFIG_CONNECTIVITY | OVK_DOMAIN_CONFIG_EXCHANGE);

  MPI_Comm SplitComm;
  MPI_Comm_split(MPI_COMM_WORLD, Rank % 2, Rank, &SplitComm);

  ovk_domain *OtherDomain = NULL;
  if (Rank % 2 == 1) {
    ovkCreateDomainParams(&DomainParams, 2);
    ovkSetDomainParamName(DomainParams, "OtherDomain");
    ovkSetDomainParamComm(DomainParams, SplitComm);
    ovkCreateDomain(Context, &OtherDomain, DomainParams);
    ovkDestroyDomain(Context, &OtherDomain);
    ovkDestroyDomainParams(&DomainParams);
  }

  MPI_Comm_free(&SplitComm);

  for (iGrid = 0; iGrid < 2; ++iGrid) {
    int GridID = iGrid+1;
    input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, GridID);
    if (InputGrid) {
      ovk_grid_params *GridParams;
      ovkCreateGridParams(&GridParams, 2);
      ovkSetGridParamName(GridParams, InputGrid->name);
      ovkSetGridParamComm(GridParams, InputGrid->comm);
      ovkSetGridParamSize(GridParams, InputGrid->global_size);
      ovk_range LocalRange;
      ovkSetRange(&LocalRange, 2, InputGrid->is, InputGrid->ie);
      ovkSetGridParamLocalRange(GridParams, &LocalRange);
      ovkCreateGridLocal(Domain, GridID, GridParams);
      ovkDestroyGridParams(&GridParams);
    } else {
      ovkCreateGridRemote(Domain, GridID);
    }
  }

//   for (iGrid = 0; iGrid < 2; ++iGrid) {
//     int GridID = iGrid+1;
//     input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, GridID);
//     ovk_grid_params *GridParams = NULL;
//     if (InputGrid) {
//       ovkCreateGridParams(&GridParams, 2);
//       ovkSetGridParamName(GridParams, InputGrid->name);
//       ovkSetGridParamComm(GridParams, InputGrid->comm);
//       ovkSetGridParamSize(GridParams, InputGrid->global_size);
//       ovk_range LocalRange;
//       ovkSetRange(&LocalRange, 2, InputGrid->is, InputGrid->ie);
//       ovkSetGridParamLocalRange(GridParams, &LocalRange);
//     }
//     ovkCreateGrid(Domain, GridID, GridParams);
//     if (InputGrid) {
//       ovkDestroyGridParams(&GridParams);
//     }
//   }

//   int GridIDs[2] = {1, 2};
//   ovk_grid_params *GridParams[2] = {NULL, NULL};
//   for (iGrid = 0; iGrid < 2; ++iGrid) {
//     input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, GridIDs[iGrid]);
//     if (InputGrid) {
//       ovkCreateGridParams(&GridParams[iGrid], 2);
//       ovkSetGridParamName(GridParams[iGrid], InputGrid->name);
//       ovkSetGridParamComm(GridParams[iGrid], InputGrid->comm);
//       ovkSetGridParamSize(GridParams[iGrid], InputGrid->global_size);
//       ovk_range LocalRange;
//       ovkSetRange(&LocalRange, 2, InputGrid->is, InputGrid->ie);
//       ovkSetGridParamLocalRange(GridParams[iGrid], &LocalRange);
//     }
//   }

//   ovkCreateGrids(Domain, 2, GridIDs, GridParams);

//   for (iGrid = 0; iGrid < 2; ++iGrid) {
//     input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, GridIDs[iGrid]);
//     if (InputGrid) {
//       ovkDestroyGridParams(&GridParams[iGrid]);
//     }
//   }

//   ovk_overlap *Overlap;
//   ovkEditOverlap(Domain, DonorGridID, ReceiverGridID, &Overlap);
//   ovk_overlap_d *OverlapD;
//   ovkEditOverlapDonorSide(Overlap, &OverlapD);
//   int *OverlappingCells[MAX_DIMS];
//   double *OverlappingCoords[2];
//   for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
//     ovkEditOverlappingCells(OverlapD, iDim, &OverlappingCells[d]);
//     ovkReleaseOverlappingCells(OverlapD, iDim, &OverlappingCells[d]);
//   }
//   for (iDim = 0; iDim < 2; ++iDim) {
//     ovkEditOverlappingCoords(OverlapD, iDim, &OverlappingCoords[d]);
//     ovkReleaseOverlappingCoords(OverlapD, iDim, &OverlappingCoords[d]);
//   }
//   ovkReleaseOverlapDonorSide(Overlap, &OverlapD);
//   
//   ovk_overlap_r *OverlapR;
//   ovkEditOverlapReceiverSide(Overlap, &OverlapR);
//   int *OverlappedPoints[MAX_DIMS];
//   for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
//     ovkEditOverlappedPoints(OverlapR, iDim, &OverlappedPoints);
//     ovkReleaseOverlappedPoints(OverlapR, iDim, &OverlappedPoints);
//   }
//   ovkReleaseOverlapReceiverSide(Overlap, &OverlapR);
//   ovkReleaseOverlap(Domain, DonorGridID, ReceiverGridID, &Overlap);



//   ovkGetLocalConnectivityCount(Domain, &NumLocalConnectivities);
//   ovk_connectivity **Connectivities = malloc(NumLocalConnectivities*sizeof(ovk_connectivity *));
//   ovkEditAllConnectivities(Domain, Connectivities);

//   int DonorGridIDs[] = {1, 2};
//   int ReceiverGridIDs[] = {2, 1};
//   ovkCreateConnectivities(Domain, 2, DonorGridIDs, ReceiverGridIDs);
//   ovkCreateExchanges(Domain, 2, DonorGridIDs, ReceiverGridIDs);

  ovk_connectivity *Connectivities[2] = {NULL, NULL};
  iConnectivity = 0;
  for (jGrid = 0; jGrid < 2; ++jGrid) {
    for (iGrid = 0; iGrid < 2; ++iGrid) {
      int DonorGridID = iGrid+1;
      int ReceiverGridID = jGrid+1;
      if (ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID)) {
        if (ovkRankHasConnectivity(Domain, DonorGridID, ReceiverGridID)) {
          ovkEditConnectivityLocal(Domain, DonorGridID, ReceiverGridID,
            &Connectivities[iConnectivity]);
          ++iConnectivity;
        } else {
          ovkEditConnectivityRemote(Domain, DonorGridID, ReceiverGridID);
        }
      }
    }
  }

  iConnectivity = 0;
  for (jGrid = 0; jGrid < 2; ++jGrid) {
    for (iGrid = 0; iGrid < 2; ++iGrid) {
      int DonorGridID = iGrid+1;
      int ReceiverGridID = jGrid+1;
      if (ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID)) {
        ovk_connectivity *Connectivity = Connectivities[iConnectivity];
        ovk_connectivity_d *Donors = NULL;
        if (ovkRankHasGrid(Domain, DonorGridID)) {
          ovkEditConnectivityDonorSideLocal(Connectivity, &Donors);
        } else {
          ovkEditConnectivityDonorSideRemote(Connectivity);
        }
        ovk_connectivity_r *Receivers = NULL;
        if (ovkRankHasGrid(Domain, ReceiverGridID)) {
          ovkEditConnectivityReceiverSideLocal(Connectivity, &Receivers);
        } else {
          ovkEditConnectivityReceiverSideRemote(Connectivity);
        }
        if (Donors) {
          input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, DonorGridID);
          bool RightEdgeOfLeftGrid = InputGrid->id == 1 && InputGrid->ie[0] == N;
          bool LeftEdgeOfRightGrid = InputGrid->id == 2 && InputGrid->is[0] == 0;
          size_t NumDonors = 0;
          if (RightEdgeOfLeftGrid || LeftEdgeOfRightGrid) {
            NumDonors = InputGrid->ie[1] - InputGrid->is[1];
          }
          ovkResizeDonors(Donors, NumDonors, 1);
          int *Extents[2][2];
          double *Coords[2];
          double *InterpCoefs[2][1];
          int *Destinations[2];
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkEditDonorExtents(Donors, iDim, &Extents[0][iDim], &Extents[1][iDim]);
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkEditDonorCoords(Donors, iDim, &Coords[iDim]);
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            for (iCoef = 0; iCoef < 1; ++iCoef) {
              ovkEditDonorInterpCoefs(Donors, iDim, iCoef, &InterpCoefs[iDim][iCoef]);
            }
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkEditDonorDestinations(Donors, iDim, &Destinations[iDim]);
          }
          if (RightEdgeOfLeftGrid || LeftEdgeOfRightGrid) {
            iDonor = 0;
            for (j = InputGrid->is[1]; j < InputGrid->ie[1]; ++j) {
              Extents[0][0][iDonor] = RightEdgeOfLeftGrid ? N-1 : 0;
              Extents[0][1][iDonor] = j;
              Extents[1][0][iDonor] = Extents[0][0][iDonor]+1;
              Extents[1][1][iDonor] = Extents[0][1][iDonor]+1;
              Coords[0][iDonor] = 0.;
              Coords[1][iDonor] = 0.;
              InterpCoefs[0][0][iDonor] = 1.;
              InterpCoefs[1][0][iDonor] = 1.;
              Destinations[0][iDonor] = RightEdgeOfLeftGrid ? 0 : N-1;
              Destinations[1][iDonor] = j;
              ++iDonor;
            }
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkReleaseDonorExtents(Donors, iDim, &Extents[0][iDim], &Extents[1][iDim]);
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkReleaseDonorCoords(Donors, iDim, &Coords[iDim]);
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            for (iCoef = 0; iCoef < 1; ++iCoef) {
              ovkReleaseDonorInterpCoefs(Donors, iDim, iCoef, &InterpCoefs[iDim][iCoef]);
            }
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkReleaseDonorDestinations(Donors, iDim, &Destinations[iDim]);
          }
        }
        if (Receivers) {
          input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, ReceiverGridID);
          bool RightEdgeOfLeftGrid = InputGrid->id == 1 && InputGrid->ie[0] == N;
          bool LeftEdgeOfRightGrid = InputGrid->id == 2 && InputGrid->is[0] == 0;
          size_t NumReceivers = 0;
          if (RightEdgeOfLeftGrid || LeftEdgeOfRightGrid) {
            NumReceivers = InputGrid->ie[1] - InputGrid->is[1];
          }
          ovkResizeReceivers(Receivers, NumReceivers);
          int *ReceiverPoints[2];
          int *Sources[2];
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkEditReceiverPoints(Receivers, iDim, &ReceiverPoints[iDim]);
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkEditReceiverSources(Receivers, iDim, &Sources[iDim]);
          }
          if (RightEdgeOfLeftGrid || LeftEdgeOfRightGrid) {
            iReceiver = 0;
            for (j = InputGrid->is[1]; j < InputGrid->ie[1]; ++j) {
              ReceiverPoints[0][iReceiver] = RightEdgeOfLeftGrid ? N-1 : 0;
              ReceiverPoints[1][iReceiver] = j;
              Sources[0][iReceiver] = RightEdgeOfLeftGrid ? 0 : N-1;
              Sources[1][iReceiver] = j;
              ++iReceiver;
            }
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkReleaseReceiverPoints(Receivers, iDim, &ReceiverPoints[iDim]);
          }
          for (iDim = 0; iDim < 2; ++iDim) {
            ovkReleaseReceiverSources(Receivers, iDim, &Sources[iDim]);
          }
        }
        if (ovkRankHasGrid(Domain, DonorGridID)) {
          ovkReleaseConnectivityDonorSideLocal(Connectivity, &Donors);
        } else {
          ovkReleaseConnectivityDonorSideRemote(Connectivity);
        }
        if (ovkRankHasGrid(Domain, ReceiverGridID)) {
          ovkReleaseConnectivityReceiverSideLocal(Connectivity, &Receivers);
        } else {
          ovkReleaseConnectivityReceiverSideRemote(Connectivity);
        }
        ++iConnectivity;
      }
    }
  }

  iConnectivity = 0;
  for (jGrid = 0; jGrid < 2; ++jGrid) {
    for (iGrid = 0; iGrid < 2; ++iGrid) {
      int DonorGridID = iGrid+1;
      int ReceiverGridID = jGrid+1;
      if (ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID)) {
        if (ovkRankHasConnectivity(Domain, DonorGridID, ReceiverGridID)) {
          ovkReleaseConnectivityLocal(Domain, DonorGridID, ReceiverGridID,
            &Connectivities[iConnectivity]);
          ++iConnectivity;
        } else {
          ovkReleaseConnectivityRemote(Domain, DonorGridID, ReceiverGridID);
        }
      }
    }
  }

  ovk_assembly_options *Options;
  int GridIDs[2] = {1, 2};
  ovkCreateAssemblyOptions(&Options, 2, 2, GridIDs);

  ovkAssemble(Domain, Options);

  ovkDestroyAssemblyOptions(&Options);

  // Exchange variant #1 -- basic exchange

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

  for (OtherRank = 0; OtherRank < NumProcs; ++OtherRank) {
    if (OtherRank == Rank) {
      for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
        input_grid *InputGrid = InputGrids+iLocalGrid;
        input_state *InputState = InputStates+iLocalGrid;
        bool RightEdgeOfLeftGrid = InputGrid->id == 1 && InputGrid->ie[0] == N;
        bool LeftEdgeOfRightGrid = InputGrid->id == 2 && InputGrid->is[0] == 0;
        if (RightEdgeOfLeftGrid) {
          i = InputGrid->ie[0]-1;
          for (j = InputGrid->is[1]; j != InputGrid->ie[1]; ++j) {
            int iSecondLast = (InputGrid->ie[0]-InputGrid->is[0])*(j-InputGrid->is[1]) +
              (i-1-InputGrid->is[0]);
            int iLast = (InputGrid->ie[0]-InputGrid->is[0])*(j-InputGrid->is[1]) +
              (i-InputGrid->is[0]);
            Printf("%lf, %lf\n", InputState->values[iSecondLast], InputState->values[iLast]);
          }
          Printf("\n");
        } else if (LeftEdgeOfRightGrid) {
          i = 0;
          for (j = InputGrid->is[1]; j != InputGrid->ie[1]; ++j) {
            int iFirst = (InputGrid->ie[0]-InputGrid->is[0])*(j-InputGrid->is[1]) +
              (i-InputGrid->is[0]);
            int iSecond = (InputGrid->ie[0]-InputGrid->is[0])*(j-InputGrid->is[1]) +
              (i+1-InputGrid->is[0]);
            Printf("%lf, %lf\n", InputState->values[iFirst], InputState->values[iSecond]);
          }
          Printf("\n");
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for (iSend = 0; iSend < NumSends; ++iSend) {
    free(SendBuffers[iSend]);
  }
  free(SendBuffers);

  for (iReceive = 0; iReceive < NumReceives; ++iReceive) {
    free(ReceiveBuffers[iReceive]);
  }
  free(ReceiveBuffers);

#ifdef NOPE

  // Exchange variant #2 -- interleaved receive/disperse

  double **SendBuffers = malloc(2*NumLocalGrids*sizeof(double *));
  double **ReceiveBuffers = malloc(2*NumLocalGrids*sizeof(double *));
  l = 0;
  for (n = 0; n < NumLocalGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      size_t NumDonors = ovkGetDonorCount(Domain, Grids[n].id, m+1);
      size_t NumReceivers = ovkGetReceiverCount(Domain, Grids[n].id, m+1);
      SendBuffers[l] = malloc(NumDonors*sizeof(double));
      ReceiveBuffers[l] = malloc(NumReceivers*sizeof(double));
      ++l;
    }
  }

  l = 0;
  for (n = 0; n < NumLocalGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      ovkCollect(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, OVK_COLLECT_INTERPOLATE,
        InputStates[n]->values, SendBuffers[l]);
      ++l;
    }
  }

  ovk_request *SendRequests = malloc(2*NumLocalGrids*sizeof(ovk_request *));
  ovk_request *ReceiveRequests = malloc(2*NumLocalGrids*sizeof(ovk_request *));
  l = 0;
  for (n = 0; n < NumLocalGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      ovkReceive(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, &ReceiveBuffers[l], 1, ReceiveRequests[l]);
      ovkSend(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, &SendBuffers[l], 1, SendRequests[l]);
      ++l;
    }
  }

  while (true) {
    ovkWaitAny(Domain, 2*NumLocalGrids, ReceiveRequests, &l);
    if (l < 0) break;
    m = l % 2;
    n = l/2 % NumLocalGrids;
    ovkDisperse(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, OVK_DISPERSE_OVERWRITE, ReceiveBuffers[l],
      InputState[n]->values);
  }
  ovkWaitAll(2*NumLocalGrids, SendRequests);

  l = 0;
  for (n = 0; n < NumLocalGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      free(SendBuffers[l]);
      free(ReceiveBuffers[l]);
      ++l;
    }
  }
  free(SendBuffers);
  free(ReceiverBuffers);

  // Exchange variant #3 -- most efficient way?

  double **ReceiveBuffers = malloc(2*NumLocalGrids*sizeof(double *));
  ovk_request *ReceiveRequests = malloc(2*NumLocalGrids*sizeof(ovk_request *));
  l = 0;
  for (n = 0; n < NumLocalGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      size_t NumReceivers = ovkGetReceiverCount(Domain, Grids[n].id, m+1);
      ReceiveBuffers[l] = malloc(NumReceivers*sizeof(double));
      ovkReceive(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, &ReceiveBuffers[l], 1, ReceiveRequests[l]);
      ++l;
    }
  }

  double **SendBuffers = malloc(2*NumLocalGrids*sizeof(double *));
  ovk_request *SendRequests = malloc(2*NumLocalGrids*sizeof(ovk_request *));
  l = 0;
  for (n = 0; n < NumLocalGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      size_t NumDonors = ovkGetDonorCount(Domain, Grids[n].id, m+1);
      SendBuffers[l] = malloc(NumDonors*sizeof(double));
      ovkCollect(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, OVK_COLLECT_INTERPOLATE,
        InputStates[n]->values, SendBuffers[l]);
      ovkSend(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, &SendBuffers[l], 1, SendRequests[l]);
      ++l;
    }
  }

  while (true) {
    ovkWaitAny(Domain, 2*NumLocalGrids, ReceiveRequests, &l);
    if (l < 0) break;
    m = l % 2;
    n = l/2 % NumLocalGrids;
    ovkDisperse(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, OVK_DISPERSE_OVERWRITE, ReceiveBuffers[l],
      InputState[n]->values);
  }
  ovkWaitAll(2*NumLocalGrids, SendRequests);

  l = 0;
  for (n = 0; n < NumLocalGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      free(SendBuffers[l]);
      free(ReceiveBuffers[l]);
      ++l;
    }
  }
  free(SendBuffers);
  free(ReceiverBuffers);

#endif

  ovkDestroyContext(&Context);

  DestroyInputs(NumLocalGrids, &InputGrids, &InputStates);

  MPI_Finalize();

}

void AssembleTest(int argc, char **argv) {

  int OtherRank;
//   int iDim;
  int iGrid, iLocalGrid;
//   size_t iPoint;
//   int i, j;

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

  int N = 100;

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
  ovkCreateDomainParams(&DomainParams, 2);
  ovkSetDomainParamName(DomainParams, "Domain");
  ovkSetDomainParamComm(DomainParams, MPI_COMM_WORLD);

  ovk_domain *Domain;
  ovkCreateDomain(Context, &Domain, DomainParams);

  ovkDestroyDomainParams(&DomainParams);

  ovkConfigureDomain(Domain, OVK_DOMAIN_CONFIG_GEOMETRY | OVK_DOMAIN_CONFIG_OVERLAP |
    OVK_DOMAIN_CONFIG_CONNECTIVITY);

  for (iGrid = 0; iGrid < 2; ++iGrid) {
    int GridID = iGrid+1;
    input_grid *InputGrid = FindLocalGrid(NumLocalGrids, InputGrids, GridID);
    if (InputGrid) {
      ovk_grid_params *GridParams;
      ovkCreateGridParams(&GridParams, 2);
      ovkSetGridParamName(GridParams, InputGrid->name);
      ovkSetGridParamComm(GridParams, InputGrid->comm);
      ovkSetGridParamSize(GridParams, InputGrid->global_size);
      ovk_range LocalRange;
      ovkSetRange(&LocalRange, 2, InputGrid->is, InputGrid->ie);
      ovkSetGridParamLocalRange(GridParams, &LocalRange);
      ovkCreateGridLocal(Domain, GridID, GridParams);
      ovkDestroyGridParams(&GridParams);
    } else {
      ovkCreateGridRemote(Domain, GridID);
    }
  }

  ovk_grid **Grids = malloc(NumLocalGrids*sizeof(ovk_grid *));
  iLocalGrid = 0;
  for (iGrid = 0; iGrid < 2; ++iGrid) {
    int GridID = iGrid+1;
    if (ovkRankHasGrid(Domain, GridID)) {
      ovkEditGridLocal(Domain, GridID, &Grids[iLocalGrid]);
      ++iLocalGrid;
    } else {
      ovkEditGridRemote(Domain, GridID);
    }
  }

#ifdef NOPE
  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    input_grid *InputGrid = InputGrids + iLocalGrid;
    ovk_grid *Grid = Grids[iLocalGrid];
    ovk_range LocalRange;
    ovkGetGridLocalRange(Grid, &LocalRange);
    double *XYZ[2];
    for (iDim = 0; iDim < 2; ++iDim) {
      ovkEditGridCoords(Grid, iDim, &XYZ[iDim]);
    }
    for (j = Grid->is[1]; j < Grid->ie[1]; ++j) {
      for (i = Grid->is[0]; i < Grid->ie[0]; ++i) {
        int Point[2] = {i, j};
        size_t iInputPoint = (i-Grid->is[0]) + (j-Grid->is[1]) * Grid->local_size[0];
        ovkRangeTupleToIndex(&LocalRange, OVK_GRID_LAYOUT, Point, &iPoint);
        XYZ[0][iPoint] = Grid->xyz[iInputPoint];
        YYZ[1][iPoint] = Grid->xyz[iInputPoint+Grid->local_count];
      }
    }
    for (iDim = 0; iDim < 2; ++iDim) {
      ovkReleaseGridCoords(Grid, iDim, &XYZ[iDim]);
    }
  }
#endif

  iLocalGrid = 0;
  for (iGrid = 0; iGrid < 2; ++iGrid) {
    int GridID = iGrid+1;
    if (ovkRankHasGrid(Domain, GridID)) {
      ovkReleaseGridLocal(Domain, GridID, &Grids[iLocalGrid]);
      ++iLocalGrid;
    } else {
      ovkReleaseGridRemote(Domain, GridID);
    }
  }
  free(Grids);

#ifdef NOPE
  ovk_assembly_options *Options;
  int GridIDs[2] = {1, 2};
  ovkCreateAssemblyOptions(&Options, 2, 2, GridIDs);
  ovkSetAssemblyOptionOverlappable(Options, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_TRUE);
  ovkSetAssemblyOptionTolerance(Options, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 0.1);
  ovkSetAssemblyOptionAccelQualityAdjust(Options, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 0.);
  ovkSetAssemblyOptionInferBoundaries(Options, OVK_ALL_GRIDS, OVK_TRUE);
  ovkSetAssemblyOptionCutBoundaryHoles(Options, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_TRUE);
  ovkSetAssemblyOptionOccludes(Options, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_OCCLUDES_COARSE);
  ovkSetAssemblyOptionEdgePadding(Options, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 2);
  ovkSetAssemblyOptionEdgeSmoothing(Options, OVK_ALL_GRIDS, 2);
  ovkSetAssemblyOptionConnectionType(Options, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_CUBIC);
  ovkSetAssemblyOptionFringeSize(Options, OVK_ALL_GRIDS, 2);
  ovkSetAssemblyOptionMinimizeOverlap(Options, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_TRUE);

  ovkAssemble(Domain, Options);

  ovkDestroyAssemblyOptions(&Options);
#endif

  ovkDestroyContext(&Context);

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
          int Dims[2] = {0, 0};
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
