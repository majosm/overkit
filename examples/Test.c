// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "overkit.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

int min(a, b) {
  return a < b ? a : b;
}

typedef struct {
  int id;
  char name[64];
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  int comm_dims[2];
  int comm_coords[2];
  int num_neighbors;
  int neighbor_ranks[8];
  int global_size[2];
  int local_size[2];
  int local_count;
  int is[2], ie[2];
  double *xyz;
} grid;

typedef struct {
  double *values;
} state;

void CreateInputs(int *NumGrids, grid **Grids, state **States);
void DestroyInputs(int NumGrids, grid **Grids, state **States);

int FindLocalGrid(int NumGrids, grid *Grids, int GridID);

void ExchangeTest(int argc, char **argv);
void AssembleTest(int argc, char **argv);

int main(int argc, char **argv) {

  ExchangeTest(argc, argv);
//   AssembleTest(argc, argv);

  return 0;

}

void ExchangeTest(int argc, char **argv) {

  int m, p;

  MPI_Init(&argc, &argv);

  int NumProcs, Rank;
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

  if (Rank == 0) {
    printf("Running with %i processes.\n", NumProcs);
  }

  for (p = 0; p < NumProcs; ++p) {
    if (p == Rank) {
      printf("Rank %i reporting for duty.\n", Rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  int NumGrids;
  grid *Grids;
  state *States;
  CreateInputs(&NumGrids, &Grids, &States);

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

  ovkConfigureDomain(Domain, OVK_DOMAIN_CONFIG_CONNECTIVITY);

  MPI_Comm SplitComm;
  MPI_Comm_split(MPI_COMM_WORLD, Rank % 2, Rank, &SplitComm);

  ovk_domain *OtherDomain = NULL;
  if (Rank % 2 == 1) {
    ovkCreateDomainParams(Context, &DomainParams);
    ovkSetDomainParamName(DomainParams, "OtherDomain");
    ovkSetDomainParamDimension(DomainParams, 2);
    ovkSetDomainParamComm(DomainParams, SplitComm);
    ovkCreateDomain(Context, &OtherDomain, DomainParams);
    ovkDestroyDomainParams(Context, &DomainParams);
    ovkDestroyDomain(Context, &OtherDomain);
  }

  MPI_Comm_free(&SplitComm);

  for (m = 0; m < 2; ++m) {
    p = FindLocalGrid(NumGrids, Grids, m+1);
    if (p != NumGrids) {
      ovk_grid_params *GridParams;
      ovkCreateGridParams(Domain, &GridParams);
      ovkSetGridParamID(GridParams, Grids[p].id);
      ovkSetGridParamName(GridParams, Grids[p].name);
      ovkSetGridParamComm(GridParams, Grids[p].comm);
      ovkSetGridParamGlobalSize(GridParams, Grids[p].global_size);
      ovkSetGridParamLocalStart(GridParams, Grids[p].is);
      ovkSetGridParamLocalEnd(GridParams, Grids[p].ie);
      ovkSetGridParamNeighborRanks(GridParams, Grids[p].num_neighbors, Grids[p].neighbor_ranks);
      ovkCreateGridLocal(Domain, NULL, GridParams);
      ovkDestroyGridParams(Domain, &GridParams);
    } else {
      ovkCreateGridRemote(Domain, NULL);
    }
  }

#ifdef NOPE

//   ovk_grid *EditGrid[NumGrids];
//   for (m = 0; m < 2; ++m) {
//     int GridID = m+1;
//     p = FindLocalGrid(NumGrids, Grids, GridID);
//     if (p != NumGrids) {
//       ovkEditGridLocal(Domain, GridID, EditGrid[p]);
//     } else {
//       ovkEditGridRemote(Domain, GridID);
//     }
//   }

//   for (m = 0; m < NumGrids; ++m) {
//     ovkEditGridCoords(EditGrid[m], Coords);
//     // Do stuff
//     ovkReleaseGridCoords(EditGrid[m], Coords);
//   }

//   for (m = 0; m < 2; ++m) {
//     int GridID = m+1;
//     p = FindLocalGrid(NumGrids, Grids, GridID);
//     if (p != NumGrids) {
//       ovkReleaseGridLocal(Domain, GridID, EditGrids[p]);
//     } else {
//       ovkReleaseGridRemote(Domain, GridID);
//     }
//   }

  for (n = 0; n < 2; ++n) {
    q = FindLocalGrid(NumGrids, Grids, n+1);
    for (m = 0; m < 2; ++m) {
      p = FindLocalGrid(NumGrids, Grids, m+1);
      ovkEditConnectivity(Domain, m, n, Connectivity);
      ovkEditConnectivityProperties(Connectivity, ConnectivityProperties);
      ovkSetConnectivityPropertyMaxDonorSize(Connectivity, 1);
      ovkReleaseConnectivityProperties(Connectivity, ConnectivityProperties);
      // Donors
      if (p != NumGrids) {
        if ((Grids[p].id == 1 && Grids[p].ie[0] == Grids[p].global_size[0]-1) ||
          (Grids[p].id == 2 && Grids[p].is[0] == 0)) {
          ovkResizeConnectivity(Connectivity, Grids[p].ie[1]-Grids[p].is[1]);
        }
      }
      // Receivers (same as donor points in this example)
      if (q != NumGrids) {
        if ((Grids[q].id == 1 && Grids[q].ie[0] == Grids[q].global_size[0]-1) ||
          (Grids[q].id == 2 && Grids[q].is[0] == 0)) {
        ovkResizeConnectivity(Connectivity, Grids[q].ie[1]-Grids[q].is[1]);
      }
      int *ConnectionID;
      int *DonorCellStart[2], *DonorCellEnd[2];
      double *DonorCellCoords[2];
      double *DonorInterpCoefs[2][1];
      int *DestinationGridID;
      int *DestinationReceiverID;
      int *ReceiverPoint[2];
      int *SourceGridID;
      int *SourceDonorID;
      // TODO: Add connectivity data from below here
      ovkReleaseConnectivity(Domain, Connectivity);
    }
  }

//   for (m = 0; m < NumGrids; ++m) {
//     ovkEditDonors(Domain, Grids[m].id, &Donors);
//     ovkEditDonorID(Donors, &DonorID);
//     ovkEditDonorCellStart(Donors, &DonorCellStart);
//     ovkEditDonorCellEnd(Donors, &DonorCellEnd);
//     ovkEditDonorCellCoords(Donors, &DonorCellCoords);
//     ovkEditDonorInterpCoefs(Donors, &DonorInterpCoefs);
//     ovkEditDonorDestinationGridID(Donors, &DestinationGridID);
//     ovkEditDonorDestinationReceiverID(Donors, &DestinationReceiverID);
//     ovkEditReceivers(Domain, Grids[m].id, &Receivers);
//     ovkEditReceiverID(Receivers, &ReceiverID);
//     ovkEditReceiverPoint(Receivers, &ReceiverPoint);
//     ovkEditReceiverSourceGridID(Receivers, &SourceGridID);
//     ovkEditReceiverSourceDonorID(Receivers, &SourceDonorID);
//     if ((Grids[m].id == 1 && Grids[m].ie[0] == Grids[m].global_size[0]-1) ||
//       (Grids[m].id == 2 && Grids[m].is[0] == 0)) {
//       int OtherGridID;
//       switch (Grids[m].id) {
//       case 1:
//         i = Grids[m].ie[0]-1;
//         OtherGridID = 2;
//         break;
//       case 2:
//         i = 0;
//         OtherGridID = 1;
//         break;
//       }
//       for (j = Grids[m].is[1]; j < Grids[m].ie[1]; ++j) {
//         DonorID[l] = j;
//         DonorCellStart[0][l] = i;
//         DonorCellStart[1][l] = j;
//         DonorCellEnd[0][l] = i;
//         DonorCellEnd[1][l] = j;
//         DonorCellCoords[0][l] = 0.;
//         DonorCellCoords[1][l] = 0.;
//         DonorInterpCoefs[0][0][l] = 1.;
//         DonorInterpCoefs[1][0][l] = 1.;
//         DestinationGridID[l] = OtherGridID;
//         DestinationReceiverID[l] = j;
//         ReceiverID[l] = j;
//         ReceiverPoint[0][l] = i;
//         ReceiverPoint[1][l] = j;
//         SourceGridID[l] = OtherGridID;
//         SourceDonorID[l] = j;
//         ++l;
//       }
//     }
//     ovkReleaseDonorID(Donors, &DonorID);
//     ovkReleaseDonorCellStart(Donors, &DonorCellStart);
//     ovkReleaseDonorCellEnd(Donors, &DonorCellEnd);
//     ovkReleaseDonorCellCoords(Donors, &DonorCellCoords);
//     ovkReleaseDonorInterpCoefs(Donors, &DonorInterpCoefs);
//     ovkReleaseDonorDestinationGridID(Donors, &DestinationGridID);
//     ovkReleaseDonorDestinationReceiverID(Donors, &DestinationReceiverID);
//     ovkReleaseDonors(Domain, &Donors);
//     ovkReleaseReceiverID(Receivers, &ReceiverID);
//     ovkReleaseReceiverPoints(Receivers, &ReceiverPoints);
//     ovkReleaseReceiverSourceGridID(Receivers, &SourceGridID);
//     ovkReleaseReceiverSourceDonorID(Receivers, &SourceDonorID);
//     ovkReleaseReceivers(Domain, &Receivers);
//   }

  // Exchange variant #1 -- simple two-way exchange

  ovk_request *Requests = malloc(2*NumGrids*sizeof(ovk_request));
  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      ovkExchange(Domain, Grids[n].id, m+1, 1, OVK_DOUBLE, OVK_COLLECT_INTERPOLATE,
        OVK_DISPERSE_OVERWRITE, States[n]->values, 1, Requests[l]);
      ++l;
    }
  }
  ovkWaitAll(Domain, Requests);

  // Exchange variant #2 -- one-way exchange

  ovk_request *Requests = malloc(2*NumGrids*sizeof(ovk_request));
  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    if (Grids[n].id == 1) {
      ovkDonate(Domain, Grids[n].id, 2, 1, OVK_DOUBLE, OVK_COLLECT_INTERPOLATE,
        OVK_DISPERSE_OVERWRITE, States[n]->values, 1, Requests[l]);
    } else {
      ovkReceive(Domain, Grids[n].id, 1, 1, OVK_DOUBLE, OVK_COLLECT_INTERPOLATE,
        OVK_DISPERSE_OVERWRITE, States[n]->values, 1, Requests[l]);
    }
    ++l;
  }
  ovkWaitAll(Domain, Requests);

  // Exchange variant #3 -- manual collect-send-receive-disperse

  double **SendBuffers = malloc(2*NumGrids*sizeof(double *));
  double **ReceiveBuffers = malloc(2*NumGrids*sizeof(double *));
  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      size_t NumDonors = ovkGetDonorCount(Domain, Grids[n].id, m+1);
      size_t NumReceivers = ovkGetReceiverCount(Domain, Grids[n].id, m+1);
      SendBuffers[l] = malloc(NumDonors*sizeof(double));
      ReceiveBuffers[l] = malloc(NumReceivers*sizeof(double));
      ++l;
    }
  }

  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      ovkCollect(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, OVK_COLLECT_INTERPOLATE,
        States[n]->values, SendBuffers[l]);
      ++l;
    }
  }

  ovk_request *Requests = malloc(4*NumGrids*sizeof(ovk_request));
  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      ovkTake(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, ReceiveBuffers[l], 1, Requests[2*l]);
      ovkGive(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, SendBuffers[l], 1, Requests[2*l+1]);
      l += 2;
    }
  }
  ovkWaitAll(Domain, 4*NumGrids, Requests);

  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      ovkDisperse(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, OVK_DISPERSE_OVERWRITE,
        ReceiveBuffers[l], State[n]->values);
      ++l;
    }
  }

  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      free(SendBuffers[l]);
      free(ReceiveBuffers[l]);
      ++l;
    }
  }
  free(SendBuffers);
  free(ReceiverBuffers);

  // Exchange variant #4 -- manual collect-send-receive-disperse with interleaved receive/disperse

  double **SendBuffers = malloc(2*NumGrids*sizeof(double *));
  double **ReceiveBuffers = malloc(2*NumGrids*sizeof(double *));
  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      size_t NumDonors = ovkGetDonorCount(Domain, Grids[n].id, m+1);
      size_t NumReceivers = ovkGetReceiverCount(Domain, Grids[n].id, m+1);
      SendBuffers[l] = malloc(NumDonors*sizeof(double));
      ReceiveBuffers[l] = malloc(NumReceivers*sizeof(double));
      ++l;
    }
  }

  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      ovkCollect(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, OVK_COLLECT_INTERPOLATE,
        States[n]->values, SendBuffers[l]);
      ++l;
    }
  }

  ovk_request *SendRequests = malloc(2*NumGrids*sizeof(ovk_request));
  ovk_request *ReceiveRequests = malloc(2*NumGrids*sizeof(ovk_request));
  l = 0;
  for (n = 0; n < NumGrids; ++n) {
    for (m = 0; m < 2; ++m) {
      ovkTake(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, ReceiveBuffers[l], 1, ReceiveRequests[l]);
      ovkGive(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, SendBuffers[l], 1, SendRequests[l]);
      ++l;
    }
  }

  ovkWaitAny(Domain, 2*NumGrids, ReceiveRequests, &l);
  while (l != OVK_UNDEFINED_INDEX) {
    ovkWait(Domain, SendRequests[l]);
    m = l % 2;
    n = l/2 % NumGrids;
    ovkDisperse(Domain, Grids[n].id, m+1, OVK_DOUBLE, 1, OVK_DISPERSE_OVERWRITE, ReceiveBuffers[l],
      State[n]->values);
    ovkWaitAny(Domain, 2*NumGrids, ReceiveRequests, &l);
  }

  l = 0;
  for (n = 0; n < NumGrids; ++n) {
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

  DestroyInputs(NumGrids, &Grids, &States);

  MPI_Finalize();

}

void AssembleTest(int argc, char **argv) {

  int p;

  MPI_Init(&argc, &argv);

  int NumProcs, Rank;
  MPI_Comm_size(MPI_COMM_WORLD, &NumProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

  if (Rank == 0) {
    printf("Running with %i processes.\n", NumProcs);
  }

  for (p = 0; p < NumProcs; ++p) {
    if (p == Rank) {
      printf("Rank %i reporting for duty.\n", Rank);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  int NumGrids;
  grid *Grids;
  state *States;
  CreateInputs(&NumGrids, &Grids, &States);

#ifdef NOPE

  ovk_domain_params *DomainParams;
  ovkCreateDomainParam(&DomainParams, 2);
  ovkSetDomainParamComm(DomainParams, MPI_COMM_WORLD);
  ovkSetDomainParamLogLevel(DomainParams, OVK_LOG_ALL);
  ovkSetDomainParamErrorHandler(DomainParams, OVK_ERROR_HANDLER_ABORT);

  ovk_domain *Domain;
  ovkCreateDomain(&Domain, DomainParams);

  ovkDestroyDomainParams(&DomainParams);

  ovkConfigureDomain(Domain, OVK_DOMAIN_CONFIG_GEOMETRY | OVK_DOMAIN_CONFIG_OVERLAP |
    OVK_DOMAIN_CONFIG_CONNECTIVITY);

  for (m = 0; m < NumGrids; ++m) {
    ovk_grid_params *GridParams;
    ovkCreateGridParams(&GridParams, 2);
    ovkSetGridParamComm(GridParams, Grids[m].comm);
    ovkSetGridParamGlobalSize(GridParams, Grids[m].global_size);
    ovkSetGridParamLocalStart(GridParams, Grids[m].is);
    ovkSetGridParamLocalEnd(GridParams, Grids[m].ie);
    ovkSetGridParamNeighborRanks(GridParams, Grids[m].num_neighbors, Grids[m].neighbor_ranks);
    ovkCreateGrid(Domain, OVK_ID_MANUAL, &Grids[m].id, GridParams);
    ovkDestroyGridParams(&GridParams);
  }

  ovk_domain_properties *Properties;
  ovkEditDomainProperties(Domain, &Properties);
  ovkSetDomainPropertyInferBoundaries(Properties, OVK_ALL_GRIDS, true);
  ovkSetDomainPropertyOverlappable(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, true);
  ovkSetDomainPropertyOverlapTolerance(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 1.e-10);
  ovkSetDomainPropertyBoundaryHoleCutting(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, false);
  ovkSetDomainPropertyOverlapHoleCutting(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, false);
  ovkSetDomainPropertyConnectionType(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_FRINGE);
  ovkSetDomainPropertyFringeSize(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 1);
  ovkSetDomainPropertyFringePadding(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 0);
  ovkReleaseDomainProperties(Domain, &Properties);

  for (m = 0; m < NumGrids; ++m) {
    ovk_grid *Grid;
    ovkEditGrid(Domain, Grids[m].id, &Grid);
    ovk_cart Cart;
    ovkGetGridCart(Grid, &Cart);
    double *XYZ[2];
    ovkEditGridCoords(Grid, &XYZ);
    for (j = Grids[m].is[1]; j < Grids[m].ie[1]; ++j) {
      for (i = Grids[m].is[0]; i < Grids[m].ie[0]; ++i) {
        int Point[2] = {i,j};
        l = (i-Grids[m].is[0]) + (j-Grids[m].is[1]) * Grids[m].local_size[0];
        p = ovkCartTupleToIndex(Cart, Point);
        XYZ[0][p] = Grids[m].xyz[l];
        YYZ[1][p] = Grids[m].xyz[l+Grids[m].local_count];
      }
    }
    ovkReleaseGridCoords(Grid, &XYZ);
    ovkReleaseGrid(Domain, &Grid);
  }

  ovkDetectOverlap(Domain);
  ovkCutHoles(Domain);
  ovkGenerateConnectivity(Domain);

  ovkDestroyDomain(&Domain);

#endif

  DestroyInputs(NumGrids, &Grids, &States);

  MPI_Finalize();

}

void CreateInputs(int *NumGrids, grid **Grids, state **States) {

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

  *NumGrids = 0;
  for (m = 0; m < 2; ++m) {
    if (Rank >= GridToProcMap[2*m] && Rank < GridToProcMap[2*m+1]) {
      *NumGrids += 1;
    }
  }

  *Grids = malloc((*NumGrids)*sizeof(**Grids));

  grid *Grids_ = *Grids;

  for (m = 0, n = 0; m < 2; ++m) {
    if (Rank >= GridToProcMap[2*m] && Rank < GridToProcMap[2*m+1]) {
      Grids_[n].id = m+1;
      ++n;
    }
  }

  for (n = 0; n < *NumGrids; ++n) {
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
      for (n = 0; n < *NumGrids; ++n) {
        grid *Grid = &Grids_[n];
        if (Grid->id == m+1) {
          int Dims[2] = {0, 0};
          MPI_Dims_create(NumProcsPerGrid[m], 2, Dims);
          int Periods[2] = {0, 0};
          MPI_Cart_create(GridComm, 2, Dims, Periods, 0, &GridCartComm);
          Grid->comm = GridCartComm;
          MPI_Comm_size(Grid->comm, &Grid->comm_size);
          MPI_Comm_rank(Grid->comm, &Grid->comm_rank);
          MPI_Cart_get(Grid->comm, 2, Grid->comm_dims, Periods, Grid->comm_coords);
          Grid->num_neighbors = 0;
          for (j = -1; j <= 1; ++j) {
            for (i = -1; i <= 1; ++i) {
              int Coords[2] = {Grid->comm_coords[0]+i,Grid->comm_coords[1]+j};
              bool ValidNeighbor = i != 0 || j != 0;
              ValidNeighbor = ValidNeighbor && Coords[0] >= 0 && Coords[0] < Grid->comm_dims[0];
              ValidNeighbor = ValidNeighbor && Coords[1] >= 0 && Coords[1] < Grid->comm_dims[1];
              if (ValidNeighbor) {
                MPI_Cart_rank(Grid->comm, Coords, &Grid->neighbor_ranks[Grid->num_neighbors]);
                ++Grid->num_neighbors;
              }
            }
          }
        }
      }
      MPI_Comm_free(&GridComm);
    }
    MPI_Group_free(&GridGroup);
    free(Ranks);
  }

  MPI_Group_free(&WorldGroup);

  for (n = 0; n < *NumGrids; ++n) {
    grid *Grid = &Grids_[n];
    Grid->global_size[0] = 100;
    Grid->global_size[1] = 100;
    for (d = 0; d < 2; ++d) {
      int N = Grid->global_size[d]/Grid->comm_dims[d];
      int R = Grid->global_size[d] - Grid->comm_dims[d] * N;
      Grid->is[d] = N*Grid->comm_coords[d] + min(R, Grid->comm_coords[d]);
      Grid->ie[d] = N*(Grid->comm_coords[d]+1) + min(R, (Grid->comm_coords[d]+1));
      Grid->local_size[d] = Grid->ie[d] - Grid->is[d];
    }
    Grid->local_count = Grid->local_size[0] * Grid->local_size[1];
  }

  for (n = 0; n < *NumGrids; ++n) {
    grid *Grid = &Grids_[n];
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
    for (j = Grid->is[1]; j < Grid->is[1]; ++j) {
      for (i = Grid->is[0]; i < Grid->ie[0]; ++i) {
        l = (i-Grid->is[0]) + (j-Grid->is[1]) * Grid->local_size[0];
        double u = (double)i/(double)(Grid->global_size[0]-1);
        double v = (double)j/(double)(Grid->global_size[1]-1);
        Grid->xyz[l] = (1.-u)*xs + u*xe;
        Grid->xyz[l+Grid->local_count] = (1.-v)*ys + v*ye;
      }
    }
  }

  *States = malloc((*NumGrids)*sizeof(**States));

  state *States_ = *States;

  for (n = 0; n < *NumGrids; ++n) {
    grid *Grid = &Grids_[n];
    state *State = &States_[n];
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
    for (j = Grid->is[1]; j < Grid->is[1]; ++j) {
      for (i = Grid->is[0]; i < Grid->ie[0]; ++i) {
        l = (i-Grid->is[0]) + (j-Grid->is[1]) * Grid->local_size[0];
        double v = (double)j/(double)(Grid->global_size[1]-1);
        State->values[l] = (1.-v)*fs + v*fe;
      }
    }
  }

}

void DestroyInputs(int NumGrids, grid **Grids, state **States) {

  int n;

  for (n = 0; n < NumGrids; ++n) {
    grid *Grid = &(*Grids)[n];
    free(Grid->xyz);
    Grid->xyz = NULL;
    MPI_Comm_free(&Grid->comm);
  }
  free(*Grids);
  *Grids = NULL;

  for (n = 0; n < NumGrids; ++n) {
    state *State = &(*States)[n];
    free(State->values);
    State->values = NULL;
  }
  free(*States);
  *States = NULL;

}

int FindLocalGrid(int NumGrids, grid *Grids, int GridID) {

  int LocalIndex;

  for (LocalIndex = 0; LocalIndex < NumGrids; ++LocalIndex) {
    if (Grids[LocalIndex].id == GridID) break;
  }

  return LocalIndex;

}
