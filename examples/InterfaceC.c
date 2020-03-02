// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <overkit.h>

#include "examples/Common.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define min(a, b) ovk_min(a, b)
#define max(a, b) ovk_max(a, b)

#define CreateCartesianDecompDims examples_CreateCartesianDecompDims
#define CartesianDecomp examples_CartesianDecomp
#define command_args examples_command_args
#define command_args_parser examples_command_args_parser
#define command_args_error examples_command_args_error
#define CreateCommandArgsParser examples_CreateCommandArgsParser
#define DestroyCommandArgsParser examples_DestroyCommandArgsParser
#define SetCommandArgsParserHelpUsage examples_SetCommandArgsParserHelpUsage
#define SetCommandArgsParserHelpDescription examples_SetCommandArgsParserHelpDescription
#define AddCommandArgsParserOption examples_AddCommandArgsParserOption
#define COMMAND_ARGS_VALUE_TYPE_BOOL EXAMPLES_COMMAND_ARGS_VALUE_TYPE_BOOL
#define COMMAND_ARGS_VALUE_TYPE_INT EXAMPLES_COMMAND_ARGS_VALUE_TYPE_INT
#define ParseCommandArgs examples_ParseCommandArgs
#define DestroyCommandArgs examples_DestroyCommandArgs
#define GetCommandOptionIfPresent examples_GetCommandOptionIfPresent
#ifdef OVK_HAVE_XDMF
#define xdmf examples_xdmf
#define xdmf_grid_meta examples_xdmf_grid_meta
#define xdmf_attribute_meta examples_xdmf_attribute_meta
#define xdmf_attribute_type examples_xdmf_attribute_type
#define XDMF_ATTRIBUTE_TYPE_INT EXAMPLES_XDMF_ATTRIBUTE_TYPE_INT
#define XDMF_ATTRIBUTE_TYPE_DOUBLE EXAMPLES_XDMF_ATTRIBUTE_TYPE_DOUBLE
#define xdmf_error examples_xdmf_error
#define CreateXDMFGridMeta examples_CreateXDMFGridMeta
#define CreateXDMFAttributeMeta examples_CreateXDMFAttributeMeta
#define CreateXDMF examples_CreateXDMF
#define OpenXDMF examples_OpenXDMF
#define CloseXDMF examples_CloseXDMF
#define WriteXDMFGeometry examples_WriteXDMFGeometry
#define WriteXDMFAttribute examples_WriteXDMFAttribute
#endif

static int GetCommandLineArguments(int argc, char **argv, bool *Help, int *N);
static int Interface(int N);

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int WorldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

  int Error;

  bool Help;
  int N;

  Error = GetCommandLineArguments(argc, argv, &Help, &N);
  if (Error) goto check_error;

  if (!Help) {
    Error = Interface(N);
    if (Error) goto check_error;
  }

  check_error:
    if (Error) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (WorldRank == 0) {
        fprintf(stderr, "Error occurred.\n"); fflush(stderr);
      }
    }

  MPI_Finalize();

  return 0;

}

int GetCommandLineArguments(int argc, char **argv, bool *Help, int *N) {

  int WorldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

  int Error = 0;

  command_args_error CommandArgsError;

  command_args_parser *CommandArgsParser;
  CreateCommandArgsParser(&CommandArgsParser, WorldRank == 0);
  SetCommandArgsParserHelpUsage(CommandArgsParser, "InterfaceC [<options> ...]");
  SetCommandArgsParserHelpDescription(CommandArgsParser, "Generates an overset mesh consisting of "
    "two grids overlapping along an interface.");
  AddCommandArgsParserOption(CommandArgsParser, "size", 'N', COMMAND_ARGS_VALUE_TYPE_INT,
    "Characteristic size of grids [ Default: 64 ]");

  command_args *CommandArgs;
  ParseCommandArgs(CommandArgsParser, argc, argv, &CommandArgs, &CommandArgsError);
  if (CommandArgsError) {
    Error = 1;
    goto cleanup1;
  }

  *Help = false;
  GetCommandOptionIfPresent(CommandArgs, "help", COMMAND_ARGS_VALUE_TYPE_BOOL, Help,
    &CommandArgsError);
  if (CommandArgsError) {
    Error = 1;
    goto cleanup2;
  }

  *N = 64;
  GetCommandOptionIfPresent(CommandArgs, "size", COMMAND_ARGS_VALUE_TYPE_INT, N, &CommandArgsError);
  if (CommandArgsError) {
    Error = 1;
    goto cleanup2;
  }

  cleanup2:
    DestroyCommandArgs(&CommandArgs);
  cleanup1:
    DestroyCommandArgsParser(&CommandArgsParser);

  return Error;

}

typedef struct {
  MPI_Comm Comm;
  int Size[3];
  int LocalRange[6];
  long long NumLocalPoints;
  int ExtendedRange[6];
  long long NumExtendedPoints;
} grid_data;

static void CreateGridData(grid_data *Data) {

  Data->Comm = MPI_COMM_NULL;
  memset(Data->Size, 0, 3);
  Data->Size[2] = 1;
  memset(Data->LocalRange, 0, 6);
  Data->LocalRange[5] = 1;
  Data->NumLocalPoints = 0;
  memset(Data->ExtendedRange, 0, 6);
  Data->ExtendedRange[5] = 1;
  Data->NumExtendedPoints = 0;

}

static void DestroyGridData(grid_data *Data) {

  if (Data->Comm != MPI_COMM_NULL) {
    MPI_Comm_free(&Data->Comm);
  }

}

#ifdef OVK_HAVE_XDMF
static int *CreateOutputState(const ovk_domain *Domain, int ConnectivityComponentID, int GridID) {

  const ovk_connectivity_component *ConnectivityComponent;
  ovkGetComponent(Domain, ConnectivityComponentID, OVK_COMPONENT_TYPE_CONNECTIVITY,
    &ConnectivityComponent);
  const ovk_grid *Grid;
  ovkGetGrid(Domain, GridID, &Grid);
  long long NumLocalPoints;
  ovkGetGridLocalCount(Grid, &NumLocalPoints);
  int LocalRange[6];
  ovkGetGridLocalRange(Grid, LocalRange, LocalRange+3);

  int *OutputState = malloc(NumLocalPoints*sizeof(int));

  long long iPoint;
  for (iPoint = 0; iPoint < NumLocalPoints; ++iPoint) {
    OutputState[iPoint] = 1;
  }
  int NumLocalConnectivityNs = ovkLocalConnectivityNCount(ConnectivityComponent);
  int *LocalConnectivityNMGridIDs = malloc(NumLocalConnectivityNs*sizeof(int));
  int *LocalConnectivityNNGridIDs = malloc(NumLocalConnectivityNs*sizeof(int));
  ovkGetLocalConnectivityNIDs(ConnectivityComponent, LocalConnectivityNMGridIDs,
    LocalConnectivityNNGridIDs);
  long long iConnectivityN;
  for (iConnectivityN = 0; iConnectivityN < NumLocalConnectivityNs; ++iConnectivityN) {
    int MGridID = LocalConnectivityNMGridIDs[iConnectivityN];
    int NGridID = LocalConnectivityNNGridIDs[iConnectivityN];
    if (NGridID != GridID) continue;
    const ovk_connectivity_n *ConnectivityN;
    ovkGetConnectivityN(ConnectivityComponent, MGridID, NGridID, &ConnectivityN);
    long long NumReceivers = ovkGetConnectivityNSize(ConnectivityN);
    const int *Points[3];
    int iDim;
    for (iDim = 0; iDim < 3; ++iDim) {
      ovkGetConnectivityNPoints(ConnectivityN, iDim, &Points[iDim]);
    }
    long long iReceiver;
    for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      int i = Points[0][iReceiver];
      int j = Points[1][iReceiver];
      int k = Points[2][iReceiver];
      long long l = (i-LocalRange[0]) + (LocalRange[3]-LocalRange[0])*((j-LocalRange[1]) +
        (LocalRange[4]-LocalRange[1])*(k-LocalRange[2]));
      OutputState[l] = -MGridID;
    }
  }
  free(LocalConnectivityNMGridIDs);
  free(LocalConnectivityNNGridIDs);

  return OutputState;

}
#endif

static int Interface(int N) {

  int iDim, iCoef;
  int i, j, l;
  long long iDonor, iReceiver;

  int NumWorldProcs, WorldRank;
  MPI_Comm_size(MPI_COMM_WORLD, &NumWorldProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

  ovk_context_params *ContextParams;
  ovkCreateContextParams(&ContextParams);
  ovkSetContextParamComm(ContextParams, MPI_COMM_WORLD);
  ovkSetContextParamStatusLoggingThreshold(ContextParams, 4);

  ovk_context *Context;
  ovk_error CreateError;
  ovkCreateContext(&Context, &ContextParams, &CreateError);
  if (CreateError != OVK_ERROR_NONE) return 1;

  ovk_shared_context *SharedContext;
  ovkShareContext(&Context, &SharedContext);

  ovk_domain_params *DomainParams;
  ovkCreateDomainParams(&DomainParams);
  ovkSetDomainParamDimension(DomainParams, 2);
  ovkSetDomainParamComm(DomainParams, MPI_COMM_WORLD);

  ovk_domain *Domain;
  ovkCreateDomain(&Domain, SharedContext, &DomainParams);

  int Size[] = {N,N,1};

  int GridIDs[] = {1, 2};
  ovk_grid_params *MaybeGridParams[] = {NULL, NULL};

  bool LeftIsLocal = WorldRank < max(NumWorldProcs/2, 1);
  bool RightIsLocal = WorldRank >= NumWorldProcs/2;

  int LeftSize[] = {(Size[0]+2)/2, Size[1], Size[2]};
  int RightSize[] = {Size[0]+2-(Size[0]+2)/2, Size[1], Size[2]};

  grid_data LeftData, RightData;
  CreateGridData(&LeftData);
  CreateGridData(&RightData);

  if (LeftIsLocal) {
    grid_data *Data = &LeftData;
    MPI_Comm TempComm;
    MPI_Comm_split(MPI_COMM_WORLD, 0, WorldRank, &TempComm);
    int NumGridProcs;
    MPI_Comm_size(TempComm, &NumGridProcs);
    int CartDims[] = {0,0,1};
    CreateCartesianDecompDims(NumGridProcs, 2, CartDims);
    int CartPeriods[] = {0,0,0};
    MPI_Cart_create(TempComm, 2, CartDims, CartPeriods, 1, &Data->Comm);
    MPI_Comm_free(&TempComm);
    Data->Size[0] = LeftSize[0];
    Data->Size[1] = LeftSize[1];
    Data->Size[2] = LeftSize[2];
    CartesianDecomp(2, Data->Size, Data->Comm, Data->LocalRange);
    Data->NumLocalPoints =
      (long long)(Data->LocalRange[3] - Data->LocalRange[0]) *
      (long long)(Data->LocalRange[4] - Data->LocalRange[1]) *
      (long long)(Data->LocalRange[5] - Data->LocalRange[2]);
    // Pretend we have a halo
    for (int iDim = 0; iDim < 3; ++iDim) {
      Data->ExtendedRange[iDim] = Data->LocalRange[iDim];
      Data->ExtendedRange[3+iDim] = Data->LocalRange[3+iDim];
      if (Data->LocalRange[iDim] > 0) --Data->ExtendedRange[iDim];
      if (Data->LocalRange[3+iDim] < Data->Size[iDim]) ++Data->ExtendedRange[3+iDim];
    }
    Data->NumExtendedPoints =
      (long long)(Data->ExtendedRange[3] - Data->ExtendedRange[0]) *
      (long long)(Data->ExtendedRange[4] - Data->ExtendedRange[1]) *
      (long long)(Data->ExtendedRange[5] - Data->ExtendedRange[2]);
  } else {
    MPI_Comm DummyComm;
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, WorldRank, &DummyComm);
  }

  if (RightIsLocal) {
    grid_data *Data = &RightData;
    MPI_Comm TempComm;
    MPI_Comm_split(MPI_COMM_WORLD, 0, WorldRank, &TempComm);
    int NumGridProcs;
    MPI_Comm_size(TempComm, &NumGridProcs);
    int CartDims[] = {0,0,1};
    CreateCartesianDecompDims(NumGridProcs, 2, CartDims);
    int CartPeriods[] = {0,0,0};
    MPI_Cart_create(TempComm, 2, CartDims, CartPeriods, 1, &Data->Comm);
    MPI_Comm_free(&TempComm);
    Data->Size[0] = RightSize[0];
    Data->Size[1] = RightSize[1];
    Data->Size[2] = RightSize[2];
    CartesianDecomp(2, Data->Size, Data->Comm, Data->LocalRange);
    Data->NumLocalPoints =
      (long long)(Data->LocalRange[3] - Data->LocalRange[0]) *
      (long long)(Data->LocalRange[4] - Data->LocalRange[1]) *
      (long long)(Data->LocalRange[5] - Data->LocalRange[2]);
    // Pretend we have a halo
    for (int iDim = 0; iDim < 3; ++iDim) {
      Data->ExtendedRange[iDim] = Data->LocalRange[iDim];
      Data->ExtendedRange[3+iDim] = Data->LocalRange[3+iDim];
      if (Data->LocalRange[iDim] > 0) --Data->ExtendedRange[iDim];
      if (Data->LocalRange[3+iDim] < Data->Size[iDim]) ++Data->ExtendedRange[3+iDim];
    }
    Data->NumExtendedPoints =
      (long long)(Data->ExtendedRange[3] - Data->ExtendedRange[0]) *
      (long long)(Data->ExtendedRange[4] - Data->ExtendedRange[1]) *
      (long long)(Data->ExtendedRange[5] - Data->ExtendedRange[2]);
  } else {
    MPI_Comm DummyComm;
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, WorldRank, &DummyComm);
  }

  if (LeftIsLocal) {
    ovkCreateGridParams(&MaybeGridParams[0]);
    ovk_grid_params *GridParams = MaybeGridParams[0];
    ovkSetGridParamName(GridParams, "Left");
    ovkSetGridParamDimension(GridParams, 2);
    ovkSetGridParamComm(GridParams, LeftData.Comm);
    const int Zero[3] = {0,0,0};
    ovkSetGridParamGlobalRange(GridParams, Zero, LeftData.Size);
    ovkSetGridParamLocalRange(GridParams, LeftData.LocalRange, LeftData.LocalRange+3);
  }

  if (RightIsLocal) {
    ovkCreateGridParams(&MaybeGridParams[1]);
    ovk_grid_params *GridParams = MaybeGridParams[1];
    ovkSetGridParamName(GridParams, "Right");
    ovkSetGridParamDimension(GridParams, 2);
    ovkSetGridParamComm(GridParams, RightData.Comm);
    const int Zero[3] = {0,0,0};
    ovkSetGridParamGlobalRange(GridParams, Zero, RightData.Size);
    ovkSetGridParamLocalRange(GridParams, RightData.LocalRange, RightData.LocalRange+3);
  }

  ovkCreateGrids(Domain, 2, GridIDs, MaybeGridParams);

  const int CONNECTIVITY_ID = 1;
  ovkCreateComponent(Domain, CONNECTIVITY_ID, OVK_COMPONENT_TYPE_CONNECTIVITY, NULL);

  {

    ovk_connectivity_component *ConnectivityComponent;
    ovkEditComponent(Domain, CONNECTIVITY_ID, OVK_COMPONENT_TYPE_CONNECTIVITY,
      &ConnectivityComponent);

    int MGridIDs[] = {1, 2};
    int NGridIDs[] = {2, 1};

    ovkCreateConnectivities(ConnectivityComponent, 2, MGridIDs, NGridIDs);

    if (LeftIsLocal) {

      const ovk_grid *Grid;
      ovkGetGrid(Domain, 1, &Grid);

      int GlobalRange[6];
      ovkGetGridGlobalRange(Grid, GlobalRange, GlobalRange+3);

      int LocalRange[6];
      ovkGetGridLocalRange(Grid, LocalRange, LocalRange+3);

      bool HasInterface = LocalRange[3] == GlobalRange[3];

      ovk_connectivity_m *ConnectivityM;
      ovkEditConnectivityM(ConnectivityComponent, 1, 2, &ConnectivityM);

      long long NumDonors = HasInterface ? (LocalRange[4]-LocalRange[1]) : 0;
      ovkResizeConnectivityM(ConnectivityM, NumDonors, 1);

      int *Extents[2][2];
      double *Coords[2];
      double *InterpCoefs[2][1];
      int *Destinations[2];
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityMExtents(ConnectivityM, iDim, &Extents[0][iDim], &Extents[1][iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityMCoords(ConnectivityM, iDim, &Coords[iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        for (iCoef = 0; iCoef < 1; ++iCoef) {
          ovkEditConnectivityMInterpCoefs(ConnectivityM, iDim, iCoef, &InterpCoefs[iDim][iCoef]);
        }
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityMDestinations(ConnectivityM, iDim, &Destinations[iDim]);
      }

      if (HasInterface) {
        iDonor = 0;
        for (j = LocalRange[1]; j < LocalRange[4]; ++j) {
          Extents[0][0][iDonor] = GlobalRange[3]-2;
          Extents[0][1][iDonor] = j;
          Extents[1][0][iDonor] = Extents[0][0][iDonor]+1;
          Extents[1][1][iDonor] = Extents[0][1][iDonor]+1;
          Coords[0][iDonor] = 0.;
          Coords[1][iDonor] = 0.;
          InterpCoefs[0][0][iDonor] = 1.;
          InterpCoefs[1][0][iDonor] = 1.;
          Destinations[0][iDonor] = 0;
          Destinations[1][iDonor] = j;
          ++iDonor;
        }
      }

      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityMExtents(ConnectivityM, iDim, &Extents[0][iDim], &Extents[1][iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityMCoords(ConnectivityM, iDim, &Coords[iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        for (iCoef = 0; iCoef < 1; ++iCoef) {
          ovkRestoreConnectivityMInterpCoefs(ConnectivityM, iDim, iCoef,
            &InterpCoefs[iDim][iCoef]);
        }
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityMDestinations(ConnectivityM, iDim, &Destinations[iDim]);
      }

      ovkRestoreConnectivityM(ConnectivityComponent, 1, 2, &ConnectivityM);

      ovk_connectivity_n *ConnectivityN;
      ovkEditConnectivityN(ConnectivityComponent, 2, 1, &ConnectivityN);

      long long NumReceivers = HasInterface ? (LocalRange[4]-LocalRange[1]) : 0;
      ovkResizeConnectivityN(ConnectivityN, NumReceivers);

      int *Points[2];
      int *Sources[2];
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityNPoints(ConnectivityN, iDim, &Points[iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityNSources(ConnectivityN, iDim, &Sources[iDim]);
      }

      if (HasInterface) {
        iReceiver = 0;
        for (j = LocalRange[1]; j < LocalRange[4]; ++j) {
          Points[0][iReceiver] = GlobalRange[3]-1;
          Points[1][iReceiver] = j;
          Sources[0][iReceiver] = 1;
          Sources[1][iReceiver] = j;
          ++iReceiver;
        }
      }

      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityNPoints(ConnectivityN, iDim, &Points[iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityNSources(ConnectivityN, iDim, &Sources[iDim]);
      }

      ovkRestoreConnectivityN(ConnectivityComponent, 2, 1, &ConnectivityN);

    }

    if (RightIsLocal) {

      const ovk_grid *Grid;
      ovkGetGrid(Domain, 2, &Grid);

      int GlobalRange[6];
      ovkGetGridGlobalRange(Grid, GlobalRange, GlobalRange+3);

      int LocalRange[6];
      ovkGetGridLocalRange(Grid, LocalRange, LocalRange+3);

      bool HasInterface = LocalRange[0] == GlobalRange[0];

      ovk_connectivity_m *ConnectivityM;
      ovkEditConnectivityM(ConnectivityComponent, 2, 1, &ConnectivityM);

      long long NumDonors = HasInterface ? (LocalRange[4]-LocalRange[1]) : 0;
      ovkResizeConnectivityM(ConnectivityM, NumDonors, 1);

      int *Extents[2][2];
      double *Coords[2];
      double *InterpCoefs[2][1];
      int *Destinations[2];
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityMExtents(ConnectivityM, iDim, &Extents[0][iDim], &Extents[1][iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityMCoords(ConnectivityM, iDim, &Coords[iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        for (iCoef = 0; iCoef < 1; ++iCoef) {
          ovkEditConnectivityMInterpCoefs(ConnectivityM, iDim, iCoef, &InterpCoefs[iDim][iCoef]);
        }
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityMDestinations(ConnectivityM, iDim, &Destinations[iDim]);
      }

      if (HasInterface) {
        iDonor = 0;
        for (j = LocalRange[1]; j < LocalRange[4]; ++j) {
          Extents[0][0][iDonor] = 1;
          Extents[0][1][iDonor] = j;
          Extents[1][0][iDonor] = Extents[0][0][iDonor]+1;
          Extents[1][1][iDonor] = Extents[0][1][iDonor]+1;
          Coords[0][iDonor] = 0.;
          Coords[1][iDonor] = 0.;
          InterpCoefs[0][0][iDonor] = 1.;
          InterpCoefs[1][0][iDonor] = 1.;
          Destinations[0][iDonor] = LeftSize[0]-1;
          Destinations[1][iDonor] = j;
          ++iDonor;
        }
      }

      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityMExtents(ConnectivityM, iDim, &Extents[0][iDim], &Extents[1][iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityMCoords(ConnectivityM, iDim, &Coords[iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        for (iCoef = 0; iCoef < 1; ++iCoef) {
          ovkRestoreConnectivityMInterpCoefs(ConnectivityM, iDim, iCoef,
            &InterpCoefs[iDim][iCoef]);
        }
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityMDestinations(ConnectivityM, iDim, &Destinations[iDim]);
      }

      ovkRestoreConnectivityM(ConnectivityComponent, 2, 1, &ConnectivityM);

      ovk_connectivity_n *ConnectivityN;
      ovkEditConnectivityN(ConnectivityComponent, 1, 2, &ConnectivityN);

      long long NumReceivers = HasInterface ? (LocalRange[4]-LocalRange[1]) : 0;
      ovkResizeConnectivityN(ConnectivityN, NumReceivers);

      int *Points[2];
      int *Sources[2];
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityNPoints(ConnectivityN, iDim, &Points[iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkEditConnectivityNSources(ConnectivityN, iDim, &Sources[iDim]);
      }

      if (HasInterface) {
        iReceiver = 0;
        for (j = LocalRange[1]; j < LocalRange[4]; ++j) {
          Points[0][iReceiver] = 0;
          Points[1][iReceiver] = j;
          Sources[0][iReceiver] = LeftSize[0]-2;
          Sources[1][iReceiver] = j;
          ++iReceiver;
        }
      }

      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityNPoints(ConnectivityN, iDim, &Points[iDim]);
      }
      for (iDim = 0; iDim < 2; ++iDim) {
        ovkRestoreConnectivityNSources(ConnectivityN, iDim, &Sources[iDim]);
      }

      ovkRestoreConnectivityN(ConnectivityComponent, 1, 2, &ConnectivityN);

    }

    ovkRestoreComponent(Domain, CONNECTIVITY_ID, OVK_COMPONENT_TYPE_CONNECTIVITY,
      &ConnectivityComponent);

  }

#ifdef OVK_HAVE_XDMF
  xdmf_grid_meta *XDMFGrids[2];
  CreateXDMFGridMeta(&XDMFGrids[0], "Left", LeftSize);
  CreateXDMFGridMeta(&XDMFGrids[1], "Right", RightSize);

  xdmf_attribute_meta *XDMFAttributes[3];
  CreateXDMFAttributeMeta(&XDMFAttributes[0], "State", XDMF_ATTRIBUTE_TYPE_INT);
  CreateXDMFAttributeMeta(&XDMFAttributes[1], "BeforeExchange", XDMF_ATTRIBUTE_TYPE_DOUBLE);
  CreateXDMFAttributeMeta(&XDMFAttributes[2], "AfterExchange", XDMF_ATTRIBUTE_TYPE_DOUBLE);

  xdmf *XDMF;
  xdmf_error XDMFError;
  CreateXDMF(&XDMF, "InterfaceC.xmf", 2, MPI_COMM_WORLD, 2, XDMFGrids, 3, XDMFAttributes,
    &XDMFError);
  CloseXDMF(&XDMF);

  if (LeftIsLocal) {
    const grid_data *Data = &LeftData;
    OpenXDMF(&XDMF, "InterfaceC.xmf", Data->Comm, &XDMFError);
    double *Coords[2];
    Coords[0] = malloc(Data->NumLocalPoints*sizeof(double));
    Coords[1] = malloc(Data->NumLocalPoints*sizeof(double));
    for (j = Data->LocalRange[1]; j < Data->LocalRange[4]; ++j) {
      for (i = Data->LocalRange[0]; i < Data->LocalRange[3]; ++i) {
        long long l = (i-Data->LocalRange[0]) + (Data->LocalRange[3]-Data->LocalRange[0])*
          (j-Data->LocalRange[1]);
        double U = (double)i/(double)(Size[0]-1);
        double V = (double)j/(double)(LeftSize[1]-1);
        Coords[0][l] = 2.*(U-0.5);
        Coords[1][l] = 2.*(V-0.5);
      }
    }
    int *OutputState = CreateOutputState(Domain, CONNECTIVITY_ID, 1);
    for (int iDim = 0; iDim < 2; ++iDim) {
      WriteXDMFGeometry(XDMF, "Left", iDim, Coords[iDim], Data->LocalRange, Data->LocalRange+3);
    }
    WriteXDMFAttribute(XDMF, "Left", "State", OutputState, Data->LocalRange, Data->LocalRange+3);
    free(OutputState);
    free(Coords[0]);
    free(Coords[1]);
    CloseXDMF(&XDMF);
  }

  if (RightIsLocal) {
    const grid_data *Data = &RightData;
    OpenXDMF(&XDMF, "InterfaceC.xmf", Data->Comm, &XDMFError);
    double *Coords[2];
    Coords[0] = malloc(Data->NumLocalPoints*sizeof(double));
    Coords[1] = malloc(Data->NumLocalPoints*sizeof(double));
    for (j = Data->LocalRange[1]; j < Data->LocalRange[4]; ++j) {
      for (i = Data->LocalRange[0]; i < Data->LocalRange[3]; ++i) {
        long long l = (i-Data->LocalRange[0]) + (Data->LocalRange[3]-Data->LocalRange[0])*
          (j-Data->LocalRange[1]);
        double U = (double)(i+LeftSize[0]-2)/(double)(Size[0]-1);
        double V = (double)j/(double)(RightSize[1]-1);
        Coords[0][l] = 2.*(U-0.5);
        Coords[1][l] = 2.*(V-0.5);
      }
    }
    int *OutputState = CreateOutputState(Domain, CONNECTIVITY_ID, 2);
    for (int iDim = 0; iDim < 2; ++iDim) {
      WriteXDMFGeometry(XDMF, "Right", iDim, Coords[iDim], Data->LocalRange, Data->LocalRange+3);
    }
    WriteXDMFAttribute(XDMF, "Right", "State", OutputState, Data->LocalRange, Data->LocalRange+3);
    free(OutputState);
    free(Coords[0]);
    free(Coords[1]);
    CloseXDMF(&XDMF);
  }
#endif

  ovk_exchanger *Exchanger;
  ovkCreateExchanger(&Exchanger, SharedContext, NULL);

  ovk_exchanger_bindings *ExchangerBindings;
  ovkCreateExchangerBindings(&ExchangerBindings);
  ovkSetExchangerBindingsConnectivityComponentID(ExchangerBindings, CONNECTIVITY_ID);
  ovkBindExchanger(Exchanger, Domain, &ExchangerBindings);

  const ovk_connectivity_component *ConnectivityComponent;
  ovkGetComponent(Domain, CONNECTIVITY_ID, OVK_COMPONENT_TYPE_CONNECTIVITY, &ConnectivityComponent);

  double *LeftDonorValues, *LeftReceiverValues;
  if (LeftIsLocal) {
    grid_data *Data = &LeftData;
    const ovk_connectivity_m *ConnectivityM;
    ovkGetConnectivityM(ConnectivityComponent, 1, 2, &ConnectivityM);
    ovkCreateExchangerCollect(Exchanger, 1, 2, 1, OVK_COLLECT_INTERPOLATE,
      OVK_DOUBLE, 1, Data->ExtendedRange, Data->ExtendedRange+3, OVK_ROW_MAJOR);
    ovkCreateExchangerSend(Exchanger, 1, 2, 1, OVK_DOUBLE, 1, 1);
    long long NumDonors = ovkGetConnectivityMSize(ConnectivityM);
    LeftDonorValues = malloc(NumDonors*sizeof(double));
    const ovk_connectivity_n *ConnectivityN;
    ovkGetConnectivityN(ConnectivityComponent, 2, 1, &ConnectivityN);
    ovkCreateExchangerReceive(Exchanger, 2, 1, 1, OVK_DOUBLE, 1, 1);
    ovkCreateExchangerDisperse(Exchanger, 2, 1, 1, OVK_DISPERSE_OVERWRITE, OVK_DOUBLE, 1,
      Data->ExtendedRange, Data->ExtendedRange+3, OVK_ROW_MAJOR);
    long long NumReceivers = ovkGetConnectivityNSize(ConnectivityN);
    LeftReceiverValues = malloc(NumReceivers*sizeof(double));
  }

  double *RightDonorValues, *RightReceiverValues;
  if (RightIsLocal) {
    grid_data *Data = &RightData;
    const ovk_connectivity_m *ConnectivityM;
    ovkGetConnectivityM(ConnectivityComponent, 2, 1, &ConnectivityM);
    ovkCreateExchangerCollect(Exchanger, 2, 1, 1, OVK_COLLECT_INTERPOLATE,
      OVK_DOUBLE, 1, Data->ExtendedRange, Data->ExtendedRange+3, OVK_ROW_MAJOR);
    ovkCreateExchangerSend(Exchanger, 2, 1, 1, OVK_DOUBLE, 1, 1);
    long long NumDonors = ovkGetConnectivityMSize(ConnectivityM);
    RightDonorValues = malloc(NumDonors*sizeof(double));
    const ovk_connectivity_n *ConnectivityN;
    ovkGetConnectivityN(ConnectivityComponent, 1, 2, &ConnectivityN);
    ovkCreateExchangerReceive(Exchanger, 1, 2, 1, OVK_DOUBLE, 1, 1);
    ovkCreateExchangerDisperse(Exchanger, 1, 2, 1, OVK_DISPERSE_OVERWRITE, OVK_DOUBLE, 1,
      Data->ExtendedRange, Data->ExtendedRange+3, OVK_ROW_MAJOR);
    long long NumReceivers = ovkGetConnectivityNSize(ConnectivityN);
    RightReceiverValues = malloc(NumReceivers*sizeof(double));
  }

  double *LeftFieldValues;
  if (LeftIsLocal) {
    grid_data *Data = &LeftData;
    LeftFieldValues = malloc(Data->NumExtendedPoints*sizeof(double));
    for (j = Data->ExtendedRange[1]; j < Data->ExtendedRange[4]; ++j) {
      for (i = Data->ExtendedRange[0]; i < Data->ExtendedRange[3]; ++i) {
        l = (Data->ExtendedRange[4]-Data->ExtendedRange[1])*(i-Data->ExtendedRange[0]) +
          (j-Data->ExtendedRange[1]);
        double U = (double)i/(double)(LeftSize[0]-1);
        double V = (double)j/(double)(LeftSize[1]-1);
        LeftFieldValues[l] = U*V;
      }
    }
  }

  double *RightFieldValues;
  if (RightIsLocal) {
    grid_data *Data = &RightData;
    RightFieldValues = malloc(Data->NumExtendedPoints*sizeof(double));
    for (j = Data->ExtendedRange[1]; j < Data->ExtendedRange[4]; ++j) {
      for (i = Data->ExtendedRange[0]; i < Data->ExtendedRange[3]; ++i) {
        l = (Data->ExtendedRange[4]-Data->ExtendedRange[1])*(i-Data->ExtendedRange[0]) +
          (j-Data->ExtendedRange[1]);
        double U = (double)i/(double)(RightSize[0]-1);
        double V = (double)j/(double)(RightSize[1]-1);
        RightFieldValues[l] = (1.-U)*(1.-V);
      }
    }
  }

#ifdef OVK_HAVE_XDMF
  if (LeftIsLocal) {
    const grid_data *Data = &LeftData;
    OpenXDMF(&XDMF, "InterfaceC.xmf", Data->Comm, &XDMFError);
    double *ValuesTransposed = malloc(Data->NumLocalPoints*sizeof(double));
    for (j = Data->LocalRange[1]; j < Data->LocalRange[4]; ++j) {
      for (i = Data->LocalRange[0]; i < Data->LocalRange[3]; ++i) {
        long long l_s = (Data->ExtendedRange[4]-Data->ExtendedRange[1])*(i-Data->ExtendedRange[0]) +
          (j-Data->ExtendedRange[1]);
        long long l_d = (i-Data->LocalRange[0]) + (Data->LocalRange[3]-Data->LocalRange[0])*
          (j-Data->LocalRange[1]);
        ValuesTransposed[l_d] = LeftFieldValues[l_s];
      }
    }
    WriteXDMFAttribute(XDMF, "Left", "BeforeExchange", ValuesTransposed, Data->LocalRange,
      Data->LocalRange+3);
    free(ValuesTransposed);
    CloseXDMF(&XDMF);
  }

  if (RightIsLocal) {
    const grid_data *Data = &RightData;
    OpenXDMF(&XDMF, "InterfaceC.xmf", Data->Comm, &XDMFError);
    double *ValuesTransposed = malloc(Data->NumLocalPoints*sizeof(double));
    for (j = Data->LocalRange[1]; j < Data->LocalRange[4]; ++j) {
      for (i = Data->LocalRange[0]; i < Data->LocalRange[3]; ++i) {
        long long l_s = (Data->ExtendedRange[4]-Data->ExtendedRange[1])*(i-Data->ExtendedRange[0]) +
          (j-Data->ExtendedRange[1]);
        long long l_d = (i-Data->LocalRange[0]) + (Data->LocalRange[3]-Data->LocalRange[0])*
          (j-Data->LocalRange[1]);
        ValuesTransposed[l_d] = RightFieldValues[l_s];
      }
    }
    WriteXDMFAttribute(XDMF, "Right", "BeforeExchange", ValuesTransposed, Data->LocalRange,
      Data->LocalRange+3);
    free(ValuesTransposed);
    CloseXDMF(&XDMF);
  }
#endif

  ovk_request **Requests = malloc(4*sizeof(ovk_request *));
  int NumRequests = 0;

  if (LeftIsLocal) {
    ovkExchangerReceive(Exchanger, 2, 1, 1, &LeftReceiverValues, &Requests[NumRequests]);
    ++NumRequests;
  }

  if (RightIsLocal) {
    ovkExchangerReceive(Exchanger, 1, 2, 1, &RightReceiverValues, &Requests[NumRequests]);
    ++NumRequests;
  }

  if (LeftIsLocal) {
    ovkExchangerCollect(Exchanger, 1, 2, 1, &LeftFieldValues, &LeftDonorValues);
    ovkExchangerSend(Exchanger, 1, 2, 1, &LeftDonorValues, &Requests[NumRequests]);
    ++NumRequests;
  }

  if (RightIsLocal) {
    ovkExchangerCollect(Exchanger, 2, 1, 1, &RightFieldValues, &RightDonorValues);
    ovkExchangerSend(Exchanger, 2, 1, 1, &RightDonorValues, &Requests[NumRequests]);
    ++NumRequests;
  }

  ovkWaitAll(NumRequests, Requests);

  free(Requests);

  if (LeftIsLocal) {
    ovkExchangerDisperse(Exchanger, 2, 1, 1, &LeftReceiverValues, &LeftFieldValues);
  }

  if (RightIsLocal) {
    ovkExchangerDisperse(Exchanger, 1, 2, 1, &RightReceiverValues, &RightFieldValues);
  }

#ifdef OVK_HAVE_XDMF
  if (LeftIsLocal) {
    const grid_data *Data = &LeftData;
    OpenXDMF(&XDMF, "InterfaceC.xmf", Data->Comm, &XDMFError);
    double *ValuesTransposed = malloc(Data->NumLocalPoints*sizeof(double));
    for (j = Data->LocalRange[1]; j < Data->LocalRange[4]; ++j) {
      for (i = Data->LocalRange[0]; i < Data->LocalRange[3]; ++i) {
        long long l_s = (Data->ExtendedRange[4]-Data->ExtendedRange[1])*(i-Data->ExtendedRange[0]) +
          (j-Data->ExtendedRange[1]);
        long long l_d = (i-Data->LocalRange[0]) + (Data->LocalRange[3]-Data->LocalRange[0])*
          (j-Data->LocalRange[1]);
        ValuesTransposed[l_d] = LeftFieldValues[l_s];
      }
    }
    WriteXDMFAttribute(XDMF, "Left", "AfterExchange", ValuesTransposed, Data->LocalRange,
      Data->LocalRange+3);
    free(ValuesTransposed);
    CloseXDMF(&XDMF);
  }

  if (RightIsLocal) {
    const grid_data *Data = &RightData;
    OpenXDMF(&XDMF, "InterfaceC.xmf", Data->Comm, &XDMFError);
    double *ValuesTransposed = malloc(Data->NumLocalPoints*sizeof(double));
    for (j = Data->LocalRange[1]; j < Data->LocalRange[4]; ++j) {
      for (i = Data->LocalRange[0]; i < Data->LocalRange[3]; ++i) {
        long long l_s = (Data->ExtendedRange[4]-Data->ExtendedRange[1])*(i-Data->ExtendedRange[0]) +
          (j-Data->ExtendedRange[1]);
        long long l_d = (i-Data->LocalRange[0]) + (Data->LocalRange[3]-Data->LocalRange[0])*
          (j-Data->LocalRange[1]);
        ValuesTransposed[l_d] = RightFieldValues[l_s];
      }
    }
    WriteXDMFAttribute(XDMF, "Right", "AfterExchange", ValuesTransposed, Data->LocalRange,
      Data->LocalRange+3);
    free(ValuesTransposed);
    CloseXDMF(&XDMF);
  }
#endif

  if (LeftIsLocal) {
    free(LeftDonorValues);
    free(LeftReceiverValues);
    free(LeftFieldValues);
  }

  if (RightIsLocal) {
    free(RightDonorValues);
    free(RightReceiverValues);
    free(RightFieldValues);
  }

  ovkDestroyExchanger(&Exchanger);
  ovkDestroyDomain(&Domain);

  ovkResetSharedContext(&SharedContext);

  DestroyGridData(&LeftData);
  DestroyGridData(&RightData);

  return 0;

}
