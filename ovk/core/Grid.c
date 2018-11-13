// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Grid.h"

#include "ovk/core/Cart.h"
#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Logger.h"
#include "ovk/core/MPIUtils.h"
#include "ovk/core/OrderedMap.h"
#include "ovk/core/PartitionHash.h"
#include "ovk/core/Range.h"
#include "ovk/core/TextUtils.h"

static void CreateNeighbors(ovk_grid *Grid);
static void DestroyNeighbors(ovk_grid *Grid);
static void PrintGridSummary(const ovk_grid *Grid);
static void PrintGridDecomposition(const ovk_grid *Grid);
static void DefaultProperties(ovk_grid_properties *Properties, int NumDims);

void PRIVATE(CreateGrid)(ovk_grid **Grid_, int ID, const ovk_grid_params *Params, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  int iDim;

  int NumDims = Params->num_dims;

  MPI_Comm Comm;
  MPI_Comm_dup(Params->comm, &Comm);

  int CommSize, CommRank;
  MPI_Comm_size(Comm, &CommSize);
  MPI_Comm_rank(Comm, &CommRank);

  MPI_Barrier(Comm);

  *Grid_ = malloc(sizeof(ovk_grid));
  ovk_grid *Grid = *Grid_;

  DefaultProperties(&Grid->properties, NumDims);

  Grid->properties.id = ID;

  if (strlen(Params->name) > 0) {
    strncpy(Grid->properties.name, Params->name, OVK_NAME_LENGTH);
  } else {
    sprintf(Grid->properties.name, "Grid%i", Grid->properties.id);
  }

  Grid->properties.comm = Comm;
  Grid->properties.comm_size = CommSize;
  Grid->properties.comm_rank = CommRank;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Grid->properties.periodic[iDim] = Params->periodic[iDim];
    Grid->properties.periodic_length[iDim] = Params->periodic_length[iDim];
  }

  Grid->properties.periodic_storage = Params->periodic_storage;
  Grid->properties.geometry_type = Params->geometry_type;

  ovkDefaultRange(&Grid->properties.global_range, NumDims);
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Grid->properties.global_range.b[iDim] = 0;
    Grid->properties.global_range.e[iDim] = Params->size[iDim];
  }

  Grid->properties.local_range = Params->local_range;

  Grid->logger = Logger;
  Grid->error_handler = ErrorHandler;

  ovkDefaultCart(&Grid->cart, NumDims);
  for (iDim = 0; iDim < NumDims; ++iDim) {
    if (Grid->properties.periodic[iDim] && Grid->properties.periodic_storage ==
      OVK_PERIODIC_STORAGE_DUPLICATED) {
      Grid->cart.size[iDim] = Grid->properties.global_range.e[iDim]-1;
    } else {
      Grid->cart.size[iDim] = Grid->properties.global_range.e[iDim];
    }
    Grid->cart.periodic[iDim] = Grid->properties.periodic[iDim];
  }

  CreatePartitionHash(&Grid->partition_hash, Grid->properties.num_dims, Grid->properties.comm,
    &Grid->properties.global_range, &Grid->properties.local_range);

  CreateNeighbors(Grid);

  if (Grid->properties.comm_rank == 0) {
    PrintGridSummary(Grid);
  }

  if (OVK_DEBUG) {
    PrintGridDecomposition(Grid);
  }

  MPI_Barrier(Grid->properties.comm);

}

void PRIVATE(DestroyGrid)(ovk_grid **Grid_) {

  ovk_grid *Grid = *Grid_;

  MPI_Barrier(Grid->properties.comm);

  DestroyPartitionHash(&Grid->partition_hash);

  DestroyNeighbors(Grid);

  t_logger *Logger = Grid->logger;
  MPI_Comm Comm = Grid->properties.comm;
  bool IsRoot = Grid->properties.comm_rank == 0;
  char Name[OVK_NAME_LENGTH];
  strncpy(Name, Grid->properties.name, OVK_NAME_LENGTH);

  free_null(Grid_);

  MPI_Barrier(Comm);

  MPI_Comm_free(&Comm);

  LogStatus(Logger, IsRoot, 0, "Destroyed grid %s.", Name);

}

void PRIVATE(CreateGridInfo)(ovk_grid_info **Info_, const ovk_grid *Grid, MPI_Comm Comm,
  int CommRank) {

  int iDim;

  *Info_ = malloc(sizeof(ovk_grid_info));
  ovk_grid_info *Info = *Info_;

  bool IsLocal = Grid != NULL;
  bool IsRoot = false;
  if (IsLocal) {
    IsRoot = Grid->properties.comm_rank == 0;
  }

  int RootRank;
  if (IsRoot) RootRank = CommRank;
  BroadcastAnySource(&RootRank, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) {
    Info->id = Grid->properties.id;
    strcpy(Info->name, Grid->properties.name);
    Info->num_dims = Grid->properties.num_dims;
  }
  MPI_Bcast(&Info->id, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info->name, OVK_NAME_LENGTH, MPI_CHAR, RootRank, Comm);
  MPI_Bcast(&Info->num_dims, 1, MPI_INT, RootRank, Comm);

  Info->root_rank = RootRank;

  if (IsRoot) {
    Info->global_range = Grid->properties.global_range;
  }
  MPI_Bcast(&Info->global_range.nd, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info->global_range.b, MAX_DIMS, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info->global_range.e, MAX_DIMS, MPI_INT, RootRank, Comm);

  int PeriodicInt[MAX_DIMS];
  if (IsRoot) {
    Info->cart = Grid->cart;
    PeriodicInt[0] = (int)Info->cart.periodic[0];
    PeriodicInt[1] = (int)Info->cart.periodic[1];
    PeriodicInt[2] = (int)Info->cart.periodic[2];
  }
  MPI_Bcast(&Info->cart.nd, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info->cart.size, MAX_DIMS, MPI_INT, RootRank, Comm);
  MPI_Bcast(&PeriodicInt, MAX_DIMS, MPI_INT, RootRank, Comm);
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Info->cart.periodic[iDim] = (bool)PeriodicInt[iDim];
  }

  if (IsRoot) {
    Info->periodic_length[0] = Grid->properties.periodic_length[0];
    Info->periodic_length[1] = Grid->properties.periodic_length[1];
    Info->periodic_length[2] = Grid->properties.periodic_length[2];
    Info->geometry_type = Grid->properties.geometry_type;
  }
  MPI_Bcast(&Info->periodic_length, MAX_DIMS, MPI_DOUBLE, RootRank, Comm);
  MPI_Bcast(&Info->geometry_type, 1, MPI_INT, RootRank, Comm);

}

void PRIVATE(DestroyGridInfo)(ovk_grid_info **Info) {

  free_null(Info);

}

void ovkGetGridProperties(const ovk_grid *Grid, const ovk_grid_properties **Properties) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");

  *Properties = &Grid->properties;

}

void ovkGetGridCart(const ovk_grid *Grid, ovk_cart *Cart) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(Cart, "Invalid cart pointer.");

  *Cart = Grid->cart;

}

static void CreateNeighbors(ovk_grid *Grid) {

  int iDim, jDim;
  int iFace;
  int i, j, k;
  size_t iPoint;
  int iNeighbor;
  t_ordered_map_entry *Entry;

  int NumDims = Grid->properties.num_dims;

  ovk_range GlobalRange = Grid->properties.global_range;
  ovk_range LocalRange = Grid->properties.local_range;
  int Periodic[MAX_DIMS] = {
    Grid->properties.periodic[0],
    Grid->properties.periodic[1],
    Grid->properties.periodic[2]
  };

  bool HasNeighborsBefore[MAX_DIMS];
  bool HasNeighborsAfter[MAX_DIMS];
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    bool HasPartitionBefore = LocalRange.b[iDim] > GlobalRange.b[iDim];
    bool HasPartitionAfter = LocalRange.e[iDim] < GlobalRange.e[iDim];
    if (Periodic[iDim]) {
      HasNeighborsBefore[iDim] = HasPartitionBefore || HasPartitionAfter;
      HasNeighborsAfter[iDim] = HasPartitionBefore || HasPartitionAfter;
    } else {
      HasNeighborsBefore[iDim] = HasPartitionBefore;
      HasNeighborsAfter[iDim] = HasPartitionAfter;
    }
  }

  int NumNeighborFaces = 0;
  ovk_range NeighborFaces[2*MAX_DIMS];
  for (iDim = 0; iDim < NumDims; ++iDim) {
    if (HasNeighborsBefore[iDim]) {
      ovk_range *Face = NeighborFaces+NumNeighborFaces;
      *Face = LocalRange;
      Face->b[iDim] -= 1;
      Face->e[iDim] = Face->b[iDim]+1;
      for (jDim = iDim+1; jDim < NumDims; ++jDim) {
        Face->b[jDim] -= (int)HasNeighborsBefore[jDim];
        Face->e[jDim] += (int)HasNeighborsAfter[jDim];
      }
      ++NumNeighborFaces;
    }
    if (HasNeighborsAfter[iDim]) {
      ovk_range *Face = NeighborFaces+NumNeighborFaces;
      *Face = LocalRange;
      Face->e[iDim] += 1;
      Face->b[iDim] = Face->e[iDim]-1;
      for (jDim = iDim+1; jDim < NumDims; ++jDim) {
        Face->b[jDim] -= (int)HasNeighborsBefore[jDim];
        Face->e[jDim] += (int)HasNeighborsAfter[jDim];
      }
      ++NumNeighborFaces;
    }
  }

  size_t NumNeighborPoints = 0;
  for (iFace = 0; iFace < NumNeighborFaces; ++iFace) {
    size_t NumPoints;
    ovkRangeCount(NeighborFaces+iFace, &NumPoints);
    NumNeighborPoints += NumPoints;
  }

  int *NeighborPoints[MAX_DIMS];
  NeighborPoints[0] = malloc(MAX_DIMS*NumNeighborPoints*sizeof(int));
  NeighborPoints[1] = NeighborPoints[0] + NumNeighborPoints;
  NeighborPoints[2] = NeighborPoints[1] + NumNeighborPoints;

  iPoint = 0;
  for (iFace = 0; iFace < NumNeighborFaces; ++iFace) {
    ovk_range *Face = NeighborFaces+iFace;
    for (k = Face->b[2]; k < Face->e[2]; ++k) {
      for (j = Face->b[1]; j < Face->e[1]; ++j) {
        for (i = Face->b[0]; i < Face->e[0]; ++i) {
          int Point[MAX_DIMS] = {i, j, k};
          ovkCartPeriodicAdjust(&Grid->cart, Point, Point);
          for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
            NeighborPoints[iDim][iPoint] = Point[iDim];
          }
          ++iPoint;
        }
      }
    }
  }

  int *NeighborPointBinIndices = malloc(NumNeighborPoints*sizeof(int));

  MapToPartitionBins(Grid->partition_hash, NumNeighborPoints, (const int **)NeighborPoints,
    NeighborPointBinIndices);

  t_ordered_map *Bins;
  OMCreate(&Bins);

  for (iPoint = 0; iPoint < NumNeighborPoints; ++iPoint) {
    int BinIndex = NeighborPointBinIndices[iPoint];
    if (!OMExists(Bins, BinIndex)) {
      OMInsert(Bins, BinIndex, NULL);
    }
  }

  RetrievePartitionBins(Grid->partition_hash, Bins);

  int *NeighborPointRanks = malloc(NumNeighborPoints*sizeof(int));

  FindPartitions(Grid->partition_hash, Bins, NumNeighborPoints, (const int **)NeighborPoints,
    NeighborPointBinIndices, NeighborPointRanks);

  ClearPartitionBins(Bins);
  OMDestroy(&Bins);

  t_ordered_map *NeighborRanks;
  OMCreate(&NeighborRanks);

  for (iPoint = 0; iPoint < NumNeighborPoints; ++iPoint) {
    int Rank = NeighborPointRanks[iPoint];
    if (!OMExists(NeighborRanks, Rank)) {
      OMInsert(NeighborRanks, Rank, NULL);
    }
  }

  free(NeighborPoints[0]);
  free(NeighborPointBinIndices);
  free(NeighborPointRanks);

  int NumNeighbors = OMSize(NeighborRanks);

  Grid->num_neighbors = NumNeighbors;
  Grid->neighbors = malloc(NumNeighbors*sizeof(t_grid_neighbor));

  int *NeighborRangesFlat = malloc(2*MAX_DIMS*NumNeighbors*sizeof(int));

  int LocalRangeFlat[2*MAX_DIMS];
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    LocalRangeFlat[iDim] = LocalRange.b[iDim];
    LocalRangeFlat[MAX_DIMS+iDim] = LocalRange.e[iDim];
  }

  MPI_Request *Requests = malloc(2*NumNeighbors*sizeof(MPI_Request));
  Entry = OMBegin(NeighborRanks);
  iNeighbor = 0;
  while (Entry != OMEnd(NeighborRanks)) {
    int Rank = OMKey(Entry);
    MPI_Irecv(NeighborRangesFlat+2*MAX_DIMS*iNeighbor, 2*MAX_DIMS, MPI_INT, Rank, 0,
      Grid->properties.comm, Requests+2*iNeighbor);
    MPI_Isend(LocalRangeFlat, 2*MAX_DIMS, MPI_INT, Rank, 0, Grid->properties.comm,
      Requests+2*iNeighbor+1);
    Entry = OMNext(Entry);
    ++iNeighbor;
  }
  MPI_Waitall(2*NumNeighbors, Requests, MPI_STATUSES_IGNORE);

  Entry = OMBegin(NeighborRanks);
  iNeighbor = 0;
  while (Entry != OMEnd(NeighborRanks)) {
    int Rank = OMKey(Entry);
    Grid->neighbors[iNeighbor].comm_rank = Rank;
    ovkDefaultRange(&Grid->neighbors[iNeighbor].local_range, Grid->properties.num_dims);
    int *NeighborRangeFlat = NeighborRangesFlat+2*MAX_DIMS*iNeighbor;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Grid->neighbors[iNeighbor].local_range.b[iDim] = NeighborRangeFlat[iDim];
      Grid->neighbors[iNeighbor].local_range.e[iDim] = NeighborRangeFlat[MAX_DIMS+iDim];
    }
    Entry = OMNext(Entry);
    ++iNeighbor;
  }

  free(Requests);
  free(NeighborRangesFlat);

  OMDestroy(&NeighborRanks);

}

static void DestroyNeighbors(ovk_grid *Grid) {

  free(Grid->neighbors);

}

static void PrintGridSummary(const ovk_grid *Grid) {

  int iDim;
  int Offset;

  char GlobalSizeString[256];
  Offset = 0;
  for (iDim = 0; iDim < Grid->properties.num_dims; ++iDim) {
    char GlobalEndString[NUMBER_STRING_LENGTH];
    IntToString(Grid->properties.global_range.e[iDim], GlobalEndString);
    Offset += sprintf(GlobalSizeString+Offset, "%s", GlobalEndString);
    if (iDim != Grid->properties.num_dims-1) {
      Offset += sprintf(GlobalSizeString+Offset, " x ");
    }
  }

  size_t TotalPoints;
  ovkRangeCount(&Grid->properties.global_range, &TotalPoints);

  char TotalPointsString[NUMBER_STRING_LENGTH+7];
  PluralizeLabel(TotalPoints, "points", "point", TotalPointsString);

  char ProcessesString[NUMBER_STRING_LENGTH+10];
  PluralizeLabel(Grid->properties.comm_size, "processes", "process", ProcessesString);

  LogStatus(Grid->logger, true, 0, "Created grid %s (ID=%i): %s (%s) on %s.",
    Grid->properties.name, Grid->properties.id, GlobalSizeString, TotalPointsString,
    ProcessesString);

}

static void PrintGridDecomposition(const ovk_grid *Grid) {

  int iDim;
  int iNeighbor;
  int Offset;

  char DimNames[3] = {'i', 'j', 'k'};
  char LocalRangeString[256];
  Offset = 0;
  for (iDim = 0; iDim < Grid->properties.num_dims; ++iDim) {
    char LocalBeginString[NUMBER_STRING_LENGTH];
    char LocalEndString[NUMBER_STRING_LENGTH];
    IntToString(Grid->properties.local_range.b[iDim], LocalBeginString);
    IntToString(Grid->properties.local_range.e[iDim], LocalEndString);
    Offset += sprintf(LocalRangeString+Offset, "%c=%s:%s", DimNames[iDim], LocalBeginString,
      LocalEndString);
    if (iDim != Grid->properties.num_dims-1) {
      Offset += sprintf(LocalRangeString+Offset, ", ");
    }
  }

  char TotalLocalPointsString[NUMBER_STRING_LENGTH+7];
  size_t TotalLocalPoints;
  ovkRangeCount(&Grid->properties.local_range, &TotalLocalPoints);
  PluralizeLabel(TotalLocalPoints, "points", "point", TotalLocalPointsString);

  char NeighborRanksString[256];
  Offset = 0;
  for (iNeighbor = 0; iNeighbor < Grid->num_neighbors; ++iNeighbor) {
    Offset += sprintf(NeighborRanksString+Offset, "%i", Grid->neighbors[iNeighbor].comm_rank);
    if (iNeighbor != Grid->num_neighbors-1) {
      Offset += sprintf(NeighborRanksString+Offset, ", ");
    }
  }

  LogStatus(Grid->logger, Grid->properties.comm_rank == 0, 0, "Grid %s decomposition info:",
    Grid->properties.name);

  MPI_Barrier(Grid->properties.comm);

  int OtherRank;
  for (OtherRank = 0; OtherRank < Grid->properties.comm_size; ++OtherRank) {
    if (OtherRank == Grid->properties.comm_rank) {
      LogStatus(Grid->logger, true, 1, "Rank %i (global rank @rank@) contains %s (%s).",
        Grid->properties.comm_rank, LocalRangeString, TotalLocalPointsString);
      if (Grid->num_neighbors > 0) {
        LogStatus(Grid->logger, true, 1, "Rank %i has neighbors: %s", Grid->properties.comm_rank,
          NeighborRanksString);
      }
    }
    MPI_Barrier(Grid->properties.comm);
  }

}

void PRIVATE(GetGridNeighborCount)(const ovk_grid *Grid, int *NumNeighbors) {

  *NumNeighbors = Grid->num_neighbors;

}

void PRIVATE(GetGridNeighborRange)(const ovk_grid *Grid, int iNeighbor, ovk_range *NeighborRange) {

  *NeighborRange = Grid->neighbors[iNeighbor].local_range;

}

void PRIVATE(GetGridNeighborRank)(const ovk_grid *Grid, int iNeighbor, int *NeighborRank) {

  *NeighborRank = Grid->neighbors[iNeighbor].comm_rank;

}

void ovkCreateGridParams(ovk_grid_params **Params_, int NumDims) {

  OVK_DEBUG_ASSERT(Params_, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(NumDims == 2 || NumDims == 3, "Invalid dimension.");

  int iDim;

  *Params_ = malloc(sizeof(ovk_grid_params));
  ovk_grid_params *Params = *Params_;

  Params->num_dims = NumDims;
  memset(Params->name, 0, OVK_NAME_LENGTH);
  Params->comm = MPI_COMM_NULL;

  for (iDim = 0; iDim < NumDims; ++iDim) {
    Params->size[iDim] = 0;
  }
  for (iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
    Params->size[iDim] = 1;
  }

  Params->periodic[0] = false;
  Params->periodic[1] = false;
  Params->periodic[2] = false;
  Params->periodic_storage = OVK_PERIODIC_STORAGE_UNIQUE;
  Params->periodic_length[0] = 0.;
  Params->periodic_length[1] = 0.;
  Params->periodic_length[2] = 0.;
  Params->geometry_type = OVK_GEOMETRY_CURVILINEAR;

  ovkDefaultRange(&Params->local_range, NumDims);

}

void ovkDestroyGridParams(ovk_grid_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  free_null(Params);

}

void ovkGetGridParamName(const ovk_grid_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  strcpy(Name, Params->name);

}

void ovkSetGridParamName(ovk_grid_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  strncpy(Params->name, Name, OVK_NAME_LENGTH);

}

void ovkGetGridParamDimension(const ovk_grid_params *Params, int *NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Params->num_dims;

}

void ovkGetGridParamComm(const ovk_grid_params *Params, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Params->comm;

}

void ovkSetGridParamComm(ovk_grid_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Params->comm = Comm;

}

void ovkGetGridParamSize(const ovk_grid_params *Params, int *Size) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  int iDim;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Size[iDim] = Params->size[iDim];
  }

}

void ovkSetGridParamSize(ovk_grid_params *Params, const int *Size) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  int iDim;

  for (iDim = 0; iDim < Params->num_dims; ++iDim) {
    OVK_DEBUG_ASSERT(Size[iDim] > 0, "Size must be greater than 0 in each dimension.");
    Params->size[iDim] = Size[iDim];
  }

}

void ovkGetGridParamPeriodic(const ovk_grid_params *Params, bool *Periodic) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  int iDim;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Periodic[iDim] = Params->periodic[iDim];
  }

}

void ovkSetGridParamPeriodic(ovk_grid_params *Params, const bool *Periodic) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  int iDim;

  for (iDim = 0; iDim < Params->num_dims; ++iDim) {
    Params->periodic[iDim] = Periodic[iDim];
  }

}

void ovkGetGridParamPeriodicStorage(const ovk_grid_params *Params,
  ovk_periodic_storage *PeriodicStorage) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(PeriodicStorage, "Invalid periodic storage pointer.");

  *PeriodicStorage = Params->periodic_storage;

}

void ovkSetGridParamPeriodicStorage(ovk_grid_params *Params, ovk_periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(ValidPeriodicStorage(PeriodicStorage), "Invalid periodic storage.");

  OVK_DEBUG_ASSERT(PeriodicStorage == OVK_PERIODIC_STORAGE_DUPLICATED, "Duplicated periodic storage "
    "is not currently supported.");

  Params->periodic_storage = PeriodicStorage;

}

void ovkGetGridParamPeriodicLength(const ovk_grid_params *Params, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  int iDim;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = Params->periodic_length[iDim];
  }

}

void ovkSetGridParamPeriodicLength(ovk_grid_params *Params, const double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  int iDim;

  for (iDim = 0; iDim < Params->num_dims; ++iDim) {
    OVK_DEBUG_ASSERT(PeriodicLength[iDim] >= 0., "Periodic length must be nonnegative.");
    Params->periodic_length[iDim] = PeriodicLength[iDim];
  }

}

void ovkGetGridParamGeometryType(const ovk_grid_params *Params, ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  *GeometryType = Params->geometry_type;

}

void ovkSetGridParamGeometryType(ovk_grid_params *Params, ovk_geometry_type GeometryType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(ValidGeometryType(GeometryType), "Invalid geometry type.");

  Params->geometry_type = GeometryType;

}

void ovkGetGridParamLocalRange(const ovk_grid_params *Params, ovk_range *LocalRange) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LocalRange, "Invalid local range pointer.");

  *LocalRange = Params->local_range;

}

void ovkSetGridParamLocalRange(ovk_grid_params *Params, const ovk_range *LocalRange) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LocalRange, "Invalid local range pointer.");

  Params->local_range = *LocalRange;

}

static void DefaultProperties(ovk_grid_properties *Properties, int NumDims) {

  Properties->id = -1;

  memset(Properties->name, 0, OVK_NAME_LENGTH);

  Properties->num_dims = NumDims;
  Properties->comm = MPI_COMM_NULL;
  Properties->comm_size = 0;
  Properties->comm_rank = 0;
  Properties->periodic[0] = false;
  Properties->periodic[1] = false;
  Properties->periodic[2] = false;
  Properties->periodic_storage = OVK_PERIODIC_STORAGE_UNIQUE;
  Properties->periodic_length[0] = 0.;
  Properties->periodic_length[1] = 0.;
  Properties->periodic_length[2] = 0.;
  Properties->geometry_type = OVK_GEOMETRY_CURVILINEAR;

  ovkDefaultRange(&Properties->global_range, NumDims);
  ovkDefaultRange(&Properties->local_range, NumDims);

}

void ovkGetGridPropertyID(const ovk_grid_properties *Properties, int *ID) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(ID, "Invalid ID pointer.");

  *ID = Properties->id;

}

void ovkGetGridPropertyName(const ovk_grid_properties *Properties, char *Name) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  strcpy(Name, Properties->name);

}

void ovkGetGridPropertyDimension(const ovk_grid_properties *Properties, int *NumDims) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Properties->num_dims;

}

void ovkGetGridPropertyComm(const ovk_grid_properties *Properties, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Properties->comm;

}

void ovkGetGridPropertyCommSize(const ovk_grid_properties *Properties, int *CommSize) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  *CommSize = Properties->comm_size;

}

void ovkGetGridPropertyCommRank(const ovk_grid_properties *Properties, int *CommRank) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  *CommRank = Properties->comm_rank;

}

void ovkGetGridPropertySize(const ovk_grid_properties *Properties, int *Size) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  int iDim;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Size[iDim] = Properties->global_range.e[iDim];
  }

}

void ovkGetGridPropertyPeriodic(const ovk_grid_properties *Properties, bool *Periodic) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  int iDim;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Periodic[iDim] = Properties->periodic[iDim];
  }

}

void ovkGetGridPropertyPeriodicStorage(const ovk_grid_properties *Properties,
  ovk_periodic_storage *PeriodicStorage) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(PeriodicStorage, "Invalid periodic storage pointer.");

  *PeriodicStorage = Properties->periodic_storage;

}

void ovkGetGridPropertyPeriodicLength(const ovk_grid_properties *Properties,
  double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  int iDim;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = Properties->periodic_length[iDim];
  }

}

void ovkGetGridPropertyGeometryType(const ovk_grid_properties *Properties,
  ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  *GeometryType = Properties->geometry_type;

}

void ovkGetGridPropertyGlobalRange(const ovk_grid_properties *Properties, ovk_range *GlobalRange) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(GlobalRange, "Invalid global range pointer.");

  *GlobalRange = Properties->global_range;

}

void ovkGetGridPropertyLocalRange(const ovk_grid_properties *Properties, ovk_range *LocalRange) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(LocalRange, "Invalid local range pointer.");

  *LocalRange = Properties->local_range;

}

void ovkGetGridInfoID(const ovk_grid_info *Info, int *ID) {

  OVK_DEBUG_ASSERT(Info, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(ID, "Invalid ID pointer.");

  *ID = Info->id;

}

void ovkGetGridInfoName(const ovk_grid_info *Info, char *Name) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  strcpy(Name, Info->name);

}

void ovkGetGridInfoDimension(const ovk_grid_info *Info, int *NumDims) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Info->num_dims;

}

void ovkGetGridInfoRootRank(const ovk_grid_info *Info, int *RootRank) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(RootRank, "Invalid root rank pointer.");

  *RootRank = Info->root_rank;

}

void ovkGetGridInfoGlobalRange(const ovk_grid_info *Info, ovk_range *GlobalRange) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(GlobalRange, "Invalid global range pointer.");

  *GlobalRange = Info->global_range;

}

void ovkGetGridInfoCart(const ovk_grid_info *Info, ovk_cart *Cart) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Cart, "Invalid cart pointer.");

  *Cart = Info->cart;

}


void ovkGetGridInfoPeriodicLength(const ovk_grid_info *Info, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  int iDim;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = Info->periodic_length[iDim];
  }

}

void ovkGetGridInfoGeometryType(const ovk_grid_info *Info, ovk_geometry_type *GeometryType) {

  OVK_DEBUG_ASSERT(Info, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(GeometryType, "Invalid geometry type pointer.");

  *GeometryType = Info->geometry_type;

}
