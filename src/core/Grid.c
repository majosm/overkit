// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "Grid.h"

#include "Cart.h"
#include "Debug.h"
#include "ErrorHandler.h"
#include "Global.h"
#include "Logger.h"
#include "Range.h"
#include "TextUtils.h"

static void CreateNeighborInfo(ovk_grid *Grid);
static void DestroyNeighborInfo(ovk_grid *Grid);
static void PrintGridSummary(const ovk_grid *Grid);
static void PrintGridDecomposition(const ovk_grid *Grid);

static void DefaultGridProperties(ovk_grid_properties *Properties);

void CreateGrid(ovk_grid **Grid_, int ID, const ovk_grid_params *Params, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  int i;

  MPI_Comm Comm;
  MPI_Comm_dup(Params->comm, &Comm);

  MPI_Barrier(Comm);

  *Grid_ = malloc(sizeof(ovk_grid));
  ovk_grid *Grid = *Grid_;

  DefaultGridProperties(&Grid->properties);

  Grid->properties.id = ID;

  if (strlen(Params->name) > 0) {
    strncpy(Grid->properties.name, Params->name, OVK_NAME_LENGTH);
  } else {
    sprintf(Grid->properties.name, "Grid%i", Grid->properties.id);
  }

  Grid->properties.num_dims = Params->num_dims;

  Grid->properties.comm = Comm;
  MPI_Comm_size(Grid->properties.comm, &Grid->properties.comm_size);
  MPI_Comm_rank(Grid->properties.comm, &Grid->properties.comm_rank);

  for (i = 0; i < MAX_DIMS; ++i) {
    Grid->properties.size[i] = Params->size[i];
    Grid->properties.periodic[i] = Params->periodic[i];
    Grid->properties.periodic_length[i] = Params->periodic_length[i];
  }

  Grid->properties.periodic_storage = Params->periodic_storage;
  Grid->properties.geometry_type = Params->geometry_type;

  Grid->properties.local_range = Params->local_range;

  Grid->properties.num_neighbors = Params->num_neighbors;
  if (Params->num_neighbors > 0) {
    Grid->properties.neighbor_ranks = malloc(Params->num_neighbors*sizeof(int));
  }

  for (i = 0; i < Params->num_neighbors; ++i) {
    Grid->properties.neighbor_ranks[i] = Params->neighbor_ranks[i];
  }

  Grid->logger = Logger;
  Grid->error_handler = ErrorHandler;

  CreateNeighborInfo(Grid);

  ovkCartDefault(&Grid->cart, Grid->properties.num_dims);
  for (i = 0; i < Grid->properties.num_dims; ++i) {
    if (Grid->properties.periodic[i] && Grid->properties.periodic_storage == OVK_OVERLAP_PERIODIC) {
      Grid->cart.size[i] = Grid->properties.size[i]-1;
    } else {
      Grid->cart.size[i] = Grid->properties.size[i];
    }
    Grid->cart.periodic[i] = Grid->properties.periodic[i];
  }

  if (Grid->properties.comm_rank == 0) {
    PrintGridSummary(Grid);
  }

  if (OVK_DEBUG) {
    PrintGridDecomposition(Grid);
  }

  MPI_Barrier(Grid->properties.comm);

}

void DestroyGrid(ovk_grid **Grid_) {

  ovk_grid *Grid = *Grid_;

  MPI_Barrier(Grid->properties.comm);

  DestroyNeighborInfo(Grid);

  t_logger *Logger = Grid->logger;
  MPI_Comm Comm = Grid->properties.comm;
  bool IsGridRoot = Grid->properties.comm_rank == 0;
  char Name[OVK_NAME_LENGTH];
  strncpy(Name, Grid->properties.name, OVK_NAME_LENGTH);

  free(Grid->properties.neighbor_ranks);

  free(*Grid_);
  *Grid_ = NULL;

  MPI_Barrier(Comm);

  MPI_Comm_free(&Comm);

  LogStatus(Logger, IsGridRoot, 0, "Destroyed grid %s.", Name);

}

void ovkGetGridProperties(const ovk_grid *Grid, const ovk_grid_properties **Properties) {

  *Properties = &Grid->properties;

}

static void CreateNeighborInfo(ovk_grid *Grid) {

  int i, j;

  int NumNeighbors = Grid->properties.num_neighbors;
  int *NeighborRanks = Grid->properties.neighbor_ranks;

  Grid->neighbors = malloc(NumNeighbors*sizeof(t_grid_neighbor_info));

  int *NeighborIndexRanges = malloc(6*NumNeighbors*sizeof(int));

  int IndexRange[6];
  IndexRange[0] = Grid->properties.local_range.b[0];
  IndexRange[1] = Grid->properties.local_range.b[1];
  IndexRange[2] = Grid->properties.local_range.b[2];
  IndexRange[3] = Grid->properties.local_range.e[0];
  IndexRange[4] = Grid->properties.local_range.e[1];
  IndexRange[5] = Grid->properties.local_range.e[2];

  MPI_Request *Requests = malloc(2*NumNeighbors*sizeof(MPI_Request));
  for (i = 0; i < NumNeighbors; ++i) {
    MPI_Request *NeighborRequests = Requests + 2*i;
    MPI_Irecv(NeighborIndexRanges+6*i, 6, MPI_INT, NeighborRanks[i], 0, Grid->properties.comm,
      NeighborRequests);
    MPI_Isend(IndexRange, 6, MPI_INT, NeighborRanks[i], 0, Grid->properties.comm,
      NeighborRequests+1);
  }
  MPI_Waitall(2*NumNeighbors, Requests, MPI_STATUSES_IGNORE);

  for (i = 0; i < NumNeighbors; ++i) {
    Grid->neighbors[i].comm_rank = NeighborRanks[i];
    ovkRangeDefault(&Grid->neighbors[i].local_range, Grid->properties.num_dims);
    int *NeighborIndexRange = NeighborIndexRanges+6*i;
    for (j = 0; j < MAX_DIMS; ++j) {
      Grid->neighbors[i].local_range.b[j] = NeighborIndexRange[j];
      Grid->neighbors[i].local_range.e[j] = NeighborIndexRanges[MAX_DIMS+j];
    }
  }

  free(Requests);
  free(NeighborIndexRanges);

}

static void DestroyNeighborInfo(ovk_grid *Grid) {

  free(Grid->neighbors);

}

static void PrintGridSummary(const ovk_grid *Grid) {

  char IDString[NUMBER_STRING_LENGTH];
  IntToString(Grid->properties.id, IDString);

  size_t TotalPoints =
    (size_t)Grid->properties.size[0] *
    (size_t)Grid->properties.size[1] *
    (size_t)Grid->properties.size[2];

  char TotalPointsString[NUMBER_STRING_LENGTH+7];
  PluralizeLabel(TotalPoints, "points", "point", TotalPointsString);

  char ProcessesString[NUMBER_STRING_LENGTH+10];
  PluralizeLabel(Grid->properties.comm_size, "processes", "process", ProcessesString);

  char ISizeString[NUMBER_STRING_LENGTH];
  char JSizeString[NUMBER_STRING_LENGTH];
  char KSizeString[NUMBER_STRING_LENGTH];
  SizeToString(Grid->properties.size[0], ISizeString);
  SizeToString(Grid->properties.size[1], JSizeString);
  SizeToString(Grid->properties.size[2], KSizeString);

  switch (Grid->properties.num_dims) {
  case 2:
    LogStatus(Grid->logger, true, 0, "Created grid %s (ID=%s): %s x %s (%s) on %s.",
      Grid->properties.name, IDString, ISizeString, JSizeString, TotalPointsString,
      ProcessesString);
    break;
  case 3:
    LogStatus(Grid->logger, true, 0, "Created grid %s (ID=%s): %s x %s x %s (%s) on %s.",
      Grid->properties.name, IDString, ISizeString, JSizeString, KSizeString, TotalPointsString,
      ProcessesString);
    break;
  }

}

static void PrintGridDecomposition(const ovk_grid *Grid) {

  int i;

  char IDString[NUMBER_STRING_LENGTH];
  IntToString(Grid->properties.id, IDString);

  char RankString[NUMBER_STRING_LENGTH];
  IntToString(Grid->properties.comm_rank, RankString);

  char ILocalBeginString[NUMBER_STRING_LENGTH], ILocalEndString[NUMBER_STRING_LENGTH];
  char JLocalBeginString[NUMBER_STRING_LENGTH], JLocalEndString[NUMBER_STRING_LENGTH];
  char KLocalBeginString[NUMBER_STRING_LENGTH], KLocalEndString[NUMBER_STRING_LENGTH];
  IntToString(Grid->properties.local_range.b[0], ILocalBeginString);
  IntToString(Grid->properties.local_range.b[1], JLocalBeginString);
  IntToString(Grid->properties.local_range.b[2], KLocalBeginString);
  IntToString(Grid->properties.local_range.e[0], ILocalEndString);
  IntToString(Grid->properties.local_range.e[1], JLocalEndString);
  IntToString(Grid->properties.local_range.e[2], KLocalEndString);

  char TotalLocalPointsString[NUMBER_STRING_LENGTH+7];
  size_t TotalLocalPoints;
  ovkRangeCount(&Grid->properties.local_range, &TotalLocalPoints);
  PluralizeLabel(TotalLocalPoints, "points", "point", TotalLocalPointsString);

  char NeighborRanksString[256];
  int Offset = 0;
  for (i = 0; i < Grid->properties.num_neighbors; ++i) {
    char NeighborRankString[NUMBER_STRING_LENGTH];
    IntToString(Grid->properties.neighbor_ranks[i], NeighborRankString);
    Offset += sprintf(NeighborRanksString+Offset, "%s", NeighborRankString);
    if (i != Grid->properties.num_neighbors-1) {
      Offset += sprintf(NeighborRanksString+Offset, ", ");
    }
  }

  for (i = 0; i < Grid->properties.comm_size; ++i) {

    if (Grid->properties.comm_rank == i) {

      LogStatus(Grid->logger, Grid->properties.comm_rank == 0, 0, "Grid %s decomposition info:",
        Grid->properties.name);

      switch (Grid->properties.num_dims) {
      case 2:
        LogStatus(Grid->logger, true, 1, "Rank %s contains i=%s:%s, j=%s:%s (%s)",
          RankString, ILocalBeginString, ILocalEndString, JLocalBeginString, JLocalEndString,
          TotalLocalPointsString);
        break;
      case 3:
        LogStatus(Grid->logger, true, 1, "Rank %s contains i=%s:%s, j=%s:%s, k=%s:%s (%s)",
          RankString, ILocalBeginString, ILocalEndString, JLocalBeginString, JLocalEndString,
          KLocalBeginString, KLocalEndString, TotalLocalPointsString);
        break;
      }

      if (Grid->properties.num_neighbors > 0) {
        LogStatus(Grid->logger, true, 1, "Rank %s has neighbors: %s", RankString,
          NeighborRanksString);
      }

    }

    MPI_Barrier(Grid->properties.comm);

  }

}

void CreateGridParams(ovk_grid_params **Params_, int NumDims, MPI_Comm DefaultComm) {

  *Params_ = malloc(sizeof(ovk_grid_params));
  ovk_grid_params *Params = *Params_;

  Params->id = -1;

  memset(Params->name, 0, OVK_NAME_LENGTH);

  Params->num_dims = NumDims;
  Params->comm = DefaultComm;
  Params->size[0] = 0;
  Params->size[1] = 0;
  Params->size[2] = 1;
  Params->periodic[0] = false;
  Params->periodic[1] = false;
  Params->periodic[2] = false;
  Params->periodic_storage = OVK_NO_OVERLAP_PERIODIC;
  Params->periodic_length[0] = 0.;
  Params->periodic_length[1] = 0.;
  Params->periodic_length[2] = 0.;
  Params->geometry_type = OVK_GEOMETRY_TYPE_CURVILINEAR;

  ovkRangeDefault(&Params->local_range, NumDims);

  Params->num_neighbors = 0;
  Params->neighbor_ranks = NULL;

}

void DestroyGridParams(ovk_grid_params **Params_) {

  ovk_grid_params *Params = *Params_;

  free(Params->neighbor_ranks);

  free(*Params_);
  *Params_ = NULL;

}

void ovkGetGridParamID(const ovk_grid_params *Params, int *ID) {

  *ID = Params->id;

}

void ovkSetGridParamID(ovk_grid_params *Params, int ID) {

  Params->id = ID;

}

void ovkGetGridParamName(const ovk_grid_params *Params, char *Name) {

  strcpy(Name, Params->name);

}

void ovkSetGridParamName(ovk_grid_params *Params, const char *Name) {

  strncpy(Params->name, Name, OVK_NAME_LENGTH);

}

void ovkGetGridParamComm(const ovk_grid_params *Params, MPI_Comm *Comm) {

  *Comm = Params->comm;

}

void ovkSetGridParamComm(ovk_grid_params *Params, MPI_Comm Comm) {

  Params->comm = Comm;

}

void ovkGetGridParamSize(const ovk_grid_params *Params, int *Size) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    Size[i] = Params->size[i];
  }

}

void ovkSetGridParamSize(ovk_grid_params *Params, const int *Size) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    Params->size[i] = Size[i];
  }

}

void ovkGetGridParamPeriodic(const ovk_grid_params *Params, bool *Periodic) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    Periodic[i] = Params->periodic[i];
  }

}

void ovkSetGridParamPeriodic(ovk_grid_params *Params, const bool *Periodic) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    Params->periodic[i] = Periodic[i];
  }

}

void ovkGetGridParamPeriodicStorage(const ovk_grid_params *Params,
  ovk_periodic_storage *PeriodicStorage) {

  *PeriodicStorage = Params->periodic_storage;

}

void ovkSetGridParamPeriodicStorage(ovk_grid_params *Params, ovk_periodic_storage PeriodicStorage) {

  Params->periodic_storage = PeriodicStorage;

}

void ovkGetGridParamPeriodicLength(const ovk_grid_params *Params, double *PeriodicLength) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    PeriodicLength[i] = Params->periodic_length[i];
  }

}

void ovkSetGridParamPeriodicLength(ovk_grid_params *Params, const double *PeriodicLength) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    Params->periodic_length[i] = PeriodicLength[i];
  }

}

void ovkGetGridParamGeometryType(const ovk_grid_params *Params, ovk_geometry_type *GeometryType) {

  *GeometryType = Params->geometry_type;

}

void ovkSetGridParamGeometryType(ovk_grid_params *Params, ovk_geometry_type GeometryType) {

  Params->geometry_type = GeometryType;

}

void ovkGetGridParamLocalBegin(const ovk_grid_params *Params, int *LocalBegin) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    LocalBegin[i] = Params->local_range.b[i];
  }

}

void ovkSetGridParamLocalBegin(ovk_grid_params *Params, const int *LocalBegin) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    Params->local_range.b[i] = LocalBegin[i];
  }

}

void ovkGetGridParamLocalEnd(const ovk_grid_params *Params, int *LocalEnd) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    LocalEnd[i] = Params->local_range.e[i];
  }

}

void ovkSetGridParamLocalEnd(ovk_grid_params *Params, const int *LocalEnd) {

  int i;

  for (i = 0; i < Params->num_dims; ++i) {
    Params->local_range.e[i] = LocalEnd[i];
  }

}

void ovkGetGridParamLocalRange(const ovk_grid_params *Params, ovk_range *LocalRange) {

  *LocalRange = Params->local_range;

}

void ovkSetGridParamLocalRange(ovk_grid_params *Params, const ovk_range *LocalRange) {

  Params->local_range = *LocalRange;

}

void ovkGetGridParamNumNeighborRanks(ovk_grid_params *Params, int *NumNeighbors) {

  *NumNeighbors = Params->num_neighbors;

}

void ovkGetGridParamNeighborRanks(ovk_grid_params *Params, int *NeighborRanks) {

  int i;

  for (i = 0; i < Params->num_neighbors; ++i) {
    NeighborRanks[i] = Params->neighbor_ranks[i];
  }

}

void ovkSetGridParamNeighborRanks(ovk_grid_params *Params, int NumNeighbors,
  const int *NeighborRanks) {

  int i;

  free(Params->neighbor_ranks);

  Params->num_neighbors = NumNeighbors;

  if (NumNeighbors > 0) {
    Params->neighbor_ranks = malloc(NumNeighbors*sizeof(int));
    for (i = 0; i < NumNeighbors; ++i) {
      Params->neighbor_ranks[i] = NeighborRanks[i];
    }
  } else {
    Params->neighbor_ranks = NULL;
  }

}

static void DefaultGridProperties(ovk_grid_properties *Properties) {

  Properties->id = -1;

  memset(Properties->name, 0, OVK_NAME_LENGTH);

  Properties->num_dims = 2;
  Properties->comm = MPI_COMM_NULL;
  Properties->comm_size = 0;
  Properties->comm_rank = 0;
  Properties->size[0] = 0;
  Properties->size[1] = 0;
  Properties->size[2] = 1;
  Properties->periodic[0] = false;
  Properties->periodic[1] = false;
  Properties->periodic[2] = false;
  Properties->periodic_storage = OVK_NO_OVERLAP_PERIODIC;
  Properties->periodic_length[0] = 0.;
  Properties->periodic_length[1] = 0.;
  Properties->periodic_length[2] = 0.;
  Properties->geometry_type = OVK_GEOMETRY_TYPE_CURVILINEAR;

  ovkRangeDefault(&Properties->local_range, 2);

  Properties->num_neighbors = 0;
  Properties->neighbor_ranks = NULL;

}

void ovkGetGridPropertyID(const ovk_grid_properties *Properties, int *ID) {

  *ID = Properties->id;

}

void ovkGetGridPropertyName(const ovk_grid_properties *Properties, char *Name) {

  strcpy(Name, Properties->name);

}

void ovkGetGridPropertyDimension(const ovk_grid_properties *Properties, int *NumDims) {

  *NumDims = Properties->num_dims;

}

void ovkGetGridPropertyComm(const ovk_grid_properties *Properties, MPI_Comm *Comm) {

  *Comm = Properties->comm;

}

void ovkGetGridPropertySize(const ovk_grid_properties *Properties, int *Size) {

  int i;

  for (i = 0; i < Properties->num_dims; ++i) {
    Size[i] = Properties->size[i];
  }

}

void ovkGetGridPropertyPeriodic(const ovk_grid_properties *Properties, bool *Periodic) {

  int i;

  for (i = 0; i < Properties->num_dims; ++i) {
    Periodic[i] = Properties->periodic[i];
  }

}

void ovkGetGridPropertyPeriodicStorage(const ovk_grid_properties *Properties,
  ovk_periodic_storage *PeriodicStorage) {

  *PeriodicStorage = Properties->periodic_storage;

}

void ovkGetGridPropertyPeriodicLength(const ovk_grid_properties *Properties,
  double *PeriodicLength) {

  int i;

  for (i = 0; i < Properties->num_dims; ++i) {
    PeriodicLength[i] = Properties->periodic_length[i];
  }

}

void ovkGetGridPropertyGeometryType(const ovk_grid_properties *Properties,
  ovk_geometry_type *GeometryType) {

  *GeometryType = Properties->geometry_type;

}

void ovkGetGridPropertyLocalBegin(const ovk_grid_properties *Properties, int *LocalBegin) {

  int i;

  for (i = 0; i < Properties->num_dims; ++i) {
    LocalBegin[i] = Properties->local_range.b[i];
  }

}

void ovkGetGridPropertyLocalEnd(const ovk_grid_properties *Properties, int *LocalEnd) {

  int i;

  for (i = 0; i < Properties->num_dims; ++i) {
    LocalEnd[i] = Properties->local_range.e[i];
  }

}

void ovkGetGridPropertyLocalRange(const ovk_grid_properties *Properties, ovk_range *LocalRange) {

  *LocalRange = Properties->local_range;

}

void ovkGetGridPropertyNumNeighbors(ovk_grid_properties *Properties, int *NumNeighbors) {

  *NumNeighbors = Properties->num_neighbors;

}

void ovkGetGridPropertyNeighborRanks(ovk_grid_properties *Properties, int *NeighborRanks) {

  int i;

  for (i = 0; i < Properties->num_neighbors; ++i) {
    NeighborRanks[i] = Properties->neighbor_ranks[i];
  }

}
