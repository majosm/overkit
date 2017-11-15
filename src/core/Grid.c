// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "Grid.h"

#include "Debug.h"
#include "ErrorHandler.h"
#include "Global.h"
#include "Logger.h"
#include "TextUtils.h"

static void PrintGridSummary(const ovk_grid *Grid);
static void PrintGridDecomposition(const ovk_grid *Grid);

static void CreateGridProperties(ovk_grid_properties **Properties_);
static void DestroyGridProperties(ovk_grid_properties **Properties_);

void CreateGrid(ovk_grid **Grid_, int ID, const ovk_grid_params *Params, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  int i;

  MPI_Comm Comm;
  MPI_Comm_dup(Params->comm, &Comm);

  MPI_Barrier(Comm);

  *Grid_ = malloc(sizeof(ovk_grid));
  ovk_grid *Grid = *Grid_;

  CreateGridProperties(&Grid->properties);

  Grid->properties->id = ID;
  Grid->properties->num_dims = Params->num_dims;

  for (i = 0; i < MAX_DIMS; ++i) {
    Grid->properties->global_size[i] = Params->global_size[i];
    Grid->properties->local_start[i] = Params->local_start[i];
    Grid->properties->local_end[i] = Params->local_end[i];
    Grid->properties->periodic[i] = Params->periodic[i];
    Grid->properties->periodic_length[i] = Params->periodic_length[i];
  }

  Grid->properties->periodic_storage = Params->periodic_storage;
  Grid->properties->geometry_type = Params->geometry_type;

  Grid->properties->comm = Comm;
  MPI_Comm_size(Grid->properties->comm, &Grid->properties->comm_size);
  MPI_Comm_rank(Grid->properties->comm, &Grid->properties->comm_rank);

  Grid->properties->num_neighbors = Params->num_neighbors;
  if (Params->num_neighbors > 0) {
    Grid->properties->neighbor_ranks = malloc(Params->num_neighbors*sizeof(int));
  }

  for (i = 0; i < Params->num_neighbors; ++i) {
    Grid->properties->neighbor_ranks[i] = Params->neighbor_ranks[i];
  }

  Grid->logger = Logger;
  Grid->error_handler = ErrorHandler;

  if (Grid->properties->comm_rank == 0) {
    PrintGridSummary(Grid);
  }

  if (OVK_DEBUG) {
    PrintGridDecomposition(Grid);
  }

  MPI_Barrier(Grid->properties->comm);

}

void DestroyGrid(ovk_grid **Grid_) {

  ovk_grid *Grid = *Grid_;

  MPI_Barrier(Grid->properties->comm);

  t_logger *Logger = Grid->logger;
  MPI_Comm Comm = Grid->properties->comm;
  bool IsGridRoot = Grid->properties->comm_rank == 0;
  int ID = Grid->properties->id;

  DestroyGridProperties(&Grid->properties);

  free(*Grid_);
  *Grid_ = NULL;

  MPI_Barrier(Comm);

  MPI_Comm_free(&Comm);

  LogStatus(Logger, IsGridRoot, 0, "Destroyed grid %i.", ID);

}

void ovkGetGridProperties(const ovk_grid *Grid, const ovk_grid_properties **Properties) {

  *Properties = Grid->properties;

}

static void PrintGridSummary(const ovk_grid *Grid) {

  size_t TotalPoints =
    (size_t)Grid->properties->global_size[0] *
    (size_t)Grid->properties->global_size[1] *
    (size_t)Grid->properties->global_size[2];

  char TotalPointsString[NUMBER_STRING_LENGTH+7];
  PluralizeLabel(TotalPoints, "points", "point", TotalPointsString);

  char ProcessesString[NUMBER_STRING_LENGTH+10];
  PluralizeLabel(Grid->properties->comm_size, "processes", "process", ProcessesString);

  char ISizeString[NUMBER_STRING_LENGTH];
  char JSizeString[NUMBER_STRING_LENGTH];
  char KSizeString[NUMBER_STRING_LENGTH];
  SizeToString(Grid->properties->global_size[0], ISizeString);
  SizeToString(Grid->properties->global_size[1], JSizeString);
  SizeToString(Grid->properties->global_size[2], KSizeString);

  switch (Grid->properties->num_dims) {
  case 2:
    LogStatus(Grid->logger, true, 0, "Created grid %i: %s x %s (%s) on %s.",
      Grid->properties->id, ISizeString, JSizeString, TotalPointsString, ProcessesString);
    break;
  case 3:
    LogStatus(Grid->logger, true, 0, "Created grid %i: %s x %s x %s (%s) on %s.",
      Grid->properties->id, ISizeString, JSizeString, KSizeString, TotalPointsString,
      ProcessesString);
    break;
  }

}

static void PrintGridDecomposition(const ovk_grid *Grid) {

  int i;

  char IDString[NUMBER_STRING_LENGTH];
  IntToString(Grid->properties->id, IDString);

  char RankString[NUMBER_STRING_LENGTH];
  IntToString(Grid->properties->comm_rank, RankString);

  char ILocalStartString[NUMBER_STRING_LENGTH], ILocalEndString[NUMBER_STRING_LENGTH];
  char JLocalStartString[NUMBER_STRING_LENGTH], JLocalEndString[NUMBER_STRING_LENGTH];
  char KLocalStartString[NUMBER_STRING_LENGTH], KLocalEndString[NUMBER_STRING_LENGTH];
  IntToString(Grid->properties->local_start[0], ILocalStartString);
  IntToString(Grid->properties->local_start[1], JLocalStartString);
  IntToString(Grid->properties->local_start[2], KLocalStartString);
  IntToString(Grid->properties->local_end[0], ILocalEndString);
  IntToString(Grid->properties->local_end[1], JLocalEndString);
  IntToString(Grid->properties->local_end[2], KLocalEndString);

  char TotalLocalPointsString[NUMBER_STRING_LENGTH+7];
  size_t TotalLocalPoints;
  switch (Grid->properties->num_dims) {
  case 2:
    TotalLocalPoints =
      (size_t)(Grid->properties->local_end[0] - Grid->properties->local_start[0]) *
      (size_t)(Grid->properties->local_end[1] - Grid->properties->local_start[1]);
    break;
  case 3:
    TotalLocalPoints =
      (size_t)(Grid->properties->local_end[0] - Grid->properties->local_start[0]) *
      (size_t)(Grid->properties->local_end[1] - Grid->properties->local_start[1]) *
      (size_t)(Grid->properties->local_end[2] - Grid->properties->local_start[2]);
    break;
  }
  PluralizeLabel(TotalLocalPoints, "points", "point", TotalLocalPointsString);

  char NeighborRanksString[256];
  int Offset = 0;
  for (int i = 0; i < Grid->properties->num_neighbors; ++i) {
    char NeighborRankString[NUMBER_STRING_LENGTH];
    IntToString(Grid->properties->neighbor_ranks[i], NeighborRankString);
    Offset += sprintf(NeighborRanksString+Offset, "%s", NeighborRankString);
    if (i != Grid->properties->num_neighbors-1) { 
      Offset += sprintf(NeighborRanksString+Offset, ", ");
    }
  }

  for (i = 0; i < Grid->properties->comm_size; ++i) {

    if (Grid->properties->comm_rank == i) {

      LogStatus(Grid->logger, Grid->properties->comm_rank == 0, 0, "Grid %s decomposition info:",
        IDString);

      switch (Grid->properties->num_dims) {
      case 2:
        LogStatus(Grid->logger, true, 1, "Rank %s contains i=%s:%s, j=%s:%s (%s)",
          RankString, ILocalStartString, ILocalEndString, JLocalStartString, JLocalEndString,
          TotalLocalPointsString);
        break;
      case 3:
        LogStatus(Grid->logger, true, 1, "Rank %s contains i=%s:%s, j=%s:%s, k=%s:%s (%s)",
          RankString, ILocalStartString, ILocalEndString, JLocalStartString, JLocalEndString,
          KLocalStartString, KLocalEndString, TotalLocalPointsString);
        break;
      }

      if (Grid->properties->num_neighbors > 0) {
        LogStatus(Grid->logger, true, 1, "Rank %s has neighbors: %s", RankString,
          NeighborRanksString);
      }

    }

    MPI_Barrier(Grid->properties->comm);

  }

}

void CreateGridParams(ovk_grid_params **Params_, int NumDims, MPI_Comm DefaultComm) {

  *Params_ = malloc(sizeof(ovk_grid_params));
  ovk_grid_params *Params = *Params_;

  Params->num_dims = NumDims;
  Params->global_size[0] = 0;
  Params->global_size[1] = 0;
  Params->global_size[2] = 1;
  Params->local_start[0] = 0;
  Params->local_start[1] = 0;
  Params->local_start[2] = 1;
  Params->local_end[0] = 0;
  Params->local_end[1] = 0;
  Params->local_end[2] = 1;
  Params->periodic[0] = false;
  Params->periodic[1] = false;
  Params->periodic[2] = false;
  Params->periodic_storage = OVK_NO_OVERLAP_PERIODIC;
  Params->periodic_length[0] = 0.;
  Params->periodic_length[1] = 0.;
  Params->periodic_length[2] = 0.;
  Params->geometry_type = OVK_GEOMETRY_TYPE_CURVILINEAR;
  Params->comm = DefaultComm;
  Params->num_neighbors = 0;
  Params->neighbor_ranks = NULL;

}

void DestroyGridParams(ovk_grid_params **Params_) {

  ovk_grid_params *Params = *Params_;

  if (Params->num_neighbors > 0) {
    free(Params->neighbor_ranks);
  }

  free(*Params_);
  *Params_ = NULL;

}

void ovkGetGridParamID(const ovk_grid_params *Params, int *ID) {

  *ID = Params->id;

}

void ovkSetGridParamID(ovk_grid_params *Params, int ID) {

  Params->id = ID;

}

void ovkGetGridParamGlobalSize(const ovk_grid_params *Params, int *GlobalSize) {

  for (int i = 0; i < Params->num_dims; ++i) {
    GlobalSize[i] = Params->global_size[i];
  }

}

void ovkSetGridParamGlobalSize(ovk_grid_params *Params, const int *GlobalSize) {

  for (int i = 0; i < Params->num_dims; ++i) {
    Params->global_size[i] = GlobalSize[i];
  }

}

void ovkGetGridParamLocalStart(const ovk_grid_params *Params, int *LocalStart) {

  for (int i = 0; i < Params->num_dims; ++i) {
    LocalStart[i] = Params->local_start[i];
  }

}

void ovkSetGridParamLocalStart(ovk_grid_params *Params, const int *LocalStart) {

  for (int i = 0; i < Params->num_dims; ++i) {
    Params->local_start[i] = LocalStart[i];
  }

}

void ovkGetGridParamLocalEnd(const ovk_grid_params *Params, int *LocalEnd) {

  for (int i = 0; i < Params->num_dims; ++i) {
    LocalEnd[i] = Params->local_end[i];
  }

}

void ovkSetGridParamLocalEnd(ovk_grid_params *Params, const int *LocalEnd) {

  for (int i = 0; i < Params->num_dims; ++i) {
    Params->local_end[i] = LocalEnd[i];
  }

}

void ovkGetGridParamPeriodic(const ovk_grid_params *Params, bool *Periodic) {

  for (int i = 0; i < Params->num_dims; ++i) {
    Periodic[i] = Params->periodic[i];
  }

}

void ovkSetGridParamPeriodic(ovk_grid_params *Params, const bool *Periodic) {

  for (int i = 0; i < Params->num_dims; ++i) {
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

  for (int i = 0; i < Params->num_dims; ++i) {
    PeriodicLength[i] = Params->periodic_length[i];
  }

}

void ovkSetGridParamPeriodicLength(ovk_grid_params *Params, const double *PeriodicLength) {

  for (int i = 0; i < Params->num_dims; ++i) {
    Params->periodic_length[i] = PeriodicLength[i];
  }

}

void ovkGetGridParamGeometryType(const ovk_grid_params *Params, ovk_geometry_type *GeometryType) {

  *GeometryType = Params->geometry_type;

}

void ovkSetGridParamGeometryType(ovk_grid_params *Params, ovk_geometry_type GeometryType) {

  Params->geometry_type = GeometryType;

}

void ovkGetGridParamComm(const ovk_grid_params *Params, MPI_Comm *Comm) {

  *Comm = Params->comm;

}

void ovkSetGridParamComm(ovk_grid_params *Params, MPI_Comm Comm) {

  Params->comm = Comm;

}

void ovkGetGridParamNumNeighborRanks(ovk_grid_params *Params, int *NumNeighbors) {

  *NumNeighbors = Params->num_neighbors;

}

void ovkGetGridParamNeighborRanks(ovk_grid_params *Params, int *NeighborRanks) {

  for (int i = 0; i < Params->num_neighbors; ++i) {
    NeighborRanks[i] = Params->neighbor_ranks[i];
  }

}

void ovkSetGridParamNeighborRanks(ovk_grid_params *Params, int NumNeighbors,
  const int *NeighborRanks) {

  if (Params->num_neighbors > 0) {
    free(Params->neighbor_ranks);
  }

  Params->num_neighbors = NumNeighbors;

  if (NumNeighbors > 0) {
    Params->neighbor_ranks = malloc(NumNeighbors*sizeof(int));
    for (int i = 0; i < NumNeighbors; ++i) {
      Params->neighbor_ranks[i] = NeighborRanks[i];
    }
  } else {
    Params->neighbor_ranks = NULL;
  }

}

static void CreateGridProperties(ovk_grid_properties **Properties_) {
  
  *Properties_ = malloc(sizeof(ovk_grid_properties));
  ovk_grid_properties *Properties = *Properties_;

  Properties->id = -1;
  Properties->num_dims = 2;
  Properties->global_size[0] = 0;
  Properties->global_size[1] = 0;
  Properties->global_size[2] = 1;
  Properties->local_start[0] = 0;
  Properties->local_start[1] = 0;
  Properties->local_start[2] = 1;
  Properties->local_end[0] = 0;
  Properties->local_end[1] = 0;
  Properties->local_end[2] = 1;
  Properties->periodic[0] = false;
  Properties->periodic[1] = false;
  Properties->periodic[2] = false;
  Properties->periodic_storage = OVK_NO_OVERLAP_PERIODIC;
  Properties->periodic_length[0] = 0.;
  Properties->periodic_length[1] = 0.;
  Properties->periodic_length[2] = 0.;
  Properties->geometry_type = OVK_GEOMETRY_TYPE_CURVILINEAR;
  Properties->comm = MPI_COMM_NULL;
  Properties->comm_size = 0;
  Properties->comm_rank = 0;

}

static void DestroyGridProperties(ovk_grid_properties **Properties_) {

  ovk_grid_properties *Properties = *Properties_;

  if (Properties->num_neighbors > 0) {
    free(Properties->neighbor_ranks);
  }

  free(*Properties_);
  *Properties_ = NULL;

}

void ovkGetGridPropertyDimension(const ovk_grid_properties *Properties, int *NumDims) {

  *NumDims = Properties->num_dims;

}

void ovkGetGridPropertyGlobalSize(const ovk_grid_properties *Properties, int *GlobalSize) {

  for (int i = 0; i < Properties->num_dims; ++i) {
    GlobalSize[i] = Properties->global_size[i];
  }

}

void ovkGetGridPropertyLocalStart(const ovk_grid_properties *Properties, int *LocalStart) {

  for (int i = 0; i < Properties->num_dims; ++i) {
    LocalStart[i] = Properties->local_start[i];
  }

}

void ovkGetGridPropertyLocalEnd(const ovk_grid_properties *Properties, int *LocalEnd) {

  for (int i = 0; i < Properties->num_dims; ++i) {
    LocalEnd[i] = Properties->local_end[i];
  }

}

void ovkGetGridPropertyPeriodic(const ovk_grid_properties *Properties, bool *Periodic) {

  for (int i = 0; i < Properties->num_dims; ++i) {
    Periodic[i] = Properties->periodic[i];
  }

}

void ovkGetGridPropertyPeriodicStorage(const ovk_grid_properties *Properties,
  ovk_periodic_storage *PeriodicStorage) {

  *PeriodicStorage = Properties->periodic_storage;

}

void ovkGetGridPropertyPeriodicLength(const ovk_grid_properties *Properties,
  double *PeriodicLength) {

  for (int i = 0; i < Properties->num_dims; ++i) {
    PeriodicLength[i] = Properties->periodic_length[i];
  }

}

void ovkGetGridPropertyGeometryType(const ovk_grid_properties *Properties,
  ovk_geometry_type *GeometryType) {

  *GeometryType = Properties->geometry_type;

}

void ovkGetGridPropertyComm(const ovk_grid_properties *Properties, MPI_Comm *Comm) {

  *Comm = Properties->comm;

}

void ovkGetGridPropertyNumNeighbors(ovk_grid_properties *Properties, int *NumNeighbors) {

  *NumNeighbors = Properties->num_neighbors;

}

void ovkGetGridPropertyNeighborRanks(ovk_grid_properties *Properties, int *NeighborRanks) {

  for (int i = 0; i < Properties->num_neighbors; ++i) {
    NeighborRanks[i] = Properties->neighbor_ranks[i];
  }

}
