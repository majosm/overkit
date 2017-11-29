// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GRID_INCLUDED
#define OVK_CORE_GRID_INCLUDED

#include "ovkGrid.h"

#include "ErrorHandler.h"
#include "Global.h"
#include "Logger.h"

struct ovk_grid_params {
  int id;
  char name[OVK_NAME_LENGTH];
  int num_dims;
  int global_size[MAX_DIMS];
  int local_start[MAX_DIMS];
  int local_end[MAX_DIMS];
  bool periodic[MAX_DIMS];
  ovk_periodic_storage periodic_storage;
  double periodic_length[MAX_DIMS];
  ovk_geometry_type geometry_type;
  MPI_Comm comm;
  int num_neighbors;
  int *neighbor_ranks;
};

typedef struct {
  int comm_rank;
  int local_start[MAX_DIMS];
  int local_end[MAX_DIMS];
} t_grid_neighbor_info;

struct ovk_grid {
  ovk_grid_properties *properties;
  t_logger *logger;
  t_error_handler *error_handler;
  t_grid_neighbor_info *neighbors;
};

struct ovk_grid_properties {
  int id;
  char name[OVK_NAME_LENGTH];
  int num_dims;
  int global_size[MAX_DIMS];
  int local_start[MAX_DIMS];
  int local_end[MAX_DIMS];
  bool periodic[MAX_DIMS];
  ovk_periodic_storage periodic_storage;
  double periodic_length[MAX_DIMS];
  ovk_geometry_type geometry_type;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  int num_neighbors;
  int *neighbor_ranks;
};

void CreateGridParams(ovk_grid_params **Params, int NumDims, MPI_Comm DefaultComm);
void DestroyGridParams(ovk_grid_params **Params);

void CreateGrid(ovk_grid **Grid, int ID, const ovk_grid_params *Params, t_logger *Logger,
  t_error_handler *ErrorHandler);
void DestroyGrid(ovk_grid **Grid);

#endif
