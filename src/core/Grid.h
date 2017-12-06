// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GRID_INCLUDED
#define OVK_CORE_GRID_INCLUDED

#include "ovkGrid.h"

#include "Cart.h"
#include "ErrorHandler.h"
#include "Global.h"
#include "Logger.h"
#include "Range.h"

struct ovk_grid_params {
  int id;
  char name[OVK_NAME_LENGTH];
  int num_dims;
  MPI_Comm comm;
  int size[MAX_DIMS];
  bool periodic[MAX_DIMS];
  ovk_periodic_storage periodic_storage;
  double periodic_length[MAX_DIMS];
  ovk_geometry_type geometry_type;
  ovk_range local_range;
  int num_neighbors;
  int *neighbor_ranks;
};

struct ovk_grid_properties {
  int id;
  char name[OVK_NAME_LENGTH];
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  int size[MAX_DIMS];
  bool periodic[MAX_DIMS];
  ovk_periodic_storage periodic_storage;
  double periodic_length[MAX_DIMS];
  ovk_geometry_type geometry_type;
  ovk_range local_range;
  int num_neighbors;
  int *neighbor_ranks;
};

typedef struct {
  int comm_rank;
  ovk_range local_range;
} t_grid_neighbor_info;

struct ovk_grid {
  ovk_grid_properties properties;
  t_logger *logger;
  t_error_handler *error_handler;
  t_grid_neighbor_info *neighbors;
  ovk_cart cart;
};

void CreateGridParams(ovk_grid_params **Params, int NumDims, MPI_Comm DefaultComm);
void DestroyGridParams(ovk_grid_params **Params);

void CreateGrid(ovk_grid **Grid, int ID, const ovk_grid_params *Params, t_logger *Logger,
  t_error_handler *ErrorHandler);
void DestroyGrid(ovk_grid **Grid);

#endif