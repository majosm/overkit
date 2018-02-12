// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GRID_INCLUDED
#define OVK_CORE_GRID_INCLUDED

#include "ovk/core/ovkGrid.h"

#include "ovk/core/Cart.h"
#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Logger.h"
#include "ovk/core/PartitionHash.h"
#include "ovk/core/Range.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_grid_params {
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
  bool periodic[MAX_DIMS];
  ovk_periodic_storage periodic_storage;
  double periodic_length[MAX_DIMS];
  ovk_geometry_type geometry_type;
  ovk_range global_range;
  ovk_range local_range;
  int num_neighbors;
  int *neighbor_ranks;
};

typedef struct {
  int comm_rank;
  ovk_range local_range;
} t_grid_neighbor;

struct ovk_grid {
  ovk_grid_properties properties;
  t_logger *logger;
  t_error_handler *error_handler;
  ovk_cart cart;
  t_grid_neighbor *neighbors;
  t_partition_hash *partition_hash;
};

struct ovk_grid_info {
  int id;
  char name[OVK_NAME_LENGTH];
  int num_dims;
  int root_rank;
  ovk_range global_range;
  ovk_cart cart;
  double periodic_length[MAX_DIMS];
  ovk_geometry_type geometry_type;
};

void PRIVATE(CreateGrid)(ovk_grid **Grid, int ID, const ovk_grid_params *Params, t_logger *Logger,
  t_error_handler *ErrorHandler);
#define CreateGrid(...) PRIVATE(CreateGrid)(__VA_ARGS__)
void PRIVATE(DestroyGrid)(ovk_grid **Grid);
#define DestroyGrid(...) PRIVATE(DestroyGrid)(__VA_ARGS__)

static inline void GetGridLogger(const ovk_grid *Grid, t_logger **Logger) {
  *Logger = (t_logger *)Grid->logger;
}
static inline void GetGridErrorHandler(const ovk_grid *Grid, t_error_handler **ErrorHandler) {
  *ErrorHandler = (t_error_handler *)Grid->error_handler;
}

void PRIVATE(GetGridNeighborRange)(const ovk_grid *Grid, int iNeighbor, ovk_range *NeighborRange);
#define GetGridNeighborRange(...) PRIVATE(GetGridNeighborRange)(__VA_ARGS__)
void PRIVATE(GetGridNeighborRank)(const ovk_grid *Grid, int iNeighbor, int *NeighborRank);
#define GetGridNeighborRank(...) PRIVATE(GetGridNeighborRank)(__VA_ARGS__)

static inline void GetGridPartitionHash(const ovk_grid *Grid, const t_partition_hash **Hash) {
  *Hash = Grid->partition_hash;
}

void PRIVATE(CreateGridParams)(ovk_grid_params **Params, int NumDims, MPI_Comm DefaultComm);
#define CreateGridParams(...) PRIVATE(CreateGridParams)(__VA_ARGS__)
void PRIVATE(DestroyGridParams)(ovk_grid_params **Params);
#define DestroyGridParams(...) PRIVATE(DestroyGridParams)(__VA_ARGS__)

void PRIVATE(CreateGridInfo)(ovk_grid_info **Info, const ovk_grid *Grid, MPI_Comm Comm,
  int CommRank);
#define CreateGridInfo(...) PRIVATE(CreateGridInfo)(__VA_ARGS__)
void PRIVATE(DestroyGridInfo)(ovk_grid_info **Info);
#define DestroyGridInfo(...) PRIVATE(DestroyGridInfo)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
