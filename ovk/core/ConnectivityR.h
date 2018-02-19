// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONNECTIVITY_R_INCLUDED
#define OVK_CORE_CONNECTIVITY_R_INCLUDED

#include "ovk/core/ovkConnectivityR.h"

#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/Logger.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity_r_properties {
  int grid_id;
  int source_grid_id;
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  size_t num_receivers;
};

typedef struct {
  bool num_receivers;
  bool points;
  bool sources;
} t_connectivity_r_edits;

struct ovk_connectivity_r {
  ovk_connectivity_r_properties properties;
  int properties_edit_ref_count;
  t_logger *logger;
  t_error_handler *error_handler;
  const ovk_grid *grid;
  t_connectivity_r_edits edits;
  int *points[MAX_DIMS];
  int points_edit_ref_count;
  int *sources[MAX_DIMS];
  int sources_edit_ref_count;
  int *source_ranks;
  int source_ranks_edit_ref_count;
};

void PRIVATE(CreateConnectivityReceiverSide)(ovk_connectivity_r **Receivers,
  const ovk_grid *Grid, int SourceGridID, t_logger *Logger, t_error_handler *ErrorHandler);
#define CreateConnectivityReceiverSide(...) PRIVATE(CreateConnectivityReceiverSide)(__VA_ARGS__)
void PRIVATE(DestroyConnectivityReceiverSide)(ovk_connectivity_r **Receivers);
#define DestroyConnectivityReceiverSide(...) PRIVATE(DestroyConnectivityReceiverSide)(__VA_ARGS__)

void PRIVATE(GetConnectivityReceiverSideEdits)(const ovk_connectivity_r *Receivers,
  const t_connectivity_r_edits **Edits);
#define GetConnectivityReceiverSideEdits(...) PRIVATE(GetConnectivityReceiverSideEdits)(__VA_ARGS__)
void PRIVATE(ResetConnectivityReceiverSideEdits)(ovk_connectivity_r *Receivers);
#define ResetConnectivityReceiverSideEdits(...) PRIVATE(ResetConnectivityReceiverSideEdits)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
