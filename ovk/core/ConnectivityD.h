// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONNECTIVITY_D_INCLUDED
#define OVK_CORE_CONNECTIVITY_D_INCLUDED

#include "ovk/core/ovkConnectivityD.h"

#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/Logger.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity_d_properties {
  int grid_id;
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  size_t num_donors;
  int max_donor_size;
};

typedef struct {
  bool num_donors;
  bool extents;
  bool coords;
  bool interp_coefs;
  bool destinations;
} t_connectivity_d_edits;

struct ovk_connectivity_d {
  ovk_connectivity_d_properties properties;
  int properties_edit_ref_count;
  t_logger *logger;
  t_error_handler *error_handler;
  const ovk_grid *grid;
  t_connectivity_d_edits edits;
  int *extents[2][MAX_DIMS];
  int extents_edit_ref_count;
  double *coords[MAX_DIMS];
  int coords_edit_ref_count;
  double **interp_coefs[MAX_DIMS];
  int interp_coefs_edit_ref_count;
  int *destinations[MAX_DIMS];
  int destinations_edit_ref_count;
  int *destination_ranks;
  int destination_ranks_edit_ref_count;
};

void PRIVATE(CreateConnectivityDonorSide)(ovk_connectivity_d **Donors, const ovk_grid *Grid,
  t_logger *Logger, t_error_handler *ErrorHandler);
#define CreateConnectivityDonorSide(...) PRIVATE(CreateConnectivityDonorSide)(__VA_ARGS__)
void PRIVATE(DestroyConnectivityDonorSide)(ovk_connectivity_d **Donors);
#define DestroyConnectivityDonorSide(...) PRIVATE(DestroyConnectivityDonorSide)(__VA_ARGS__)

void PRIVATE(GetConnectivityDonorSideEdits)(const ovk_connectivity_d *Donors,
  const t_connectivity_d_edits **Edits);
#define GetConnectivityDonorSideEdits(...) PRIVATE(GetConnectivityDonorSideEdits)(__VA_ARGS__)
void PRIVATE(ResetConnectivityDonorSideEdits)(ovk_connectivity_d *Donors);
#define ResetConnectivityDonorSideEdits(...) PRIVATE(ResetConnectivityDonorSideEdits)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
