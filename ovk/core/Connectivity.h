// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONNECTIVITY_INCLUDED
#define OVK_CORE_CONNECTIVITY_INCLUDED

#include "ovk/core/ovkConnectivity.h"

#include "ovk/core/ConnectivityD.h"
#include "ovk/core/ConnectivityR.h"
#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/Logger.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  bool num_donors;
  bool donor_extents;
  bool donor_coords;
  bool donor_interp_coefs;
  bool donor_destinations;
  bool num_receivers;
  bool receiver_points;
  bool receiver_sources;
} t_connectivity_edits;

typedef struct {
  ovk_grid_info *grid_info;
  int edit_ref_count;
  bool has_local_data;
  const ovk_grid *grid;
  ovk_connectivity_d *donors;
} t_connectivity_donor_side_container;

typedef struct {
  ovk_grid_info *grid_info;
  int edit_ref_count;
  bool has_local_data;
  const ovk_grid *grid;
  ovk_connectivity_r *receivers;
} t_connectivity_receiver_side_container;

struct ovk_connectivity {
  t_logger *logger;
  t_error_handler *error_handler;
  int donor_grid_id;
  int receiver_grid_id;
  char name[OVK_NAME_LENGTH];
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  t_connectivity_edits edits;
  t_connectivity_donor_side_container *donors_container;
  t_connectivity_receiver_side_container *receivers_container;
};

struct ovk_connectivity_info {
  int donor_grid_id;
  int receiver_grid_id;
  char name[OVK_NAME_LENGTH];
  int num_dims;
  int root_rank;
};

void PRIVATE(CreateConnectivity)(ovk_connectivity **Connectivity, int NumDims, MPI_Comm Comm,
  const ovk_grid *DonorGrid, const ovk_grid *ReceiverGrid, t_logger *Logger,
  t_error_handler *ErrorHandler);
#define CreateConnectivity(...) PRIVATE(CreateConnectivity)(__VA_ARGS__)
void PRIVATE(DestroyConnectivity)(ovk_connectivity **Connectivity);
#define DestroyConnectivity(...) PRIVATE(DestroyConnectivity)(__VA_ARGS__)

void PRIVATE(CreateConnectivityInfo)(ovk_connectivity_info **Info,
  const ovk_connectivity *Connectivity, MPI_Comm Comm, int CommRank);
#define CreateConnectivityInfo(...) PRIVATE(CreateConnectivityInfo)(__VA_ARGS__)
void PRIVATE(DestroyConnectivityInfo)(ovk_connectivity_info **Info);
#define DestroyConnectivityInfo(...) PRIVATE(DestroyConnectivityInfo)(__VA_ARGS__)

void PRIVATE(GetConnectivityEdits)(const ovk_connectivity *Connectivity,
  const t_connectivity_edits **Edits);
#define GetConnectivityEdits(...) PRIVATE(GetConnectivityEdits)(__VA_ARGS__)
void PRIVATE(ResetConnectivityEdits)(ovk_connectivity *Connectivity);
#define ResetConnectivityEdits(...) PRIVATE(ResetConnectivityEdits)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
