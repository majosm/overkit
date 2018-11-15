// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ConnectivityR.h"

#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/Logger.h"

static bool EditingPoints(const ovk_connectivity_r *Receivers);
static bool EditingSources(const ovk_connectivity_r *Receivers);
static bool EditingSourceRanks(const ovk_connectivity_r *Receivers);

static void DefaultEdits(t_connectivity_r_edits *Edits);

void PRIVATE(CreateConnectivityReceiverSide)(ovk_connectivity_r **Receivers_, const ovk_grid *Grid,
  int SourceGridID, t_logger *Logger, t_error_handler *ErrorHandler) {

  int iDim;

  int GridID;
  int NumDims;
  MPI_Comm Comm;
  int CommSize, CommRank;
  ovkGetGridID(Grid, &GridID);
  ovkGetGridDimension(Grid, &NumDims);
  ovkGetGridComm(Grid, &Comm);
  ovkGetGridCommSize(Grid, &CommSize);
  ovkGetGridCommRank(Grid, &CommRank);

  MPI_Barrier(Comm);

  *Receivers_ = malloc(sizeof(ovk_connectivity_r));
  ovk_connectivity_r *Receivers = *Receivers_;

  Receivers->logger = Logger;
  Receivers->error_handler = ErrorHandler;

  Receivers->grid_id = GridID;
  Receivers->source_grid_id = SourceGridID;
  Receivers->num_dims = NumDims;
  Receivers->comm = Comm;
  Receivers->comm_size = CommSize;
  Receivers->comm_rank = CommRank;
  Receivers->count = 0;

  Receivers->grid = Grid;

  DefaultEdits(&Receivers->edits);

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers->points[iDim] = NULL;
  }
  Receivers->points_edit_ref_count = 0;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers->sources[iDim] = NULL;
  }
  Receivers->sources_edit_ref_count = 0;

  Receivers->source_ranks = NULL;
  Receivers->source_ranks_edit_ref_count = 0;

  MPI_Barrier(Receivers->comm);

}

void PRIVATE(DestroyConnectivityReceiverSide)(ovk_connectivity_r **Receivers_) {

  ovk_connectivity_r *Receivers = *Receivers_;

  MPI_Barrier(Receivers->comm);

  free(Receivers->points[0]);
  free(Receivers->sources[0]);
  free(Receivers->source_ranks);

  MPI_Comm Comm = Receivers->comm;

  free_null(Receivers_);

  MPI_Barrier(Comm);

}

void ovkGetConnectivityReceiverSideGridID(const ovk_connectivity_r *Receivers, int *GridID) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  *GridID = Receivers->grid_id;

}

void ovkGetConnectivityReceiverSideSourceGridID(const ovk_connectivity_r *Receivers, int
  *SourceGridID) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(SourceGridID, "Invalid source grid ID pointer.");

  *SourceGridID = Receivers->source_grid_id;

}

void ovkGetConnectivityReceiverSideDimension(const ovk_connectivity_r *Receivers, int *NumDims) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Receivers->num_dims;

}

void ovkGetConnectivityReceiverSideComm(const ovk_connectivity_r *Receivers, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Receivers->comm;

}

void ovkGetConnectivityReceiverSideCommSize(const ovk_connectivity_r *Receivers, int *CommSize) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  *CommSize = Receivers->comm_size;

}

void ovkGetConnectivityReceiverSideCommRank(const ovk_connectivity_r *Receivers, int *CommRank) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  *CommRank = Receivers->comm_rank;

}

void ovkGetConnectivityReceiverSideCount(const ovk_connectivity_r *Receivers, long long *NumReceivers) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(NumReceivers, "Invalid num receivers pointer.");

  *NumReceivers = Receivers->count;

}

void ovkResizeReceivers(ovk_connectivity_r *Receivers, long long NumReceivers) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(NumReceivers >= 0, "Invalid receiver count.");

  MPI_Barrier(Receivers->comm);

  OVK_DEBUG_ASSERT(!EditingPoints(Receivers), "Cannot resize receivers while editing points.");
  OVK_DEBUG_ASSERT(!EditingSources(Receivers), "Cannot resize receivers while editing sources.");
  OVK_DEBUG_ASSERT(!EditingSourceRanks(Receivers), "Cannot resize receivers while editing source "
    "ranks.");

  int iDim;
  long long iReceiver;

  free(Receivers->points[0]);
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers->points[iDim] = NULL;
  }
  free(Receivers->sources[0]);
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers->sources[iDim] = NULL;
  }
  free_null(&Receivers->source_ranks);

  Receivers->count = NumReceivers;

  if (NumReceivers > 0) {

    Receivers->points[0] = malloc(MAX_DIMS*NumReceivers*sizeof(int));
    Receivers->points[1] = Receivers->points[0] + NumReceivers;
    Receivers->points[2] = Receivers->points[1] + NumReceivers;
    Receivers->sources[0] = malloc(MAX_DIMS*NumReceivers*sizeof(int));
    Receivers->sources[1] = Receivers->sources[0] + NumReceivers;
    Receivers->sources[2] = Receivers->sources[1] + NumReceivers;
    Receivers->source_ranks = malloc(NumReceivers*sizeof(int));

    for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Receivers->points[iDim][iReceiver] = 0;
      }
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Receivers->sources[iDim][iReceiver] = 0;
      }
      Receivers->source_ranks[iReceiver] = -1;
    }

  }

  Receivers->edits.count = true;
  Receivers->edits.points = true;
  Receivers->edits.sources = true;

  MPI_Barrier(Receivers->comm);

}

void ovkGetReceiverPoints(const ovk_connectivity_r *Receivers, int Dimension, const int **Points) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  *Points = Receivers->points[Dimension];

}

void ovkEditReceiverPoints(ovk_connectivity_r *Receivers, int Dimension, int **Points) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");

  bool StartEdit = Receivers->points_edit_ref_count == 0;
  ++Receivers->points_edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Receivers->comm);
  }

  *Points = Receivers->points[Dimension];

}

void ovkReleaseReceiverPoints(ovk_connectivity_r *Receivers, int Dimension, int **Points) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points, "Invalid points pointer.");
  OVK_DEBUG_ASSERT(*Points == Receivers->points[Dimension], "Invalid points pointer.");
  OVK_DEBUG_ASSERT(EditingPoints(Receivers), "Unable to release points; not currently being "
    "edited.");

  --Receivers->points_edit_ref_count;
  bool EndEdit = Receivers->points_edit_ref_count == 0;

  *Points = NULL;

  if (EndEdit) {
    Receivers->edits.points = true;
    MPI_Barrier(Receivers->comm);
  }

}

void ovkGetReceiverSources(const ovk_connectivity_r *Receivers, int Dimension, const int **Sources)
  {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  *Sources = Receivers->sources[Dimension];

}

void ovkEditReceiverSources(ovk_connectivity_r *Receivers, int Dimension, int **Sources) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");

  bool StartEdit = Receivers->sources_edit_ref_count == 0;
  ++Receivers->sources_edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Receivers->comm);
  }

  *Sources = Receivers->sources[Dimension];

}

void ovkReleaseReceiverSources(ovk_connectivity_r *Receivers, int Dimension, int **Sources) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources, "Invalid sources pointer.");
  OVK_DEBUG_ASSERT(*Sources == Receivers->sources[Dimension], "Invalid sources pointer.");
  OVK_DEBUG_ASSERT(EditingSources(Receivers), "Unable to release sources; not currently being "
    "edited.");

  --Receivers->sources_edit_ref_count;
  bool EndEdit = Receivers->sources_edit_ref_count == 0;

  *Sources = NULL;

  if (EndEdit) {
    Receivers->edits.sources = true;
    MPI_Barrier(Receivers->comm);
  }

}

void ovkGetReceiverSourceRanks(const ovk_connectivity_r *Receivers, const int **SourceRanks) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  *SourceRanks = Receivers->source_ranks;

}

void ovkEditReceiverSourceRanks(ovk_connectivity_r *Receivers, int **SourceRanks) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");

  bool StartEdit = Receivers->source_ranks_edit_ref_count == 0;
  ++Receivers->source_ranks_edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Receivers->comm);
  }

  *SourceRanks = Receivers->source_ranks;

}

void ovkReleaseReceiverSourceRanks(ovk_connectivity_r *Receivers, int **SourceRanks) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(SourceRanks, "Invalid source ranks pointer.");
  OVK_DEBUG_ASSERT(*SourceRanks == Receivers->source_ranks, "Invalid source ranks pointer.");
  OVK_DEBUG_ASSERT(EditingSourceRanks(Receivers), "Unable to release source ranks; not currently "
    "being edited.");

  --Receivers->source_ranks_edit_ref_count;
  bool EndEdit = Receivers->source_ranks_edit_ref_count == 0;

  *SourceRanks = NULL;

  if (EndEdit) {
    Receivers->edits.sources = true;
    MPI_Barrier(Receivers->comm);
  }

}

void ovkGetConnectivityReceiverSideGrid(const ovk_connectivity_r *Receivers,
  const ovk_grid **ReceiverGrid) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(ReceiverGrid, "Invalid receiver grid pointer.");

  *ReceiverGrid = Receivers->grid;

}

static bool EditingPoints(const ovk_connectivity_r *Receivers) {

  return Receivers->points_edit_ref_count > 0;

}

static bool EditingSources(const ovk_connectivity_r *Receivers) {

  return Receivers->sources_edit_ref_count > 0;

}

static bool EditingSourceRanks(const ovk_connectivity_r *Receivers) {

  return Receivers->source_ranks_edit_ref_count > 0;

}

void PRIVATE(GetConnectivityReceiverSideEdits)(const ovk_connectivity_r *Receivers,
  const t_connectivity_r_edits **Edits) {

  *Edits = &Receivers->edits;

}

void PRIVATE(ResetConnectivityReceiverSideEdits)(ovk_connectivity_r *Receivers) {

  MPI_Barrier(Receivers->comm);

  DefaultEdits(&Receivers->edits);

  MPI_Barrier(Receivers->comm);

}

static void DefaultEdits(t_connectivity_r_edits *Edits) {

  Edits->count = false;
  Edits->points = false;
  Edits->sources = false;

}
