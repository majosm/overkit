// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ConnectivityR.h"

#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/Logger.h"

static bool EditingProperties(const ovk_connectivity_r *Receivers);
static bool EditingPoints(const ovk_connectivity_r *Receivers);
static bool EditingSources(const ovk_connectivity_r *Receivers);
static bool EditingSourceRanks(const ovk_connectivity_r *Receivers);

static void DefaultProperties(ovk_connectivity_r_properties *Properties);
static void DefaultEdits(t_connectivity_r_edits *Edits);

void PRIVATE(CreateConnectivityReceiverSide)(ovk_connectivity_r **Receivers_, const ovk_grid *Grid,
  int SourceGridID, t_logger *Logger, t_error_handler *ErrorHandler) {

  int iDim;

  const ovk_grid_properties *GridProperties;
  ovkGetGridProperties(Grid, &GridProperties);

  int GridID;
  ovkGetGridPropertyID(GridProperties, &GridID);

  int NumDims;
  ovkGetGridPropertyDimension(GridProperties, &NumDims);

  MPI_Comm Comm;
  ovkGetGridPropertyComm(GridProperties, &Comm);

  MPI_Barrier(Comm);

  int CommSize, CommRank;
  MPI_Comm_size(Comm, &CommSize);
  MPI_Comm_rank(Comm, &CommRank);

  *Receivers_ = malloc(sizeof(ovk_connectivity_r));
  ovk_connectivity_r *Receivers = *Receivers_;

  DefaultProperties(&Receivers->properties);

  Receivers->properties.grid_id = GridID;
  Receivers->properties.source_grid_id = SourceGridID;
  Receivers->properties.num_dims = NumDims;
  Receivers->properties.comm = Comm;
  Receivers->properties.comm_size = CommSize;
  Receivers->properties.comm_rank = CommRank;

  Receivers->properties_edit_ref_count = 0;

  Receivers->logger = Logger;
  Receivers->error_handler = ErrorHandler;

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

  MPI_Barrier(Receivers->properties.comm);

}

void PRIVATE(DestroyConnectivityReceiverSide)(ovk_connectivity_r **Receivers_) {

  ovk_connectivity_r *Receivers = *Receivers_;

  MPI_Barrier(Receivers->properties.comm);

  free(Receivers->points[0]);
  free(Receivers->sources[0]);
  free(Receivers->source_ranks);

  MPI_Comm Comm = Receivers->properties.comm;

  free_null(Receivers_);

  MPI_Barrier(Comm);

}

void ovkGetConnectivityReceiverSideProperties(const ovk_connectivity_r *Receivers,
  const ovk_connectivity_r_properties **Properties) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");

  *Properties = &Receivers->properties;

}

void ovkResizeReceivers(ovk_connectivity_r *Receivers, size_t NumReceivers) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(NumReceivers >= 0, "Invalid receiver count.");

  MPI_Barrier(Receivers->properties.comm);

  OVK_DEBUG_ASSERT(!EditingProperties(Receivers), "Cannot resize receivers while editing "
    "properties.");
  OVK_DEBUG_ASSERT(!EditingPoints(Receivers), "Cannot resize receivers while editing points.");
  OVK_DEBUG_ASSERT(!EditingSources(Receivers), "Cannot resize receivers while editing sources.");
  OVK_DEBUG_ASSERT(!EditingSourceRanks(Receivers), "Cannot resize receivers while editing source "
    "ranks.");

  int iDim;
  size_t iReceiver;

  free(Receivers->points[0]);
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers->points[iDim] = NULL;
  }
  free(Receivers->sources[0]);
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers->sources[iDim] = NULL;
  }
  free_null(&Receivers->source_ranks);

  Receivers->properties.count = NumReceivers;

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

  MPI_Barrier(Receivers->properties.comm);

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
  OVK_DEBUG_ASSERT(!EditingProperties(Receivers), "Cannot edit points while editing properties.");

  bool StartEdit = Receivers->points_edit_ref_count == 0;
  ++Receivers->points_edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Receivers->properties.comm);
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
    MPI_Barrier(Receivers->properties.comm);
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
  OVK_DEBUG_ASSERT(!EditingProperties(Receivers), "Cannot edit sources while editing properties.");

  bool StartEdit = Receivers->sources_edit_ref_count == 0;
  ++Receivers->sources_edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Receivers->properties.comm);
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
    MPI_Barrier(Receivers->properties.comm);
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
  OVK_DEBUG_ASSERT(!EditingProperties(Receivers), "Cannot edit source ranks while editing "
    "properties.");

  bool StartEdit = Receivers->source_ranks_edit_ref_count == 0;
  ++Receivers->source_ranks_edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Receivers->properties.comm);
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
    MPI_Barrier(Receivers->properties.comm);
  }

}

void ovkGetConnectivityReceiverSideGrid(const ovk_connectivity_r *Receivers,
  const ovk_grid **ReceiverGrid) {

  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(ReceiverGrid, "Invalid receiver grid pointer.");

  *ReceiverGrid = Receivers->grid;

}

static bool EditingProperties(const ovk_connectivity_r *Receivers) {

  return Receivers->properties_edit_ref_count > 0;

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

  MPI_Barrier(Receivers->properties.comm);

  DefaultEdits(&Receivers->edits);

  MPI_Barrier(Receivers->properties.comm);

}

static void DefaultProperties(ovk_connectivity_r_properties *Properties) {

  Properties->grid_id = -1;
  Properties->source_grid_id = -1;
  Properties->num_dims = 2;
  Properties->comm = MPI_COMM_NULL;
  Properties->comm_size = 0;
  Properties->comm_rank = 0;
  Properties->count = 0;

}

void ovkGetConnectivityReceiverSidePropertyGridID(const ovk_connectivity_r_properties *Properties,
  int *GridID) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  *GridID = Properties->grid_id;

}

void ovkGetConnectivityReceiverSidePropertySourceGridID(const ovk_connectivity_r_properties
  *Properties, int *SourceGridID) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(SourceGridID, "Invalid source grid ID pointer.");

  *SourceGridID = Properties->source_grid_id;

}

void ovkGetConnectivityReceiverSidePropertyDimension(const ovk_connectivity_r_properties *Properties,
  int *NumDims) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Properties->num_dims;

}

void ovkGetConnectivityReceiverSidePropertyComm(const ovk_connectivity_r_properties *Properties,
  MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Properties->comm;

}

void ovkGetConnectivityReceiverSidePropertyCommSize(const ovk_connectivity_r_properties *Properties,
  int *CommSize) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  *CommSize = Properties->comm_size;

}

void ovkGetConnectivityReceiverSidePropertyCommRank(const ovk_connectivity_r_properties *Properties,
  int *CommRank) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  *CommRank = Properties->comm_rank;

}

void ovkGetConnectivityReceiverSidePropertyCount(const ovk_connectivity_r_properties *Properties,
  size_t *NumReceivers) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(NumReceivers, "Invalid num receivers pointer.");

  *NumReceivers = Properties->count;

}

static void DefaultEdits(t_connectivity_r_edits *Edits) {

  Edits->count = false;
  Edits->points = false;
  Edits->sources = false;

}
