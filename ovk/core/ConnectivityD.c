// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ConnectivityD.h"

#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/Logger.h"
#include "ovk/core/Range.h"

static bool EditingProperties(const ovk_connectivity_d *Donors);
static bool EditingExtents(const ovk_connectivity_d *Donors);
static bool EditingCoords(const ovk_connectivity_d *Donors);
static bool EditingInterpCoefs(const ovk_connectivity_d *Donors);
static bool EditingDestinations(const ovk_connectivity_d *Donors);
static bool EditingDestinationRanks(const ovk_connectivity_d *Donors);

static void DefaultProperties(ovk_connectivity_d_properties *Properties);
static void DefaultEdits(t_connectivity_d_edits *Edits);

void PRIVATE(CreateConnectivityDonorSide)(ovk_connectivity_d **Donors_, const ovk_grid *Grid,
  t_logger *Logger, t_error_handler *ErrorHandler) {

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

  *Donors_ = malloc(sizeof(ovk_connectivity_d));
  ovk_connectivity_d *Donors = *Donors_;

  DefaultProperties(&Donors->properties);

  Donors->properties.grid_id = GridID;
  Donors->properties.num_dims = NumDims;
  Donors->properties.comm = Comm;
  Donors->properties.comm_size = CommSize;
  Donors->properties.comm_rank = CommRank;

  Donors->properties_edit_ref_count = 0;

  Donors->logger = Logger;
  Donors->error_handler = ErrorHandler;

  Donors->grid = Grid;

  DefaultEdits(&Donors->edits);

  Donors->extents = malloc(2*sizeof(int **));
  Donors->extents[0] = malloc(MAX_DIMS*sizeof(int *));
  Donors->extents[1] = malloc(MAX_DIMS*sizeof(int *));
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors->extents[0][iDim] = NULL;
    Donors->extents[1][iDim] = NULL;
  }
  Donors->extents_edit_ref_count = 0;

  Donors->coords = malloc(NumDims*sizeof(double *));
  for (iDim = 0; iDim < NumDims; ++iDim) {
    Donors->coords[iDim] = NULL;
  }
  Donors->coords_edit_ref_count = 0;

  Donors->interp_coefs = malloc(NumDims*sizeof(double **));
  for (iDim = 0; iDim < NumDims; ++iDim) {
    Donors->interp_coefs[iDim] = NULL;
  }
  Donors->interp_coefs_edit_ref_count = 0;

  Donors->destinations = malloc(MAX_DIMS*sizeof(int *));
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors->destinations[iDim] = NULL;
  }
  Donors->destinations_edit_ref_count = 0;

  Donors->dest_ranks = NULL;
  Donors->dest_ranks_edit_ref_count = 0;

  MPI_Barrier(Donors->properties.comm);

}

void PRIVATE(DestroyConnectivityDonorSide)(ovk_connectivity_d **Donors_) {

  int iDim, iCoef;

  ovk_connectivity_d *Donors = *Donors_;

  MPI_Barrier(Donors->properties.comm);

  int NumDims = Donors->properties.num_dims;
  int MaxDonorSize = Donors->properties.max_donor_size;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    free(Donors->extents[0][iDim]);
    free(Donors->extents[1][iDim]);
  }
  free(Donors->extents[0]);
  free(Donors->extents[1]);
  free(Donors->extents);

  for (iDim = 0; iDim < NumDims; ++iDim) {
    free(Donors->coords[iDim]);
  }
  free(Donors->coords);

  for (iDim = 0; iDim < NumDims; ++iDim) {
    for (iCoef = 0; iCoef < MaxDonorSize; ++iCoef) {
      free(Donors->interp_coefs[iDim][iCoef]);
    }
    free(Donors->interp_coefs[iDim]);
  }
  free(Donors->interp_coefs);

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    free(Donors->destinations[iDim]);
  }
  free(Donors->destinations);

  free(Donors->dest_ranks);

  MPI_Comm Comm = Donors->properties.comm;

  free_null(Donors_);

  MPI_Barrier(Comm);

}

void ovkGetConnectivityDonorSideProperties(const ovk_connectivity_d *Donors,
  const ovk_connectivity_d_properties **Properties) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");

  *Properties = &Donors->properties;

}

void ovkResizeDonors(ovk_connectivity_d *Donors, size_t NumDonors, int MaxDonorSize) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(NumDonors >= 0, "Invalid donor count.");
  OVK_DEBUG_ASSERT(MaxDonorSize >= 0, "Invalid max donor size.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(!EditingProperties(Donors), "Cannot resize donors while editing properties.");
  OVK_DEBUG_ASSERT(!EditingExtents(Donors), "Cannot resize donors while editing extents.");
  OVK_DEBUG_ASSERT(!EditingCoords(Donors), "Cannot resize donors while editing coords.");
  OVK_DEBUG_ASSERT(!EditingInterpCoefs(Donors), "Cannot resize donors while editing interp coefs.");
  OVK_DEBUG_ASSERT(!EditingDestinations(Donors), "Cannot resize donors while editing destinations.");
  OVK_DEBUG_ASSERT(!EditingDestinationRanks(Donors), "Cannot resize donors while editing "
    "destination ranks.");

  int iDim, iCoef;
  size_t iDonor;

  int NumDims = Donors->properties.num_dims;

  int PrevMaxDonorSize = Donors->properties.max_donor_size;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    free_null(&Donors->extents[0][iDim]);
    free_null(&Donors->extents[1][iDim]);
  }
  for (iDim = 0; iDim < NumDims; ++iDim) {
    free_null(&Donors->coords[iDim]);
  }
  for (iDim = 0; iDim < NumDims; ++iDim) {
    for (iCoef = 0; iCoef < PrevMaxDonorSize; ++iCoef) {
      free_null(&Donors->interp_coefs[iDim][iCoef]);
    }
    if (MaxDonorSize != PrevMaxDonorSize) {
      free_null(&Donors->interp_coefs[iDim]);
      Donors->interp_coefs[iDim] = malloc(MaxDonorSize*sizeof(double *));
      for (iCoef = 0; iCoef < MaxDonorSize; ++iCoef) {
        Donors->interp_coefs[iDim][iCoef] = NULL;
      }
    }
  }
  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    free_null(&Donors->destinations[iDim]);
  }
  free_null(&Donors->dest_ranks);

  Donors->properties.num_donors = NumDonors;
  Donors->properties.max_donor_size = MaxDonorSize;

  if (NumDonors > 0) {

    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Donors->extents[0][iDim] = malloc(NumDonors*sizeof(int));
      Donors->extents[1][iDim] = malloc(NumDonors*sizeof(int));
    }
    for (iDim = 0; iDim < NumDims; ++iDim) {
      Donors->coords[iDim] = malloc(NumDonors*sizeof(double));
    }
    for (iDim = 0; iDim < NumDims; ++iDim) {
      for (iCoef = 0; iCoef < MaxDonorSize; ++iCoef) {
        Donors->interp_coefs[iDim][iCoef] = malloc(NumDonors*sizeof(double));
      }
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Donors->destinations[iDim] = malloc(NumDonors*sizeof(int));
    }
    Donors->dest_ranks = malloc(NumDonors*sizeof(int));

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      ovk_range EmptyRange;
      ovkDefaultRange(&EmptyRange, NumDims);
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Donors->extents[0][iDim][iDonor] = EmptyRange.b[iDim];
        Donors->extents[1][iDim][iDonor] = EmptyRange.e[iDim];
      }
      for (iDim = 0; iDim < NumDims; ++iDim) {
        Donors->coords[iDim][iDonor] = 0.;
      }
      for (iDim = 0; iDim < NumDims; ++iDim) {
        for (iCoef = 0; iCoef < MaxDonorSize; ++iCoef) {
          Donors->interp_coefs[iDim][iCoef][iDonor] = 0.;
        }
      }
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Donors->destinations[iDim][iDonor] = 0;
      }
      Donors->dest_ranks[iDonor] = -1;
    }

  }

  Donors->edits.num_donors = true;
  Donors->edits.extents = true;
  Donors->edits.coords = true;
  Donors->edits.interp_coefs = true;
  Donors->edits.destinations = true;

  MPI_Barrier(Donors->properties.comm);

}

void ovkEditDonorExtents(ovk_connectivity_d *Donors, int Dimension, int **Begins, int **Ends) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Begins, "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(Ends, "Invalid ends pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(!EditingProperties(Donors), "Cannot edit extents while editing properties.");

  ++Donors->extents_edit_ref_count;

  *Begins = Donors->extents[0][Dimension];
  *Ends = Donors->extents[1][Dimension];

  MPI_Barrier(Donors->properties.comm);

}

void ovkReleaseDonorExtents(ovk_connectivity_d *Donors, int Dimension, int **Begins, int **Ends) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Begins, "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(Ends, "Invalid ends pointer.");
  OVK_DEBUG_ASSERT(*Begins == Donors->extents[0][Dimension], "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(*Ends == Donors->extents[1][Dimension], "Invalid ends pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(EditingExtents(Donors), "Unable to release extents; not currently being edited.");

  --Donors->extents_edit_ref_count;
  bool EndEdit = Donors->extents_edit_ref_count == 0;

  *Begins = NULL;
  *Ends = NULL;

  if (EndEdit) {
    Donors->edits.extents = true;
  }

  MPI_Barrier(Donors->properties.comm);

}

void ovkEditDonorCoords(ovk_connectivity_d *Donors, int Dimension, double **Coords) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < Donors->properties.num_dims,
    "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(!EditingProperties(Donors), "Cannot edit coords while editing properties.");

  ++Donors->coords_edit_ref_count;

  *Coords = Donors->coords[Dimension];

  MPI_Barrier(Donors->properties.comm);

}

void ovkReleaseDonorCoords(ovk_connectivity_d *Donors, int Dimension, double **Coords) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < Donors->properties.num_dims,
    "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");
  OVK_DEBUG_ASSERT(*Coords == Donors->coords[Dimension], "Invalid coords pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(EditingCoords(Donors), "Unable to release coords; not currently being edited.");

  --Donors->coords_edit_ref_count;
  bool EndEdit = Donors->coords_edit_ref_count == 0;

  *Coords = NULL;

  if (EndEdit) {
    Donors->edits.coords = true;
  }

  MPI_Barrier(Donors->properties.comm);

}

void ovkEditDonorInterpCoefs(ovk_connectivity_d *Donors, int Dimension, int Index,
  double **InterpCoefs) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < Donors->properties.num_dims,
    "Invalid dimension.");
  OVK_DEBUG_ASSERT(Index >= 0 && Index < Donors->properties.max_donor_size, "Invalid index.");
  OVK_DEBUG_ASSERT(InterpCoefs, "Invalid interp coefs pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(!EditingProperties(Donors), "Cannot edit interp coefs while editing properties.");

  ++Donors->interp_coefs_edit_ref_count;

  *InterpCoefs = Donors->interp_coefs[Dimension][Index];

  MPI_Barrier(Donors->properties.comm);

}

void ovkReleaseDonorInterpCoefs(ovk_connectivity_d *Donors, int Dimension, int Index,
  double **InterpCoefs) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < Donors->properties.num_dims,
    "Invalid dimension.");
  OVK_DEBUG_ASSERT(Index >= 0 && Index < Donors->properties.max_donor_size, "Invalid index.");
  OVK_DEBUG_ASSERT(InterpCoefs, "Invalid donor coefs pointer.");
  OVK_DEBUG_ASSERT(*InterpCoefs == Donors->interp_coefs[Dimension][Index], "Invalid interp coefs "
    "pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(EditingInterpCoefs(Donors), "Unable to release interp coefs; not currently "
    "being edited.");

  --Donors->interp_coefs_edit_ref_count;
  bool EndEdit = Donors->interp_coefs_edit_ref_count == 0;

  *InterpCoefs = NULL;

  if (EndEdit) {
    Donors->edits.interp_coefs = true;
  }

  MPI_Barrier(Donors->properties.comm);

}

void ovkEditDonorDestinations(ovk_connectivity_d *Donors, int Dimension, int **Destinations) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(!EditingProperties(Donors), "Cannot edit destinations while editing properties.");

  ++Donors->destinations_edit_ref_count;

  *Destinations = Donors->destinations[Dimension];

  MPI_Barrier(Donors->properties.comm);

}

void ovkReleaseDonorDestinations(ovk_connectivity_d *Donors, int Dimension, int **Destinations) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");
  OVK_DEBUG_ASSERT(*Destinations == Donors->destinations[Dimension], "Invalid destinations pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(EditingDestinations(Donors), "Unable to release destinations; not currently "
    "being edited.");

  --Donors->destinations_edit_ref_count;
  bool EndEdit = Donors->destinations_edit_ref_count == 0;

  *Destinations = NULL;

  if (EndEdit) {

    Donors->edits.destinations = true;

  }

  MPI_Barrier(Donors->properties.comm);

}

void ovkEditDonorDestinationRanks(ovk_connectivity_d *Donors, int **DestinationRanks) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(!EditingProperties(Donors), "Cannot edit destination ranks while editing "
    "properties.");

  ++Donors->dest_ranks_edit_ref_count;

  *DestinationRanks = Donors->dest_ranks;

  MPI_Barrier(Donors->properties.comm);

}

void ovkReleaseDonorDestinationRanks(ovk_connectivity_d *Donors, int **DestinationRanks) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");
  OVK_DEBUG_ASSERT(*DestinationRanks == Donors->dest_ranks, "Invalid destination ranks pointer.");

  MPI_Barrier(Donors->properties.comm);

  OVK_DEBUG_ASSERT(EditingDestinationRanks(Donors), "Unable to release destination ranks; not "
    "currently being edited.");

  --Donors->dest_ranks_edit_ref_count;
  bool EndEdit = Donors->dest_ranks_edit_ref_count == 0;

  *DestinationRanks = NULL;

  if (EndEdit) {
    Donors->edits.destinations = true;
  }

  MPI_Barrier(Donors->properties.comm);

}

void ovkGetConnectivityDonorSideGrid(const ovk_connectivity_d *Donors, const ovk_grid **DonorGrid) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(DonorGrid, "Invalid donor grid pointer.");

  *DonorGrid = Donors->grid;

}

static bool EditingProperties(const ovk_connectivity_d *Donors) {

  return Donors->properties_edit_ref_count > 0;

}

static bool EditingExtents(const ovk_connectivity_d *Donors) {

  return Donors->extents_edit_ref_count > 0;

}

static bool EditingCoords(const ovk_connectivity_d *Donors) {

  return Donors->coords_edit_ref_count > 0;

}

static bool EditingInterpCoefs(const ovk_connectivity_d *Donors) {

  return Donors->interp_coefs_edit_ref_count > 0;

}

static bool EditingDestinations(const ovk_connectivity_d *Donors) {

  return Donors->destinations_edit_ref_count > 0;

}

static bool EditingDestinationRanks(const ovk_connectivity_d *Donors) {

  return Donors->dest_ranks_edit_ref_count > 0;

}

void PRIVATE(GetConnectivityDonorSideEdits)(const ovk_connectivity_d *Donors,
  const t_connectivity_d_edits **Edits) {

  *Edits = &Donors->edits;

}

void PRIVATE(ResetConnectivityDonorSideEdits)(ovk_connectivity_d *Donors) {

  MPI_Barrier(Donors->properties.comm);

  DefaultEdits(&Donors->edits);

  MPI_Barrier(Donors->properties.comm);

}

static void DefaultProperties(ovk_connectivity_d_properties *Properties) {

  Properties->grid_id = -1;
  Properties->num_dims = 2;
  Properties->comm = MPI_COMM_NULL;
  Properties->comm_size = 0;
  Properties->comm_rank = 0;
  Properties->num_donors = 0;
  Properties->max_donor_size = 0;

}

void ovkGetConnectivityDonorSidePropertyGridID(const ovk_connectivity_d_properties *Properties,
  int *GridID) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  *GridID = Properties->grid_id;

}

void ovkGetConnectivityDonorSidePropertyDimension(const ovk_connectivity_d_properties *Properties,
  int *NumDims) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Properties->num_dims;

}

void ovkGetConnectivityDonorSidePropertyComm(const ovk_connectivity_d_properties *Properties,
  MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Properties->comm;

}

void ovkGetConnectivityDonorSidePropertyCommSize(const ovk_connectivity_d_properties *Properties,
  int *CommSize) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  *CommSize = Properties->comm_size;

}

void ovkGetConnectivityDonorSidePropertyCommRank(const ovk_connectivity_d_properties *Properties,
  int *CommRank) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  *CommRank = Properties->comm_rank;

}

void ovkGetConnectivityDonorSidePropertyDonorCount(const ovk_connectivity_d_properties *Properties,
  size_t *NumDonors) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(NumDonors, "Invalid num donors pointer.");

  *NumDonors = Properties->num_donors;

}

void ovkGetConnectivityDonorSidePropertyMaxDonorSize(const ovk_connectivity_d_properties *Properties,
  int *MaxDonorSize) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(MaxDonorSize, "Invalid max donor size pointer.");

  *MaxDonorSize = Properties->max_donor_size;

}

static void DefaultEdits(t_connectivity_d_edits *Edits) {

  Edits->num_donors = false;
  Edits->extents = false;
  Edits->coords = false;
  Edits->interp_coefs = false;
  Edits->destinations = false;

}
