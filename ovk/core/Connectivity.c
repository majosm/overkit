// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Connectivity.h"

#include "ovk/core/Cart.h"
#include "ovk/core/ConnectivityD.h"
#include "ovk/core/ConnectivityR.h"
#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/Logger.h"
#include "ovk/core/MPIUtils.h"
#include "ovk/core/Range.h"
#include "ovk/core/TextUtils.h"

static bool EditingDonorSide(const ovk_connectivity *Connectivity);
static bool EditingReceiverSide(const ovk_connectivity *Connectivity);

static void CreateDonorSideGlobal(t_connectivity_donor_side_container **DonorsContainer,
  const ovk_grid *DonorGrid, int DestinationGridID, MPI_Comm Comm, int CommRank, t_logger *Logger,
  t_error_handler *ErrorHandler);
static void DestroyDonorSideGlobal(t_connectivity_donor_side_container **DonorsContainer);

static void CreateReceiverSideGlobal(t_connectivity_receiver_side_container **ReceiversContainer,
  const ovk_grid *ReceiverGrid, int SourceGridID, MPI_Comm Comm, int CommRank, t_logger *Logger,
  t_error_handler *ErrorHandler);
static void DestroyReceiverSideGlobal(t_connectivity_receiver_side_container **ReceiversContainer);

static void CreateDonorSideContainer(t_connectivity_donor_side_container **Container,
  const ovk_grid *Grid, ovk_connectivity_d *Donors, MPI_Comm Comm, int CommSize);
static void DestroyDonorSideContainer(t_connectivity_donor_side_container **Container);

static void CreateReceiverSideContainer(t_connectivity_receiver_side_container **Container,
  const ovk_grid *Grid, ovk_connectivity_r *Receivers, MPI_Comm Comm, int CommSize);
static void DestroyReceiverSideContainer(t_connectivity_receiver_side_container **Container);

static void EditDonorSideGlobal(ovk_connectivity *Connectivity, ovk_connectivity_d **Donors);
static void ReleaseDonorSideGlobal(ovk_connectivity *Connectivity, ovk_connectivity_d **Donors);
static void EditReceiverSideGlobal(ovk_connectivity *Connectivity, ovk_connectivity_r **Receivers);
static void ReleaseReceiverSideGlobal(ovk_connectivity *Connectivity, ovk_connectivity_r **Receivers);

static void DefaultEdits(t_connectivity_edits *Edits);

void PRIVATE(CreateConnectivity)(ovk_connectivity **Connectivity_, int NumDims, MPI_Comm Comm_,
  const ovk_grid *DonorGrid, const ovk_grid *ReceiverGrid, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  MPI_Comm Comm;
  MPI_Comm_dup(Comm_, &Comm);

  MPI_Barrier(Comm);

  int CommSize, CommRank;
  MPI_Comm_size(Comm, &CommSize);
  MPI_Comm_rank(Comm, &CommRank);

  *Connectivity_ = malloc(sizeof(ovk_connectivity));
  ovk_connectivity *Connectivity = *Connectivity_;

  Connectivity->logger = Logger;
  Connectivity->error_handler = ErrorHandler;

  int DonorGridID;
  bool IsDonorGridRoot = false;
  if (DonorGrid) {
    int GridCommRank;
    ovkGetGridCommRank(DonorGrid, &GridCommRank);
    IsDonorGridRoot = GridCommRank == 0;
    if (IsDonorGridRoot) {
      ovkGetGridID(DonorGrid, &DonorGridID);
    }
  }
  BroadcastAnySource(&DonorGridID, 1, MPI_INT, IsDonorGridRoot, Comm);

  int ReceiverGridID;
  bool IsReceiverGridRoot = false;
  if (ReceiverGrid) {
    int GridCommRank;
    ovkGetGridCommRank(ReceiverGrid, &GridCommRank);
    IsReceiverGridRoot = GridCommRank == 0;
    if (IsReceiverGridRoot) {
      ovkGetGridID(ReceiverGrid, &ReceiverGridID);
    }
  }
  BroadcastAnySource(&ReceiverGridID, 1, MPI_INT, IsReceiverGridRoot, Comm);

  t_connectivity_donor_side_container *DonorsContainer;
  CreateDonorSideGlobal(&DonorsContainer, DonorGrid, ReceiverGridID, Comm, CommRank, Logger,
    ErrorHandler);

  t_connectivity_receiver_side_container *ReceiversContainer;
  CreateReceiverSideGlobal(&ReceiversContainer, ReceiverGrid, DonorGridID, Comm, CommRank, Logger,
    ErrorHandler);

  Connectivity->donor_grid_id = DonorGridID;
  Connectivity->receiver_grid_id = ReceiverGridID;

  sprintf(Connectivity->name, "(%s,%s)", DonorsContainer->grid_info->name,
    ReceiversContainer->grid_info->name);

  Connectivity->num_dims = NumDims;

  Connectivity->comm = Comm;
  Connectivity->comm_size = CommSize;
  Connectivity->comm_rank = CommRank;

  DefaultEdits(&Connectivity->edits);

  Connectivity->donors_container = DonorsContainer;
  Connectivity->receivers_container = ReceiversContainer;

  MPI_Barrier(Connectivity->comm);

  LogStatus(Connectivity->logger, Connectivity->comm_rank == 0, 0, "Created connectivity %s.",
    Connectivity->name);

}

void PRIVATE(DestroyConnectivity)(ovk_connectivity **Connectivity_) {

  ovk_connectivity *Connectivity = *Connectivity_;

  MPI_Barrier(Connectivity->comm);

  DestroyDonorSideGlobal(&Connectivity->donors_container);
  DestroyReceiverSideGlobal(&Connectivity->receivers_container);

  t_logger *Logger = Connectivity->logger;
  MPI_Comm Comm = Connectivity->comm;
  bool IsRoot = Connectivity->comm_rank == 0;
  char Name[OVK_NAME_LENGTH];
  strncpy(Name, Connectivity->name, OVK_NAME_LENGTH);

  free_null(Connectivity_);

  MPI_Barrier(Comm);

  MPI_Comm_free(&Comm);

  LogStatus(Logger, IsRoot, 0, "Destroyed connectivity %s.", Name);

}

void PRIVATE(CreateConnectivityInfo)(ovk_connectivity_info **Info_,
  const ovk_connectivity *Connectivity, MPI_Comm Comm, int CommRank) {

  *Info_ = malloc(sizeof(ovk_connectivity_info));
  ovk_connectivity_info *Info = *Info_;

  bool IsLocal = Connectivity != NULL;
  bool IsRoot = false;
  if (IsLocal) {
    IsRoot = Connectivity->comm_rank == 0;
  }

  int RootRank;
  if (IsRoot) RootRank = CommRank;
  BroadcastAnySource(&RootRank, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) {
    Info->donor_grid_id = Connectivity->donor_grid_id;
    Info->receiver_grid_id = Connectivity->receiver_grid_id;
    strcpy(Info->name, Connectivity->name);
    Info->num_dims = Connectivity->num_dims;
  }
  MPI_Bcast(&Info->donor_grid_id, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info->receiver_grid_id, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info->name, OVK_NAME_LENGTH, MPI_CHAR, RootRank, Comm);
  MPI_Bcast(&Info->num_dims, 1, MPI_INT, RootRank, Comm);

  Info->root_rank = RootRank;

}

void PRIVATE(DestroyConnectivityInfo)(ovk_connectivity_info **Info) {

  free_null(Info);

}

static void CreateDonorSideGlobal(t_connectivity_donor_side_container **DonorsContainer,
  const ovk_grid *DonorGrid, int DestinationGridID, MPI_Comm Comm, int CommRank, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  ovk_connectivity_d *Donors = NULL;
  if (DonorGrid) {
    CreateConnectivityDonorSide(&Donors, DonorGrid, DestinationGridID, Logger, ErrorHandler);
  }

  CreateDonorSideContainer(DonorsContainer, DonorGrid, Donors, Comm, CommRank);

}

static void DestroyDonorSideGlobal(t_connectivity_donor_side_container **DonorsContainer_) {

  t_connectivity_donor_side_container *DonorsContainer = *DonorsContainer_;

  if (DonorsContainer->has_local_data) {
    DestroyConnectivityDonorSide(&DonorsContainer->donors);
  }

  DestroyDonorSideContainer(DonorsContainer_);

}

static void CreateReceiverSideGlobal(t_connectivity_receiver_side_container **ReceiversContainer,
  const ovk_grid *ReceiverGrid, int SourceGridID, MPI_Comm Comm, int CommRank, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  ovk_connectivity_r *Receivers = NULL;
  if (ReceiverGrid) {
    CreateConnectivityReceiverSide(&Receivers, ReceiverGrid, SourceGridID, Logger, ErrorHandler);
  }

  CreateReceiverSideContainer(ReceiversContainer, ReceiverGrid, Receivers, Comm, CommRank);

}

static void DestroyReceiverSideGlobal(t_connectivity_receiver_side_container **ReceiversContainer_)
  {

  t_connectivity_receiver_side_container *ReceiversContainer = *ReceiversContainer_;

  if (ReceiversContainer->has_local_data) {
    DestroyConnectivityReceiverSide(&ReceiversContainer->receivers);
  }

  DestroyReceiverSideContainer(ReceiversContainer_);

}

void ovkGetConnectivityDonorGridID(const ovk_connectivity *Connectivity, int *DonorGridID) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(DonorGridID, "Invalid donor grid ID pointer.");

  *DonorGridID = Connectivity->donor_grid_id;

}

void ovkGetConnectivityReceiverGridID(const ovk_connectivity *Connectivity, int *ReceiverGridID) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(ReceiverGridID, "Invalid receiver grid ID pointer.");

  *ReceiverGridID = Connectivity->receiver_grid_id;

}

void ovkGetConnectivityName(const ovk_connectivity *Connectivity, char *Name) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  strcpy(Name, Connectivity->name);

}

void ovkGetConnectivityDimension(const ovk_connectivity *Connectivity, int *NumDims) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Connectivity->num_dims;

}

void ovkGetConnectivityComm(const ovk_connectivity *Connectivity, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Connectivity->comm;

}

void ovkGetConnectivityCommSize(const ovk_connectivity *Connectivity, int *CommSize) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  *CommSize = Connectivity->comm_size;

}

void ovkGetConnectivityCommRank(const ovk_connectivity *Connectivity, int *CommRank) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  *CommRank = Connectivity->comm_rank;

}

void ovkGetConnectivityDonorGridInfo(const ovk_connectivity *Connectivity,
  const ovk_grid_info **DonorGridInfo) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(DonorGridInfo, "Invalid donor grid info pointer.");

  const t_connectivity_donor_side_container *Container = Connectivity->donors_container;

  *DonorGridInfo = Container->grid_info;

}

void ovkGetConnectivityReceiverGridInfo(const ovk_connectivity *Connectivity,
  const ovk_grid_info **ReceiverGridInfo) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(ReceiverGridInfo, "Invalid receiver grid info pointer.");

  const t_connectivity_receiver_side_container *Container = Connectivity->receivers_container;

  *ReceiverGridInfo = Container->grid_info;

}

bool ovkRankHasConnectivityDonorSide(const ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  const t_connectivity_donor_side_container *Container = Connectivity->donors_container;

  return Container->has_local_data;

}

void ovkGetConnectivityDonorSide(const ovk_connectivity *Connectivity,
  const ovk_connectivity_d **Donors) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");

  const t_connectivity_donor_side_container *Container = Connectivity->donors_container;
  OVK_DEBUG_ASSERT(Container->has_local_data, "Connectivity %s does not have donor-side data on "
    "rank @rank@.", Connectivity->name);

  *Donors = Container->donors;

}

void ovkEditConnectivityDonorSideLocal(ovk_connectivity *Connectivity, ovk_connectivity_d **Donors)
  {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");

  EditDonorSideGlobal(Connectivity, Donors);

}

void ovkEditConnectivityDonorSideRemote(ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  EditDonorSideGlobal(Connectivity, NULL);
  
}

static void EditDonorSideGlobal(ovk_connectivity *Connectivity, ovk_connectivity_d **Donors) {

  bool IsLocal = Donors != NULL;

  t_connectivity_donor_side_container *Container = Connectivity->donors_container;

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Container->has_local_data, "Connectivity %s does not have donor-side data on "
      "rank @rank@.", Connectivity->name);
  }

  bool StartEdit = Container->edit_ref_count == 0;
  ++Container->edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Connectivity->comm);
  }

  if (IsLocal) {
    *Donors = Container->donors;
  }

}

void ovkReleaseConnectivityDonorSideLocal(ovk_connectivity *Connectivity, ovk_connectivity_d **Donors)
  {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(*Donors, "Invalid donors pointer.");

  ReleaseDonorSideGlobal(Connectivity, Donors);

}

void ovkReleaseConnectivityDonorSideRemote(ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  ReleaseDonorSideGlobal(Connectivity, NULL);
  
}

static void ReleaseDonorSideGlobal(ovk_connectivity *Connectivity, ovk_connectivity_d **Donors_) {

  OVK_DEBUG_ASSERT(EditingDonorSide(Connectivity), "Unable to release connectivity %s donor-side "
    "data; not currently being edited.", Connectivity->name);

  bool IsLocal = Donors_ != NULL;

  t_connectivity_donor_side_container *Container = Connectivity->donors_container;

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Container->has_local_data, "Connectivity %s does not have donor-side data on "
      "rank @rank@.", Connectivity->name);
    OVK_DEBUG_ASSERT(*Donors_ == Container->donors, "Invalid donors pointer.");
  }

  --Container->edit_ref_count;
  bool EndEdit = Container->edit_ref_count == 0;

  if (IsLocal) {
    *Donors_ = NULL;
  }

  if (EndEdit) {

    MPI_Barrier(Connectivity->comm);

    const t_connectivity_d_edits *DonorsEdits;
    if (IsLocal) {
      ovk_connectivity_d *Donors = Container->donors;
      GetConnectivityDonorSideEdits(Donors, &DonorsEdits);
    }

    t_connectivity_edits *Edits = &Connectivity->edits;

    int RootRank = Container->grid_info->root_rank;
    bool IsGridRoot = Connectivity->comm_rank == RootRank;

    int EditedNumDonors = 0;
    if (IsGridRoot) EditedNumDonors = DonorsEdits->count;
    MPI_Bcast(&EditedNumDonors, 1, MPI_INT, RootRank, Connectivity->comm);
    Edits->num_donors = Edits->num_donors || EditedNumDonors;

    int EditedExtents = 0;
    if (IsGridRoot) EditedExtents = DonorsEdits->extents;
    MPI_Bcast(&EditedExtents, 1, MPI_INT, RootRank, Connectivity->comm);
    Edits->donor_extents = Edits->donor_extents || EditedExtents;

    int EditedCoords = 0;
    if (IsGridRoot) EditedCoords = DonorsEdits->coords;
    MPI_Bcast(&EditedCoords, 1, MPI_INT, RootRank, Connectivity->comm);
    Edits->donor_coords = Edits->donor_coords || EditedCoords;

    int EditedInterpCoefs = 0;
    if (IsGridRoot) EditedInterpCoefs = DonorsEdits->interp_coefs;
    MPI_Bcast(&EditedInterpCoefs, 1, MPI_INT, RootRank, Connectivity->comm);
    Edits->donor_interp_coefs = Edits->donor_interp_coefs || EditedInterpCoefs;

    int EditedDestinations = 0;
    if (IsGridRoot) EditedDestinations = DonorsEdits->destinations;
    MPI_Bcast(&EditedDestinations, 1, MPI_INT, RootRank, Connectivity->comm);
    Edits->donor_destinations = Edits->donor_destinations || EditedDestinations;

  }

}

bool ovkRankHasConnectivityReceiverSide(const ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  const t_connectivity_receiver_side_container *Container = Connectivity->receivers_container;

  return Container->has_local_data;

}

void ovkGetConnectivityReceiverSide(const ovk_connectivity *Connectivity,
  const ovk_connectivity_r **Receivers) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");

  const t_connectivity_receiver_side_container *Container = Connectivity->receivers_container;
  OVK_DEBUG_ASSERT(Container->has_local_data, "Connectivity %s does not have receiver-side data on "
    "rank @rank@.", Connectivity->name);

  *Receivers = Container->receivers;

}

void ovkEditConnectivityReceiverSideLocal(ovk_connectivity *Connectivity,
  ovk_connectivity_r **Receivers) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");

  EditReceiverSideGlobal(Connectivity, Receivers);

}

void ovkEditConnectivityReceiverSideRemote(ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  EditReceiverSideGlobal(Connectivity, NULL);
  
}

static void EditReceiverSideGlobal(ovk_connectivity *Connectivity, ovk_connectivity_r **Receivers) {

  bool IsLocal = Receivers != NULL;

  t_connectivity_receiver_side_container *Container = Connectivity->receivers_container;

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Container->has_local_data, "Connectivity %s does not have receiver-side data "
      "on rank @rank@.", Connectivity->name);
  }

  bool StartEdit = Container->edit_ref_count == 0;
  ++Container->edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Connectivity->comm);
  }

  if (IsLocal) {
    *Receivers = Container->receivers;
  }

}

void ovkReleaseConnectivityReceiverSideLocal(ovk_connectivity *Connectivity,
  ovk_connectivity_r **Receivers) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(Receivers, "Invalid receivers pointer.");
  OVK_DEBUG_ASSERT(*Receivers, "Invalid receivers pointer.");

  ReleaseReceiverSideGlobal(Connectivity, Receivers);

}

void ovkReleaseConnectivityReceiverSideRemote(ovk_connectivity *Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  ReleaseReceiverSideGlobal(Connectivity, NULL);
  
}

static void ReleaseReceiverSideGlobal(ovk_connectivity *Connectivity,
  ovk_connectivity_r **Receivers_) {

  OVK_DEBUG_ASSERT(EditingReceiverSide(Connectivity), "Unable to release connectivity %s "
    "receiver-side data; not currently being edited.", Connectivity->name);

  bool IsLocal = Receivers_ != NULL;

  t_connectivity_receiver_side_container *Container = Connectivity->receivers_container;

  if (IsLocal) {
    OVK_DEBUG_ASSERT(Container->has_local_data, "Connectivity %s does not have receiver-side data "
      "on rank @rank@.", Connectivity->name);
    OVK_DEBUG_ASSERT(*Receivers_ == Container->receivers, "Invalid receivers pointer.");
  }

  --Container->edit_ref_count;
  bool EndEdit = Container->edit_ref_count == 0;

  if (IsLocal) {
    *Receivers_ = NULL;
  }

  if (EndEdit) {

    MPI_Barrier(Connectivity->comm);

    const t_connectivity_r_edits *ReceiversEdits;
    if (IsLocal) {
      ovk_connectivity_r *Receivers = Container->receivers;
      GetConnectivityReceiverSideEdits(Receivers, &ReceiversEdits);
    }

    t_connectivity_edits *Edits = &Connectivity->edits;

    int RootRank = Container->grid_info->root_rank;
    bool IsGridRoot = Connectivity->comm_rank == RootRank;

    int EditedNumReceivers = 0;
    if (IsGridRoot) EditedNumReceivers = ReceiversEdits->count;
    MPI_Bcast(&EditedNumReceivers, 1, MPI_INT, RootRank, Connectivity->comm);
    Edits->num_receivers = Edits->num_receivers || EditedNumReceivers;

    int EditedPoints = 0;
    if (IsGridRoot) EditedPoints = ReceiversEdits->points;
    MPI_Bcast(&EditedPoints, 1, MPI_INT, RootRank, Connectivity->comm);
    Edits->receiver_points = Edits->receiver_points || EditedPoints;

    int EditedSources = 0;
    if (IsGridRoot) EditedSources = ReceiversEdits->sources;
    MPI_Bcast(&EditedSources, 1, MPI_INT, RootRank, Connectivity->comm);
    Edits->receiver_sources = Edits->receiver_sources || EditedSources;

  }

}

static bool EditingDonorSide(const ovk_connectivity *Connectivity) {

  const t_connectivity_donor_side_container *Container = Connectivity->donors_container;

  return Container->edit_ref_count > 0;

}

static bool EditingReceiverSide(const ovk_connectivity *Connectivity) {

  const t_connectivity_receiver_side_container *Container = Connectivity->receivers_container;

  return Container->edit_ref_count > 0;

}

void PRIVATE(GetConnectivityEdits)(const ovk_connectivity *Connectivity,
  const t_connectivity_edits **Edits) {

  *Edits = &Connectivity->edits;

}

void PRIVATE(ResetConnectivityEdits)(ovk_connectivity *Connectivity) {

  DefaultEdits(&Connectivity->edits);

  t_connectivity_donor_side_container *DonorsContainer = Connectivity->donors_container;
  if (DonorsContainer->has_local_data) {
    ResetConnectivityDonorSideEdits(DonorsContainer->donors);
  }

  t_connectivity_receiver_side_container *ReceiversContainer = Connectivity->receivers_container;
  if (ReceiversContainer->has_local_data) {
    ResetConnectivityReceiverSideEdits(ReceiversContainer->receivers);
  }

}

static void CreateDonorSideContainer(t_connectivity_donor_side_container **Container_,
  const ovk_grid *Grid, ovk_connectivity_d *Donors, MPI_Comm Comm, int CommRank) {

  bool IsLocal = Grid != NULL;

  *Container_ = malloc(sizeof(t_connectivity_donor_side_container));
  t_connectivity_donor_side_container *Container = *Container_;

  CreateGridInfo(&Container->grid_info, Grid, Comm, CommRank);
  Container->edit_ref_count = 0;

  Container->has_local_data = IsLocal;

  if (IsLocal) {
    Container->grid = Grid;
    Container->donors = Donors;
  }

}

static void DestroyDonorSideContainer(t_connectivity_donor_side_container **Container_) {

  t_connectivity_donor_side_container *Container = *Container_;

  DestroyGridInfo(&Container->grid_info);

  free_null(Container_);

}

static void CreateReceiverSideContainer(t_connectivity_receiver_side_container **Container_,
  const ovk_grid *Grid, ovk_connectivity_r *Receivers, MPI_Comm Comm, int CommRank) {

  bool IsLocal = Grid != NULL;

  *Container_ = malloc(sizeof(t_connectivity_receiver_side_container));
  t_connectivity_receiver_side_container *Container = *Container_;

  CreateGridInfo(&Container->grid_info, Grid, Comm, CommRank);
  Container->edit_ref_count = 0;

  Container->has_local_data = IsLocal;

  if (IsLocal) {
    Container->grid = Grid;
    Container->receivers = Receivers;
  }

}

static void DestroyReceiverSideContainer(t_connectivity_receiver_side_container **Container_) {

  t_connectivity_receiver_side_container *Container = *Container_;

  DestroyGridInfo(&Container->grid_info);

  free_null(Container_);

}

static void DefaultEdits(t_connectivity_edits *Edits) {

  Edits->num_donors = false;
  Edits->donor_extents = false;
  Edits->donor_coords = false;
  Edits->donor_interp_coefs = false;
  Edits->donor_destinations = false;
  Edits->num_receivers = false;
  Edits->receiver_points = false;
  Edits->receiver_sources = false;

}

void ovkGetConnectivityInfoDonorGridID(const ovk_connectivity_info *Info, int *DonorGridID) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(DonorGridID, "Invalid donor grid ID pointer.");

  *DonorGridID = Info->donor_grid_id;

}

void ovkGetConnectivityInfoReceiverGridID(const ovk_connectivity_info *Info, int *ReceiverGridID) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(ReceiverGridID, "Invalid receiver grid ID pointer.");

  *ReceiverGridID = Info->receiver_grid_id;

}

void ovkGetConnectivityInfoName(const ovk_connectivity_info *Info, char *Name) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  strcpy(Name, Info->name);

}

void ovkGetConnectivityInfoDimension(const ovk_connectivity_info *Info, int *NumDims) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Info->num_dims;

}

void ovkGetConnectivityInfoRootRank(const ovk_connectivity_info *Info, int *RootRank) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(RootRank, "Invalid root rank pointer.");

  *RootRank = Info->root_rank;

}
