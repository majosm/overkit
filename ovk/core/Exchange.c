// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Exchange.h"

#include "ovk/core/Connectivity.h"
#include "ovk/core/Domain.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/MiscUtils.h"
#include "ovk/core/MPIUtils.h"
#include "ovk/core/PartitionHash.h"
#include "ovk/core/Range.h"
#include "ovk/core/TextUtils.h"

typedef struct {
  const ovk_exchange *exchange;
  void **buffers;
  MPI_Request *mpi_requests;
} t_send_request_data;

typedef struct {
  const ovk_exchange *exchange;
  ovk_data_type data_type;
  int count;
  void **receiver_data;
  void **buffers;
  MPI_Request *mpi_requests;
} t_recv_request_data;

static void ResizeDonors(ovk_exchange *Exchange);
static void ResizeReceivers(ovk_exchange *Exchange);

static void UpdateCollectSendInfo(ovk_exchange *Exchange);
static void UpdateCollectReceiveInfo(ovk_exchange *Exchange);

static void UpdateSourceRanks(ovk_exchange *Exchange);
static void UpdateDestRanks(ovk_exchange *Exchange);

static void UpdateDonorsSorted(ovk_exchange *Exchange);
static void UpdateReceiversSorted(ovk_exchange *Exchange);

static void UpdateSendInfo(ovk_exchange *Exchange);
static void UpdateReceiveInfo(ovk_exchange *Exchange);

static void CompleteReceive(t_recv_request_data *RequestData);

static void CreateSendRequest(const ovk_exchange *Exchange, int DataSize, int Count,
  ovk_request **Request);
static void DestroySendRequest(ovk_request **Request);
static void CreateReceiveRequest(const ovk_exchange *Exchange, ovk_data_type DataType, int DataSize,
  int Count, void **ReceiverData, ovk_request **Request);
static void DestroyReceiveRequest(ovk_request **Request);

void PRIVATE(CreateExchange)(ovk_exchange **Exchange_, const ovk_connectivity *Connectivity,
  t_logger *Logger, t_error_handler *ErrorHandler) {

  *Exchange_ = malloc(sizeof(ovk_exchange));
  ovk_exchange *Exchange = *Exchange_;

  Exchange->connectivity = Connectivity;

  const ovk_connectivity_properties *ConnectivityProperties;
  ovkGetConnectivityProperties(Connectivity, &ConnectivityProperties);

  ovkGetConnectivityPropertyDimension(ConnectivityProperties, &Exchange->num_dims);
  ovkGetConnectivityPropertyComm(ConnectivityProperties, &Exchange->comm);
  ovkGetConnectivityPropertyCommSize(ConnectivityProperties, &Exchange->comm_size);
  ovkGetConnectivityPropertyCommRank(ConnectivityProperties, &Exchange->comm_rank);

  MPI_Barrier(Exchange->comm);

  const ovk_grid_info *DonorGridInfo;
  ovkGetConnectivityDonorGridInfo(Connectivity, &DonorGridInfo);

  const ovk_grid_info *ReceiverGridInfo;
  ovkGetConnectivityDonorGridInfo(Connectivity, &ReceiverGridInfo);

  const ovk_connectivity_d *Donors = NULL;
  const ovk_grid *DonorGrid = NULL;
  if (ovkRankHasConnectivityDonorSide(Connectivity)) {
    ovkGetConnectivityDonorSide(Connectivity, &Donors);
    ovkGetConnectivityDonorSideGrid(Donors, &DonorGrid);
  }

  const ovk_connectivity_r *Receivers = NULL;
  const ovk_grid *ReceiverGrid = NULL;
  if (ovkRankHasConnectivityReceiverSide(Connectivity)) {
    ovkGetConnectivityReceiverSide(Connectivity, &Receivers);
    ovkGetConnectivityReceiverSideGrid(Receivers, &ReceiverGrid);
  }

  Exchange->logger = Logger;
  Exchange->error_handler = ErrorHandler;

  Exchange->num_collect_sends = 0;
  Exchange->collect_send_dest_ranks = NULL;
  Exchange->num_collect_send_points = NULL;
  Exchange->collect_send_points = NULL;

  Exchange->num_collect_recvs = 0;
  Exchange->collect_recv_source_ranks = NULL;
  Exchange->num_collect_recv_points = NULL;
  Exchange->collect_recv_points = NULL;

  Exchange->num_remote_donor_points = NULL;
  Exchange->remote_donor_points = NULL;
  Exchange->remote_donor_point_collect_recv_indices = NULL;
  Exchange->remote_donor_point_collect_recv_buffer_offsets = NULL;

  Exchange->donors_sorted = NULL;
  Exchange->receivers_sorted = NULL;

  Exchange->donor_dest_ranks = NULL;
  Exchange->receiver_source_ranks = NULL;

  Exchange->num_sends = 0;
  Exchange->send_ranks = NULL;
  Exchange->send_counts = NULL;
  Exchange->donor_send_indices = NULL;

  Exchange->num_recvs = 0;
  Exchange->recv_ranks = NULL;
  Exchange->recv_counts = NULL;
  Exchange->receiver_recv_indices = NULL;

  ovk_range DonorGridGlobalRange;
  ovkGetGridInfoGlobalRange(DonorGridInfo, &DonorGridGlobalRange);

  ovk_range DonorGridLocalRange;
  if (DonorGrid) {
    const ovk_grid_properties *DonorGridProperties;
    ovkGetGridProperties(DonorGrid, &DonorGridProperties);
    ovkGetGridPropertyLocalRange(DonorGridProperties, &DonorGridLocalRange);
  } else {
    ovkDefaultRange(&DonorGridLocalRange, Exchange->num_dims);
  }

  CreatePartitionHash(&Exchange->source_hash, Exchange->num_dims, Exchange->comm,
    &DonorGridGlobalRange, &DonorGridLocalRange);

  ovk_range ReceiverGridGlobalRange;
  ovkGetGridInfoGlobalRange(ReceiverGridInfo, &ReceiverGridGlobalRange);

  ovk_range ReceiverGridLocalRange;
  if (ReceiverGrid) {
    const ovk_grid_properties *ReceiverGridProperties;
    ovkGetGridProperties(ReceiverGrid, &ReceiverGridProperties);
    ovkGetGridPropertyLocalRange(ReceiverGridProperties, &ReceiverGridLocalRange);
  } else {
    ovkDefaultRange(&ReceiverGridLocalRange, Exchange->num_dims);
  }

  CreatePartitionHash(&Exchange->destination_hash, Exchange->num_dims, Exchange->comm,
    &ReceiverGridGlobalRange, &ReceiverGridLocalRange);

  MPI_Barrier(Exchange->comm);

  LogStatus(Exchange->logger, Exchange->comm_rank == 0, 0, "Created exchange %s.",
    Connectivity->properties.name);

}

void PRIVATE(DestroyExchange)(ovk_exchange **Exchange_) {

  int iDim, iCollectSend, iCollectRecv;

  ovk_exchange *Exchange = *Exchange_;

  MPI_Barrier(Exchange->comm);

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  DestroyPartitionHash(&Exchange->source_hash);
  DestroyPartitionHash(&Exchange->destination_hash);

  free(Exchange->collect_send_dest_ranks);
  free(Exchange->num_collect_send_points);

  for (iCollectSend = 0; iCollectSend < Exchange->num_collect_sends; ++iCollectSend) {
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      free(Exchange->collect_send_points[iCollectSend][iDim]);
    }
    free(Exchange->collect_send_points[iCollectSend]);
  }
  free(Exchange->collect_send_points);

  free(Exchange->collect_recv_source_ranks);
  free(Exchange->num_collect_recv_points);

  for (iCollectRecv = 0; iCollectRecv < Exchange->num_collect_recvs; ++iCollectRecv) {
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      free(Exchange->collect_recv_points[iCollectRecv][iDim]);
    }
    free(Exchange->collect_recv_points[iCollectRecv]);
  }
  free(Exchange->collect_recv_points);

  free(Exchange->num_remote_donor_points);

  if (Exchange->remote_donor_points) {
    free(Exchange->remote_donor_points[0]);
  }
  free(Exchange->remote_donor_points);

  if (Exchange->remote_donor_point_collect_recv_indices) {
    free(Exchange->remote_donor_point_collect_recv_indices[0]);
  }
  free(Exchange->remote_donor_point_collect_recv_indices);

  if (Exchange->remote_donor_point_collect_recv_buffer_offsets) {
    free(Exchange->remote_donor_point_collect_recv_buffer_offsets[0]);
  }
  free(Exchange->remote_donor_point_collect_recv_buffer_offsets);

  free(Exchange->donors_sorted);
  free(Exchange->receivers_sorted);

  free(Exchange->donor_dest_ranks);
  free(Exchange->receiver_source_ranks);

  free(Exchange->send_ranks);
  free(Exchange->send_counts);
  free(Exchange->donor_send_indices);

  free(Exchange->recv_ranks);
  free(Exchange->recv_counts);
  free(Exchange->receiver_recv_indices);

  t_logger *Logger = Exchange->logger;
  MPI_Comm Comm = Exchange->comm;
  bool IsRoot = Exchange->comm_rank == 0;
  char Name[OVK_NAME_LENGTH];
  strncpy(Name, Connectivity->properties.name, OVK_NAME_LENGTH);

  free_null(Exchange_);

  MPI_Barrier(Comm);

  LogStatus(Logger, IsRoot, 0, "Destroyed exchange %s.", Name);

}

void PRIVATE(CreateExchangeInfo)(ovk_exchange_info **Info_, const ovk_exchange *Exchange,
  MPI_Comm Comm, int CommRank) {

  *Info_ = malloc(sizeof(ovk_exchange_info));
  ovk_exchange_info *Info = *Info_;

  bool IsLocal = Exchange != NULL;
  bool IsRoot = false;
  if (IsLocal) {
    IsRoot = Exchange->comm_rank == 0;
  }

  const ovk_connectivity *Connectivity = NULL;
  if (IsLocal) {
    Connectivity = Exchange->connectivity;
  }

  int RootRank;
  if (IsRoot) RootRank = CommRank;
  BroadcastAnySource(&RootRank, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) {
    Info->donor_grid_id = Connectivity->properties.donor_grid_id;
    Info->receiver_grid_id = Connectivity->properties.receiver_grid_id;
    strcpy(Info->name, Connectivity->properties.name);
    Info->num_dims = Connectivity->properties.num_dims;
  }
  MPI_Bcast(&Info->donor_grid_id, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info->receiver_grid_id, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info->name, OVK_NAME_LENGTH, MPI_CHAR, RootRank, Comm);
  MPI_Bcast(&Info->num_dims, 1, MPI_INT, RootRank, Comm);

  Info->root_rank = RootRank;

}

void PRIVATE(DestroyExchangeInfo)(ovk_exchange_info **Info) {

  free_null(Info);

}

bool ovkRankHasExchangeDonorSide(const ovk_exchange *Exchange) {

  OVK_DEBUG_ASSERT(Exchange, "Invalid exchange pointer.");

  return ovkRankHasConnectivityDonorSide(Exchange->connectivity);

}

bool ovkRankHasExchangeReceiverSide(const ovk_exchange *Exchange) {

  OVK_DEBUG_ASSERT(Exchange, "Invalid exchange pointer.");

  return ovkRankHasConnectivityReceiverSide(Exchange->connectivity);

}

void PRIVATE(UpdateExchange)(ovk_exchange *Exchange) {

  MPI_Barrier(Exchange->comm);

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  LogStatus(Exchange->logger, Connectivity->properties.comm_rank == 0, 0,
    "Updating exchange %s...", Connectivity->properties.name);

  const t_connectivity_edits *Edits;
  GetConnectivityEdits(Connectivity, &Edits);

  bool NeedToResizeDonors = false;
  bool NeedToResizeReceivers = false;
  bool NeedToUpdateCollectInfo = false;
  bool NeedToUpdateDonorsSorted = false;
  bool NeedToUpdateReceiversSorted = false;
  bool NeedToUpdateSourceDestRanks = false;
  bool NeedToUpdateSendRecvInfo = false;

  if (Edits->num_donors) {
    NeedToResizeDonors = true;
    NeedToUpdateCollectInfo = true;
    NeedToUpdateDonorsSorted = true;
    NeedToUpdateSourceDestRanks = true;
    NeedToUpdateSendRecvInfo = true;
  }

  if (Edits->donor_extents) {
    NeedToUpdateCollectInfo = true;
  }

  if (Edits->donor_destinations) {
    NeedToUpdateDonorsSorted = true;
    NeedToUpdateSourceDestRanks = true;
    NeedToUpdateSendRecvInfo = true;
  }

  if (Edits->num_receivers) {
    NeedToResizeReceivers = true;
    NeedToUpdateReceiversSorted = true;
    NeedToUpdateSourceDestRanks = true;
    NeedToUpdateSendRecvInfo = true;
  }

  if (Edits->receiver_sources) {
    NeedToUpdateReceiversSorted = true;
    NeedToUpdateSourceDestRanks = true;
    NeedToUpdateSendRecvInfo = true;
  }

  if (NeedToResizeDonors) ResizeDonors(Exchange);
  if (NeedToResizeReceivers) ResizeReceivers(Exchange);

  if (NeedToUpdateCollectInfo) {
    UpdateCollectSendInfo(Exchange);
    UpdateCollectReceiveInfo(Exchange);
  }

  if (NeedToUpdateDonorsSorted) UpdateDonorsSorted(Exchange);
  if (NeedToUpdateReceiversSorted) UpdateReceiversSorted(Exchange);

  if (NeedToUpdateSourceDestRanks) {
    UpdateSourceRanks(Exchange);
    UpdateDestRanks(Exchange);
  }

  if (NeedToUpdateSendRecvInfo) {
    UpdateSendInfo(Exchange);
    UpdateReceiveInfo(Exchange);
  }

  MPI_Barrier(Exchange->comm);

  LogStatus(Exchange->logger, Connectivity->properties.comm_rank == 0, 0, 
    "Done updating exchange %s.", Connectivity->properties.name);

}

static void ResizeDonors(ovk_exchange *Exchange) {

  size_t iDonor;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  free_null(&Exchange->donors_sorted);
  free_null(&Exchange->donor_dest_ranks);
  free_null(&Exchange->donor_send_indices);

  free_null(&Exchange->num_remote_donor_points);
  if (Exchange->remote_donor_points) {
    free(Exchange->remote_donor_points[0]);
  }
  free_null(&Exchange->remote_donor_points);
  if (Exchange->remote_donor_point_collect_recv_indices) {
    free(Exchange->remote_donor_point_collect_recv_indices[0]);
  }
  free_null(&Exchange->remote_donor_point_collect_recv_indices);
  if (Exchange->remote_donor_point_collect_recv_buffer_offsets) {
    free(Exchange->remote_donor_point_collect_recv_buffer_offsets[0]);
  }
  free_null(&Exchange->remote_donor_point_collect_recv_buffer_offsets);

  size_t NumDonors = 0;
  if (ovkRankHasConnectivityDonorSide(Connectivity)) {
    const ovk_connectivity_d *Donors;
    ovkGetConnectivityDonorSide(Connectivity, &Donors);
    const ovk_connectivity_d_properties *DonorsProperties;
    ovkGetConnectivityDonorSideProperties(Donors, &DonorsProperties);
    ovkGetConnectivityDonorSidePropertyCount(DonorsProperties, &NumDonors);
  }

  if (NumDonors > 0) {

    Exchange->donors_sorted = malloc(NumDonors*sizeof(size_t));
    Exchange->donor_dest_ranks = malloc(NumDonors*sizeof(int));
    Exchange->donor_send_indices = malloc(NumDonors*sizeof(int));
    Exchange->num_remote_donor_points = malloc(NumDonors*sizeof(int));
    Exchange->remote_donor_points = malloc(NumDonors*sizeof(size_t *));
    Exchange->remote_donor_points[0] = NULL;
    Exchange->remote_donor_point_collect_recv_indices = malloc(NumDonors*sizeof(int *));
    Exchange->remote_donor_point_collect_recv_indices[0] = NULL;
    Exchange->remote_donor_point_collect_recv_buffer_offsets = malloc(NumDonors*sizeof(size_t *));
    Exchange->remote_donor_point_collect_recv_buffer_offsets[0] = NULL;

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      Exchange->donors_sorted[iDonor] = -1;
      Exchange->donor_dest_ranks[iDonor] = -1;
      Exchange->donor_send_indices[iDonor] = -1;
      Exchange->num_remote_donor_points[iDonor] = 0;
      Exchange->remote_donor_points[iDonor] = NULL;
      Exchange->remote_donor_point_collect_recv_indices[iDonor] = NULL;
      Exchange->remote_donor_point_collect_recv_buffer_offsets[iDonor] = NULL;
    }
    
  }

}

static void ResizeReceivers(ovk_exchange *Exchange) {

  size_t iReceiver;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  free_null(&Exchange->receivers_sorted);
  free_null(&Exchange->receiver_source_ranks);
  free_null(&Exchange->receiver_recv_indices);

  size_t NumReceivers = 0;
  if (ovkRankHasConnectivityReceiverSide(Connectivity)) {
    const ovk_connectivity_r *Receivers;
    ovkGetConnectivityReceiverSide(Connectivity, &Receivers);
    const ovk_connectivity_r_properties *ReceiversProperties;
    ovkGetConnectivityReceiverSideProperties(Receivers, &ReceiversProperties);
    ovkGetConnectivityReceiverSidePropertyCount(ReceiversProperties, &NumReceivers);
  }

  if (NumReceivers > 0) {

    Exchange->receivers_sorted = malloc(NumReceivers*sizeof(size_t));
    Exchange->receiver_source_ranks = malloc(NumReceivers*sizeof(int));
    Exchange->receiver_recv_indices = malloc(NumReceivers*sizeof(int));

    for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      Exchange->receivers_sorted[iReceiver] = -1;
      Exchange->receiver_source_ranks[iReceiver] = -1;
      Exchange->receiver_recv_indices[iReceiver] = -1;
    }

  }

}

static void UpdateCollectSendInfo(ovk_exchange *Exchange) {

  int iDim, iCollectSend;

  free_null(&Exchange->collect_send_dest_ranks);
  free_null(&Exchange->num_collect_send_points);

  for (iCollectSend = 0; iCollectSend < Exchange->num_collect_sends; ++iCollectSend) {
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      free_null(&Exchange->collect_send_points[iCollectSend][iDim]);
    }
    free_null(&Exchange->collect_send_points[iCollectSend]);
  }
  free_null(&Exchange->collect_send_points);

  Exchange->num_collect_sends = 0;

  int NumDims = Exchange->num_dims;
  const ovk_connectivity *Connectivity = Exchange->connectivity;

  const ovk_connectivity_d *Donors;
  size_t NumDonors = 0;
  int MaxSize = 0;
  if (ovkRankHasConnectivityDonorSide(Connectivity)) {
    ovkGetConnectivityDonorSide(Connectivity, &Donors);
    const ovk_connectivity_d_properties *DonorsProperties;
    ovkGetConnectivityDonorSideProperties(Donors, &DonorsProperties);
    ovkGetConnectivityDonorSidePropertyCount(DonorsProperties, &NumDonors);
    ovkGetConnectivityDonorSidePropertyMaxSize(DonorsProperties, &MaxSize);
  }

  if (NumDonors > 0) {

    size_t iDonor;
    int iNeighbor;
    size_t iPoint, iCollectSendPoint;
    int i, j, k;

    const ovk_grid *Grid;
    ovkGetConnectivityDonorSideGrid(Donors, &Grid);

    const ovk_grid_properties *GridProperties;
    ovkGetGridProperties(Grid, &GridProperties);

    ovk_range GlobalRange, LocalRange;
    ovkGetGridPropertyGlobalRange(GridProperties, &GlobalRange);
    ovkGetGridPropertyLocalRange(GridProperties, &LocalRange);

    int NumNeighbors;
    ovkGetGridPropertyNumNeighbors(GridProperties, &NumNeighbors);

    ovk_cart Cart;
    ovkGetGridCart(Grid, &Cart);

    ovk_range *SendToNeighborDataRanges = malloc(NumNeighbors*sizeof(ovk_range));

    for (iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      ovkDefaultRange(&SendToNeighborDataRanges[iNeighbor], NumDims);
    }

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS] = {
        Donors->extents[0][0][iDonor],
        Donors->extents[0][1][iDonor],
        Donors->extents[0][2][iDonor]
      };
      int DonorEnd[MAX_DIMS] = {
        Donors->extents[1][0][iDonor],
        Donors->extents[1][1][iDonor],
        Donors->extents[1][2][iDonor]
      };
      ovk_range DonorRange;
      ovkSetRange(&DonorRange, NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = ovkRangeIncludes(&GlobalRange, &DonorRange);
      for (iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
        ovk_range NeighborRange;
        GetGridNeighborRange(Grid, iNeighbor, &NeighborRange);
        bool Overlaps;
        if (AwayFromEdge) {
          Overlaps = ovkRangeOverlaps(&NeighborRange, &DonorRange);
        } else {
          Overlaps = false;
          for (k = DonorRange.b[2]; k < DonorRange.e[2]; ++k) {
            for (j = DonorRange.b[1]; j < DonorRange.e[1]; ++j) {
              for (i = DonorRange.b[0]; i < DonorRange.e[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                int AdjustedPoint[MAX_DIMS];
                ovkCartPeriodicAdjust(&Cart, Point, AdjustedPoint);
                if (ovkRangeContains(&NeighborRange, AdjustedPoint)) {
                  Overlaps = true;
                  goto endloop1;
                }
              }
            }
          }
          endloop1: ;
        }
        if (Overlaps) {
          ovk_range IntersectRange;
          ovkRangeIntersect(&LocalRange, &DonorRange, &IntersectRange);
          ovkRangeUnion(&SendToNeighborDataRanges[iNeighbor], &IntersectRange,
            &SendToNeighborDataRanges[iNeighbor]);
        }
      }
    }

    int NumCollectSends = 0;
    for (iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!ovkRangeIsEmpty(&SendToNeighborDataRanges[iNeighbor])) {
        ++NumCollectSends;
      }
    }

    int *CollectSendIndexToNeighbor = malloc(NumCollectSends*sizeof(int));
    iCollectSend = 0;
    for (iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!ovkRangeIsEmpty(&SendToNeighborDataRanges[iNeighbor])) {
        CollectSendIndexToNeighbor[iCollectSend] = iNeighbor;
        ++iCollectSend;
      }
    }

    Exchange->num_collect_sends = NumCollectSends;

    Exchange->collect_send_dest_ranks = malloc(NumCollectSends*sizeof(int));

    for (iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
      int NeighborRank;
      GetGridNeighborRank(Grid, iNeighbor, &NeighborRank);
      Exchange->collect_send_dest_ranks[iCollectSend] = NeighborRank;
    }

    bool **CollectSendMasks = malloc(NumCollectSends*sizeof(bool *));
    for (iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
      size_t NumPoints;
      ovkRangeCount(&SendToNeighborDataRanges[iNeighbor], &NumPoints);
      CollectSendMasks[iCollectSend] = malloc(NumPoints*sizeof(bool));
      for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
        CollectSendMasks[iCollectSend][iPoint] = false;
      }
    }

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS] = {
        Donors->extents[0][0][iDonor],
        Donors->extents[0][1][iDonor],
        Donors->extents[0][2][iDonor]
      };
      int DonorEnd[MAX_DIMS] = {
        Donors->extents[1][0][iDonor],
        Donors->extents[1][1][iDonor],
        Donors->extents[1][2][iDonor]
      };
      ovk_range DonorRange;
      ovkSetRange(&DonorRange, NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = ovkRangeIncludes(&GlobalRange, &DonorRange);
      for (iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
        iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
        ovk_range NeighborRange;
        GetGridNeighborRange(Grid, iNeighbor, &NeighborRange);
        bool Overlaps;
        if (AwayFromEdge) {
          Overlaps = ovkRangeOverlaps(&NeighborRange, &DonorRange);
        } else {
          Overlaps = false;
          for (k = DonorRange.b[2]; k < DonorRange.e[2]; ++k) {
            for (j = DonorRange.b[1]; j < DonorRange.e[1]; ++j) {
              for (i = DonorRange.b[0]; i < DonorRange.e[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                int AdjustedPoint[MAX_DIMS];
                ovkCartPeriodicAdjust(&Cart, Point, AdjustedPoint);
                if (ovkRangeContains(&NeighborRange, AdjustedPoint)) {
                  Overlaps = true;
                  goto endloop2;
                }
              }
            }
          }
          endloop2: ;
        }
        if (Overlaps) {
          ovk_range IntersectRange;
          ovkRangeIntersect(&LocalRange, &DonorRange, &IntersectRange);
          for (k = IntersectRange.b[2]; k < IntersectRange.e[2]; ++k) {
            for (j = IntersectRange.b[1]; j < IntersectRange.e[1]; ++j) {
              for (i = IntersectRange.b[0]; i < IntersectRange.e[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                ovkRangeTupleToIndex(&SendToNeighborDataRanges[iNeighbor], OVK_COLUMN_MAJOR, Point,
                  &iPoint);
                CollectSendMasks[iCollectSend][iPoint] = true;
              }
            }
          }
        }
      }
    }

    Exchange->num_collect_send_points = malloc(NumCollectSends*sizeof(size_t));

    for (iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
      size_t NumPoints;
      ovkRangeCount(&SendToNeighborDataRanges[iNeighbor], &NumPoints);
      size_t NumCollectSendPoints = 0;
      for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
        if (CollectSendMasks[iCollectSend][iPoint]) {
          ++NumCollectSendPoints;
        }
      }
      Exchange->num_collect_send_points[iCollectSend] = NumCollectSendPoints;
    }

    Exchange->collect_send_points = malloc(NumCollectSends*sizeof(int **));
    for (iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      Exchange->collect_send_points[iCollectSend] = malloc(MAX_DIMS*sizeof(int *));
      size_t NumCollectSendPoints = Exchange->num_collect_send_points[iCollectSend];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Exchange->collect_send_points[iCollectSend][iDim] = malloc(NumCollectSendPoints*
          sizeof(int));
      }
      iNeighbor = CollectSendIndexToNeighbor[iCollectSend];
      size_t NumPoints;
      ovkRangeCount(&SendToNeighborDataRanges[iNeighbor], &NumPoints);
      iCollectSendPoint = 0;
      for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
        if (CollectSendMasks[iCollectSend][iPoint]) {
          int Point[MAX_DIMS];
          ovkRangeIndexToTuple(&SendToNeighborDataRanges[iNeighbor], OVK_COLUMN_MAJOR, iPoint, Point);
          for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
            Exchange->collect_send_points[iCollectSend][iDim][iCollectSendPoint] = Point[iDim];
          }
          ++iCollectSendPoint;
        }
      }
    }

    free(SendToNeighborDataRanges);
    free(CollectSendIndexToNeighbor);
    for (iCollectSend = 0; iCollectSend < NumCollectSends; ++iCollectSend) {
      free(CollectSendMasks[iCollectSend]);
    }
    free(CollectSendMasks);

  }

}

static void UpdateCollectReceiveInfo(ovk_exchange *Exchange) {

  int iDim, iCollectRecv;
  int iDonorPoint;

  free_null(&Exchange->collect_recv_source_ranks);
  free_null(&Exchange->num_collect_recv_points);

  for (iCollectRecv = 0; iCollectRecv < Exchange->num_collect_recvs; ++iCollectRecv) {
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      free_null(&Exchange->collect_recv_points[iCollectRecv][iDim]);
    }
    free_null(&Exchange->collect_recv_points[iCollectRecv]);
  }
  free_null(&Exchange->collect_recv_points);

  Exchange->num_collect_recvs = 0;

  if (Exchange->remote_donor_points) {
    free_null(&Exchange->remote_donor_points[0]);
  }

  if (Exchange->remote_donor_point_collect_recv_indices) {
    free_null(&Exchange->remote_donor_point_collect_recv_indices[0]);
  }

  if (Exchange->remote_donor_point_collect_recv_buffer_offsets) {
    free_null(&Exchange->remote_donor_point_collect_recv_buffer_offsets[0]);
  }

  int NumDims = Exchange->num_dims;
  const ovk_connectivity *Connectivity = Exchange->connectivity;

  const ovk_connectivity_d *Donors;
  size_t NumDonors = 0;
  int MaxSize = 0;
  if (ovkRankHasConnectivityDonorSide(Connectivity)) {
    ovkGetConnectivityDonorSide(Connectivity, &Donors);
    const ovk_connectivity_d_properties *DonorsProperties;
    ovkGetConnectivityDonorSideProperties(Donors, &DonorsProperties);
    ovkGetConnectivityDonorSidePropertyCount(DonorsProperties, &NumDonors);
    ovkGetConnectivityDonorSidePropertyMaxSize(DonorsProperties, &MaxSize);
  }

  if (NumDonors > 0) {

    size_t iDonor;
    int iNeighbor;
    size_t iPoint, iCollectRecvPoint;
    int i, j, k;

    const ovk_grid *Grid;
    ovkGetConnectivityDonorSideGrid(Donors, &Grid);

    const ovk_grid_properties *GridProperties;
    ovkGetGridProperties(Grid, &GridProperties);

    ovk_range GlobalRange, LocalRange;
    ovkGetGridPropertyGlobalRange(GridProperties, &GlobalRange);
    ovkGetGridPropertyLocalRange(GridProperties, &LocalRange);

    int NumNeighbors;
    ovkGetGridPropertyNumNeighbors(GridProperties, &NumNeighbors);

    ovk_cart Cart;
    ovkGetGridCart(Grid, &Cart);

    ovk_range *RecvFromNeighborDataRanges = malloc(NumNeighbors*sizeof(ovk_range));

    for (iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      ovkDefaultRange(&RecvFromNeighborDataRanges[iNeighbor], NumDims);
    }

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS] = {
        Donors->extents[0][0][iDonor],
        Donors->extents[0][1][iDonor],
        Donors->extents[0][2][iDonor]
      };
      int DonorEnd[MAX_DIMS] = {
        Donors->extents[1][0][iDonor],
        Donors->extents[1][1][iDonor],
        Donors->extents[1][2][iDonor]
      };
      ovk_range DonorRange;
      ovkSetRange(&DonorRange, NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = ovkRangeIncludes(&GlobalRange, &DonorRange);
      for (iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
        ovk_range NeighborRange;
        GetGridNeighborRange(Grid, iNeighbor, &NeighborRange);
        if (AwayFromEdge) {
          ovk_range IntersectRange;
          ovkRangeIntersect(&NeighborRange, &DonorRange, &IntersectRange);
          ovkRangeUnion(&RecvFromNeighborDataRanges[iNeighbor], &IntersectRange,
            &RecvFromNeighborDataRanges[iNeighbor]);
        } else {
          for (k = DonorRange.b[2]; k < DonorRange.e[2]; ++k) {
            for (j = DonorRange.b[1]; j < DonorRange.e[1]; ++j) {
              for (i = DonorRange.b[0]; i < DonorRange.e[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                int AdjustedPoint[MAX_DIMS];
                ovkCartPeriodicAdjust(&Cart, Point, AdjustedPoint);
                if (ovkRangeContains(&NeighborRange, AdjustedPoint)) {
                  ovkRangeExtend(&RecvFromNeighborDataRanges[iNeighbor], AdjustedPoint,
                    &RecvFromNeighborDataRanges[iNeighbor]);
                }
              }
            }
          }
        }
      }
    }

    int NumCollectRecvs = 0;
    for (iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!ovkRangeIsEmpty(&RecvFromNeighborDataRanges[iNeighbor])) {
        ++NumCollectRecvs;
      }
    }

    int *CollectRecvIndexToNeighbor = malloc(NumCollectRecvs*sizeof(int));
    iCollectRecv = 0;
    for (iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!ovkRangeIsEmpty(&RecvFromNeighborDataRanges[iNeighbor])) {
        CollectRecvIndexToNeighbor[iCollectRecv] = iNeighbor;
        ++iCollectRecv;
      }
    }

    Exchange->num_collect_recvs = NumCollectRecvs;

    Exchange->collect_recv_source_ranks = malloc(NumCollectRecvs*sizeof(int));

    for (iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      int NeighborRank;
      GetGridNeighborRank(Grid, iNeighbor, &NeighborRank);
      Exchange->collect_recv_source_ranks[iCollectRecv] = NeighborRank;
    }

    bool **CollectRecvMasks = malloc(NumCollectRecvs*sizeof(bool *));
    for (iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      size_t NumPoints;
      ovkRangeCount(&RecvFromNeighborDataRanges[iNeighbor], &NumPoints);
      CollectRecvMasks[iCollectRecv] = malloc(NumPoints*sizeof(bool));
      for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
        CollectRecvMasks[iCollectRecv][iPoint] = false;
      }
    }

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS] = {
        Donors->extents[0][0][iDonor],
        Donors->extents[0][1][iDonor],
        Donors->extents[0][2][iDonor]
      };
      int DonorEnd[MAX_DIMS] = {
        Donors->extents[1][0][iDonor],
        Donors->extents[1][1][iDonor],
        Donors->extents[1][2][iDonor]
      };
      ovk_range DonorRange;
      ovkSetRange(&DonorRange, NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = ovkRangeIncludes(&GlobalRange, &DonorRange);
      for (iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
        iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
        ovk_range NeighborRange;
        GetGridNeighborRange(Grid, iNeighbor, &NeighborRange);
        if (AwayFromEdge) {
          ovk_range IntersectRange;
          ovkRangeIntersect(&NeighborRange, &DonorRange, &IntersectRange);
          for (k = IntersectRange.b[2]; k < IntersectRange.e[2]; ++k) {
            for (j = IntersectRange.b[1]; j < IntersectRange.e[1]; ++j) {
              for (i = IntersectRange.b[0]; i < IntersectRange.e[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                ovkRangeTupleToIndex(&RecvFromNeighborDataRanges[iNeighbor], OVK_COLUMN_MAJOR, Point,
                  &iPoint);
                CollectRecvMasks[iCollectRecv][iPoint] = true;
              }
            }
          }
        } else {
          for (k = DonorRange.b[2]; k < DonorRange.e[2]; ++k) {
            for (j = DonorRange.b[1]; j < DonorRange.e[1]; ++j) {
              for (i = DonorRange.b[0]; i < DonorRange.e[0]; ++i) {
                int Point[MAX_DIMS] = {i, j, k};
                int AdjustedPoint[MAX_DIMS];
                ovkCartPeriodicAdjust(&Cart, Point, AdjustedPoint);
                if (ovkRangeContains(&NeighborRange, AdjustedPoint)) {
                  ovkRangeTupleToIndex(&RecvFromNeighborDataRanges[iNeighbor], OVK_COLUMN_MAJOR,
                    AdjustedPoint, &iPoint);
                  CollectRecvMasks[iCollectRecv][iPoint] = true;
                }
              }
            }
          }
        }
      }
    }

    size_t **CollectRecvOffsets = malloc(NumCollectRecvs*sizeof(size_t *));
    for (iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      size_t NumPoints;
      ovkRangeCount(&RecvFromNeighborDataRanges[iNeighbor], &NumPoints);
      CollectRecvOffsets[iCollectRecv] = malloc(NumPoints*sizeof(size_t));
      size_t iRemotePoint = 0;
      for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
        if (CollectRecvMasks[iCollectRecv][iPoint]) {
          CollectRecvOffsets[iCollectRecv][iPoint] = iRemotePoint;
          ++iRemotePoint;
        }
      }
    }

    Exchange->num_collect_recv_points = malloc(NumCollectRecvs*sizeof(size_t));

    for (iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      size_t NumPoints;
      ovkRangeCount(&RecvFromNeighborDataRanges[iNeighbor], &NumPoints);
      size_t NumCollectRecvPoints = 0;
      for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
        if (CollectRecvMasks[iCollectRecv][iPoint]) {
          ++NumCollectRecvPoints;
        }
      }
      Exchange->num_collect_recv_points[iCollectRecv] = NumCollectRecvPoints;
    }

    Exchange->collect_recv_points = malloc(NumCollectRecvs*sizeof(int **));
    for (iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      Exchange->collect_recv_points[iCollectRecv] = malloc(MAX_DIMS*sizeof(int *));
      size_t NumCollectRecvPoints = Exchange->num_collect_recv_points[iCollectRecv];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Exchange->collect_recv_points[iCollectRecv][iDim] = malloc(NumCollectRecvPoints*
          sizeof(int));
      }
      iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
      size_t NumPoints;
      ovkRangeCount(&RecvFromNeighborDataRanges[iNeighbor], &NumPoints);
      iCollectRecvPoint = 0;
      for (iPoint = 0; iPoint < NumPoints; ++iPoint) {
        if (CollectRecvMasks[iCollectRecv][iPoint]) {
          int Point[MAX_DIMS];
          ovkRangeIndexToTuple(&RecvFromNeighborDataRanges[iNeighbor], OVK_COLUMN_MAJOR, iPoint,
            Point);
          for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
            Exchange->collect_recv_points[iCollectRecv][iDim][iCollectRecvPoint] = Point[iDim];
          }
          ++iCollectRecvPoint;
        }
      }

    }

    size_t TotalRemoteDonorPoints = 0;
    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DonorBegin[MAX_DIMS] = {
        Donors->extents[0][0][iDonor],
        Donors->extents[0][1][iDonor],
        Donors->extents[0][2][iDonor]
      };
      int DonorEnd[MAX_DIMS] = {
        Donors->extents[1][0][iDonor],
        Donors->extents[1][1][iDonor],
        Donors->extents[1][2][iDonor]
      };
      ovk_range DonorRange;
      ovkSetRange(&DonorRange, NumDims, DonorBegin, DonorEnd);
      ovk_range IntersectRange;
      ovkRangeIntersect(&LocalRange, &DonorRange, &IntersectRange);
      size_t NumDonorPoints, NumLocalDonorPoints, NumRemoteDonorPoints;
      ovkRangeCount(&DonorRange, &NumDonorPoints);
      ovkRangeCount(&IntersectRange, &NumLocalDonorPoints);
      NumRemoteDonorPoints = NumDonorPoints - NumLocalDonorPoints;
      Exchange->num_remote_donor_points[iDonor] = (int)NumRemoteDonorPoints;
      TotalRemoteDonorPoints += NumRemoteDonorPoints;
    }

    // Flattened buffers that elements of Exchange->remote_donor_points, etc. point to
    size_t *RemoteDonorPoints = malloc(TotalRemoteDonorPoints*sizeof(size_t));
    int *RemoteDonorPointCollectRecvIndices = malloc(TotalRemoteDonorPoints*sizeof(int));
    size_t *RemoteDonorPointCollectRecvBufferOffsets = malloc(TotalRemoteDonorPoints*
      sizeof(size_t));

    size_t RemoteDonorPointOffset = 0;
    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      size_t NumRemoteDonorPoints = Exchange->num_remote_donor_points[iDonor];
      Exchange->remote_donor_points[iDonor] = RemoteDonorPoints + RemoteDonorPointOffset;
      Exchange->remote_donor_point_collect_recv_indices[iDonor] =
        RemoteDonorPointCollectRecvIndices + RemoteDonorPointOffset;
      Exchange->remote_donor_point_collect_recv_buffer_offsets[iDonor] =
        RemoteDonorPointCollectRecvBufferOffsets + RemoteDonorPointOffset;
      RemoteDonorPointOffset += NumRemoteDonorPoints;
    }

    int MaxPointsInCell = 1;
    for (iDim = 0; iDim <= NumDims; ++iDim) {
      MaxPointsInCell *= MaxSize;
    }
    int *CellCollectRecvIndices = malloc(MaxPointsInCell*sizeof(int));
    int *CellCollectRecvOffsets = malloc(MaxPointsInCell*sizeof(size_t));

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      for (iDonorPoint = 0; iDonorPoint < MaxPointsInCell; ++iDonorPoint) {
        CellCollectRecvIndices[iDonorPoint] = -1;
      }
      int DonorBegin[MAX_DIMS] = {
        Donors->extents[0][0][iDonor],
        Donors->extents[0][1][iDonor],
        Donors->extents[0][2][iDonor]
      };
      int DonorEnd[MAX_DIMS] = {
        Donors->extents[1][0][iDonor],
        Donors->extents[1][1][iDonor],
        Donors->extents[1][2][iDonor]
      };
      ovk_range DonorRange;
      ovkSetRange(&DonorRange, NumDims, DonorBegin, DonorEnd);
      for (iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
        iNeighbor = CollectRecvIndexToNeighbor[iCollectRecv];
        ovk_range NeighborRange;
        GetGridNeighborRange(Grid, iNeighbor, &NeighborRange);
        ovk_range IntersectRange;
        ovkRangeIntersect(&NeighborRange, &DonorRange, &IntersectRange);
        for (k = IntersectRange.b[2]; k < IntersectRange.e[2]; ++k) {
          for (j = IntersectRange.b[1]; j < IntersectRange.e[1]; ++j) {
            for (i = IntersectRange.b[0]; i < IntersectRange.e[0]; ++i) {
              int Point[MAX_DIMS] = {i, j, k};
              ovkRangeTupleToIndexSmall(&DonorRange, OVK_COLUMN_MAJOR, Point, &iDonorPoint);
              ovkRangeTupleToIndex(&RecvFromNeighborDataRanges[iNeighbor], OVK_COLUMN_MAJOR, Point,
                &iPoint);
              CellCollectRecvIndices[iDonorPoint] = iCollectRecv;
              CellCollectRecvOffsets[iDonorPoint] = CollectRecvOffsets[iCollectRecv][iPoint];
            }
          }
        }
      }
      int iRemoteDonorPoint = 0;
      iDonorPoint = 0;
      for (k = DonorRange.b[2]; k < DonorRange.e[2]; ++k) {
        for (j = DonorRange.b[1]; j < DonorRange.e[1]; ++j) {
          for (i = DonorRange.b[0]; i < DonorRange.e[0]; ++i) {
            if (CellCollectRecvIndices[iDonorPoint] >= 0) {
              Exchange->remote_donor_points[iDonor][iRemoteDonorPoint] = iDonorPoint;
              Exchange->remote_donor_point_collect_recv_indices[iDonor][iRemoteDonorPoint] =
                CellCollectRecvIndices[iDonorPoint];
              Exchange->remote_donor_point_collect_recv_buffer_offsets[iDonor][iRemoteDonorPoint]
                = CellCollectRecvOffsets[iDonorPoint];
              ++iRemoteDonorPoint;
            }
            ++iDonorPoint;
          }
        }
      }
      Exchange->num_remote_donor_points[iDonor] = iRemoteDonorPoint;
    }

    free(RecvFromNeighborDataRanges);
    free(CollectRecvIndexToNeighbor);
    for (iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      free(CollectRecvMasks[iCollectRecv]);
    }
    free(CollectRecvMasks);
    for (iCollectRecv = 0; iCollectRecv < NumCollectRecvs; ++iCollectRecv) {
      free(CollectRecvOffsets[iCollectRecv]);
    }
    free(CollectRecvOffsets);
    free(CellCollectRecvIndices);
    free(CellCollectRecvOffsets);

  }

}

static void UpdateDonorsSorted(ovk_exchange *Exchange) {

  size_t iDonor;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  if (ovkRankHasConnectivityDonorSide(Connectivity)) {

    const ovk_connectivity_d *Donors;
    ovkGetConnectivityDonorSide(Connectivity, &Donors);

    const ovk_connectivity_d_properties *DonorsProperties;
    ovkGetConnectivityDonorSideProperties(Donors, &DonorsProperties);

    size_t NumDonors;
    ovkGetConnectivityDonorSidePropertyCount(DonorsProperties, &NumDonors);

    const ovk_grid_info *ReceiverGridInfo;
    ovkGetConnectivityReceiverGridInfo(Connectivity, &ReceiverGridInfo);

    ovk_range ReceiverGridGlobalRange;
    ovkGetGridInfoGlobalRange(ReceiverGridInfo, &ReceiverGridGlobalRange);

    size_t *DestinationIndices = malloc(NumDonors*sizeof(size_t));

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      int DestinationPoint[MAX_DIMS] = {
        Donors->destinations[0][iDonor],
        Donors->destinations[1][iDonor],
        Donors->destinations[2][iDonor]
      };
      ovkRangeTupleToIndex(&ReceiverGridGlobalRange, OVK_COLUMN_MAJOR, DestinationPoint,
        DestinationIndices+iDonor);
    }

    bool Sorted = true;

    // Check if they're already sorted
    size_t PrevIndex = 0;
    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      if (DestinationIndices[iDonor] < PrevIndex) {
        Sorted = false;
        break;
      }
      PrevIndex = DestinationIndices[iDonor];
    }

    if (Sorted) {
      for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
        Exchange->donors_sorted[iDonor] = iDonor;
      }
    } else {
      SortPermutation_size_t(NumDonors, DestinationIndices, Exchange->donors_sorted);
    }

    free(DestinationIndices);

  }

}

static void UpdateReceiversSorted(ovk_exchange *Exchange) {

  size_t iReceiver;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  if (ovkRankHasConnectivityReceiverSide(Connectivity)) {

    const ovk_connectivity_r *Receivers;
    ovkGetConnectivityReceiverSide(Connectivity, &Receivers);

    const ovk_connectivity_r_properties *ReceiversProperties;
    ovkGetConnectivityReceiverSideProperties(Receivers, &ReceiversProperties);

    size_t NumReceivers = 0;
    ovkGetConnectivityReceiverSidePropertyCount(ReceiversProperties, &NumReceivers);

    const ovk_grid *ReceiverGrid;
    ovkGetConnectivityReceiverSideGrid(Receivers, &ReceiverGrid);

    const ovk_grid_properties *ReceiverGridProperties;
    ovkGetGridProperties(ReceiverGrid, &ReceiverGridProperties);

    ovk_range GlobalRange;
    ovkGetGridPropertyGlobalRange(ReceiverGridProperties, &GlobalRange);

    size_t *PointIndices = malloc(NumReceivers*sizeof(size_t));

    for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      int Point[MAX_DIMS] = {
        Receivers->points[0][iReceiver],
        Receivers->points[1][iReceiver],
        Receivers->points[2][iReceiver]
      };
      ovkRangeTupleToIndex(&GlobalRange, OVK_COLUMN_MAJOR, Point, PointIndices+iReceiver);
    }

    bool Sorted = true;

    // Check if they're already sorted
    size_t PrevIndex = 0;
    for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      if (PointIndices[iReceiver] < PrevIndex) {
        Sorted = false;
        break;
      }
      PrevIndex = PointIndices[iReceiver];
    }

    if (Sorted) {
      for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
        Exchange->receivers_sorted[iReceiver] = iReceiver;
      }
    } else {
      SortPermutation_size_t(NumReceivers, PointIndices, Exchange->receivers_sorted);
    }

    free(PointIndices);

  }

}

static void UpdateSourceRanks(ovk_exchange *Exchange) {

  size_t iReceiver;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  bool ReceiverGridIsLocal = ovkRankHasConnectivityReceiverSide(Connectivity);

  const ovk_connectivity_r *Receivers;
  size_t NumReceivers;
  if (ReceiverGridIsLocal) {
    ovkGetConnectivityReceiverSide(Connectivity, &Receivers);
    const ovk_connectivity_r_properties *ReceiversProperties;
    ovkGetConnectivityReceiverSideProperties(Receivers, &ReceiversProperties);
    ovkGetConnectivityReceiverSidePropertyCount(ReceiversProperties, &NumReceivers);
  }

  t_ordered_map *Bins;
  OMCreate(&Bins);

  int *SourceBinIndices = NULL;
  if (ReceiverGridIsLocal) {
    SourceBinIndices = malloc(NumReceivers*sizeof(int));
    MapToPartitionBins(Exchange->source_hash, NumReceivers, (const int **)Receivers->sources,
      SourceBinIndices);
    for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      t_ordered_map_entry *Entry = OMFind(Bins, SourceBinIndices[iReceiver]);
      if (Entry == OMEnd(Bins)) {
        OMInsert(Bins, SourceBinIndices[iReceiver], NULL);
      }
    }
  }

  RetrievePartitionBins(Exchange->source_hash, Bins);

  if (ReceiverGridIsLocal) {
    FindPartitions(Exchange->source_hash, Bins, NumReceivers, (const int **)Receivers->sources,
      SourceBinIndices, Exchange->receiver_source_ranks);
    free(SourceBinIndices);
  }

  ClearPartitionBins(Bins);
  OMDestroy(&Bins);

}

static void UpdateDestRanks(ovk_exchange *Exchange) {

  size_t iDonor;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  bool DonorGridIsLocal = ovkRankHasConnectivityDonorSide(Connectivity);

  const ovk_connectivity_d *Donors;
  size_t NumDonors;
  if (DonorGridIsLocal) {
    ovkGetConnectivityDonorSide(Connectivity, &Donors);
    const ovk_connectivity_d_properties *DonorsProperties;
    ovkGetConnectivityDonorSideProperties(Donors, &DonorsProperties);
    ovkGetConnectivityDonorSidePropertyCount(DonorsProperties, &NumDonors);
  }

  t_ordered_map *Bins;
  OMCreate(&Bins);

  int *DestinationBinIndices = NULL;
  if (DonorGridIsLocal) {
    DestinationBinIndices = malloc(NumDonors*sizeof(int));
    MapToPartitionBins(Exchange->destination_hash, NumDonors, (const int **)Donors->destinations,
      DestinationBinIndices);
    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      t_ordered_map_entry *Entry = OMFind(Bins, DestinationBinIndices[iDonor]);
      if (Entry == OMEnd(Bins)) {
        OMInsert(Bins, DestinationBinIndices[iDonor], NULL);
      }
    }
  }

  RetrievePartitionBins(Exchange->destination_hash, Bins);

  if (DonorGridIsLocal) {
    FindPartitions(Exchange->destination_hash, Bins, NumDonors, (const int **)Donors->destinations,
      DestinationBinIndices, Exchange->donor_dest_ranks);
    free(DestinationBinIndices);
  }

  ClearPartitionBins(Bins);
  OMDestroy(&Bins);

}

static void UpdateSendInfo(ovk_exchange *Exchange) {

  int iSend;
  size_t iDonor;

  free_null(&Exchange->send_ranks);
  free_null(&Exchange->send_counts);

  Exchange->num_sends = 0;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  const ovk_connectivity_d *Donors;
  size_t NumDonors = 0;
  ovk_range LocalRange;
  if (ovkRankHasConnectivityDonorSide(Connectivity)) {
    ovkGetConnectivityDonorSide(Connectivity, &Donors);
    const ovk_connectivity_d_properties *DonorsProperties;
    ovkGetConnectivityDonorSideProperties(Donors, &DonorsProperties);
    ovkGetConnectivityDonorSidePropertyCount(DonorsProperties, &NumDonors);
    const ovk_grid *DonorGrid;
    ovkGetConnectivityDonorSideGrid(Donors, &DonorGrid);
    const ovk_grid_properties *DonorGridProperties;
    ovkGetGridProperties(DonorGrid, &DonorGridProperties);
    ovkGetGridPropertyLocalRange(DonorGridProperties, &LocalRange);
  }

  if (NumDonors > 0) {

    bool *DonorCommunicates = malloc(NumDonors*sizeof(bool));
    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      bool Communicates = Exchange->donor_dest_ranks[iDonor] >= 0;
      if (Communicates) {
        int DonorCell[MAX_DIMS] = {
          Donors->extents[0][0][iDonor],
          Donors->extents[0][1][iDonor],
          Donors->extents[0][2][iDonor]
        };
        Communicates = ovkRangeContains(&LocalRange, DonorCell);
      }
      DonorCommunicates[iDonor] = Communicates;
    }

    t_ordered_map *SendCounts;
    OMCreate(&SendCounts);

    t_ordered_map_entry *Entry;

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      if (DonorCommunicates[iDonor]) {
        Entry = OMFind(SendCounts, Exchange->donor_dest_ranks[iDonor]);
        if (Entry != OMEnd(SendCounts)) {
          size_t *SendCount = OMData(Entry);
          ++(*SendCount);
        } else {
          size_t *SendCount = malloc(sizeof(size_t));
          *SendCount = 1;
          OMInsert(SendCounts, Exchange->donor_dest_ranks[iDonor], SendCount);
        }
      }
    }

    int NumSends = OMSize(SendCounts);

    Exchange->num_sends = NumSends;

    Exchange->send_ranks = malloc(NumSends*sizeof(int));
    Exchange->send_counts = malloc(NumSends*sizeof(size_t));

    Entry = OMBegin(SendCounts);
    iSend = 0;
    while (Entry != OMEnd(SendCounts)) {
      int Rank = OMKey(Entry);
      size_t *SendCount = OMData(Entry);
      Exchange->send_ranks[iSend] = Rank;
      Exchange->send_counts[iSend] = *SendCount;
      free(SendCount);
      ++iSend;
      Entry = OMNext(Entry);
    }

    OMDestroy(&SendCounts);

    t_ordered_map *RankToSendIndex;
    OMCreate(&RankToSendIndex);

    for (iSend = 0; iSend < NumSends; ++iSend) {
      int *SendIndex = malloc(sizeof(int));
      *SendIndex = iSend;
      OMInsert(RankToSendIndex, Exchange->send_ranks[iSend], SendIndex);
    }

    for (iDonor = 0; iDonor < NumDonors; ++iDonor) {
      if (DonorCommunicates[iDonor]) {
        int *SendIndex = OMData(OMFind(RankToSendIndex, Exchange->donor_dest_ranks[iDonor]));
        Exchange->donor_send_indices[iDonor] = *SendIndex;
      } else {
        Exchange->donor_send_indices[iDonor] = -1;
      }
    }

    Entry = OMBegin(RankToSendIndex);
    while (Entry != OMEnd(RankToSendIndex)) {
      int *SendIndex = OMData(Entry);
      free(SendIndex);
      Entry = OMNext(Entry);
    }

    OMDestroy(&RankToSendIndex);

    free(DonorCommunicates);

  }

}

static void UpdateReceiveInfo(ovk_exchange *Exchange) {

  int iRecv;
  size_t iReceiver;

  free_null(&Exchange->recv_ranks);
  free_null(&Exchange->recv_counts);

  Exchange->num_recvs = 0;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  const ovk_connectivity_r *Receivers;
  size_t NumReceivers = 0;
  if (ovkRankHasConnectivityReceiverSide(Connectivity)) {
    ovkGetConnectivityReceiverSide(Connectivity, &Receivers);
    const ovk_connectivity_r_properties *ReceiversProperties;
    ovkGetConnectivityReceiverSideProperties(Receivers, &ReceiversProperties);
    ovkGetConnectivityReceiverSidePropertyCount(ReceiversProperties, &NumReceivers);
  }

  if (NumReceivers > 0) {

    t_ordered_map *ReceiveCounts;
    OMCreate(&ReceiveCounts);

    t_ordered_map_entry *Entry;

    for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      if (Exchange->receiver_source_ranks[iReceiver] >= 0) {
        Entry = OMFind(ReceiveCounts, Exchange->receiver_source_ranks[iReceiver]);
        if (Entry != OMEnd(ReceiveCounts)) {
          size_t *ReceiveCount = OMData(Entry);
          ++(*ReceiveCount);
        } else {
          size_t *ReceiveCount = malloc(sizeof(size_t));
          *ReceiveCount = 1;
          OMInsert(ReceiveCounts, Exchange->receiver_source_ranks[iReceiver], ReceiveCount);
        }
      }
    }

    int NumReceives = OMSize(ReceiveCounts);

    Exchange->num_recvs = NumReceives;

    Exchange->recv_ranks = malloc(NumReceives*sizeof(int));
    Exchange->recv_counts = malloc(NumReceives*sizeof(size_t));

    Entry = OMBegin(ReceiveCounts);
    iRecv = 0;
    while (Entry != OMEnd(ReceiveCounts)) {
      int Rank = OMKey(Entry);
      size_t *ReceiveCount = OMData(Entry);
      Exchange->recv_ranks[iRecv] = Rank;
      Exchange->recv_counts[iRecv] = *ReceiveCount;
      free(ReceiveCount);
      ++iRecv;
      Entry = OMNext(Entry);
    }

    OMDestroy(&ReceiveCounts);

    t_ordered_map *RankToReceiveIndex;
    OMCreate(&RankToReceiveIndex);

    for (iRecv = 0; iRecv < NumReceives; ++iRecv) {
      int *ReceiveIndex = malloc(sizeof(int));
      *ReceiveIndex = iRecv;
      OMInsert(RankToReceiveIndex, Exchange->recv_ranks[iRecv], ReceiveIndex);
    }

    for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      if (Exchange->receiver_source_ranks[iReceiver] >= 0) {
        int *ReceiveIndex = OMData(OMFind(RankToReceiveIndex,
          Exchange->receiver_source_ranks[iReceiver]));
        Exchange->receiver_recv_indices[iReceiver] = *ReceiveIndex;
      } else {
        Exchange->receiver_recv_indices[iReceiver] = -1;
      }
    }

    Entry = OMBegin(RankToReceiveIndex);
    while (Entry != OMEnd(RankToReceiveIndex)) {
      int *ReceiveIndex = OMData(Entry);
      free(ReceiveIndex);
      Entry = OMNext(Entry);
    }

    OMDestroy(&RankToReceiveIndex);

  }

}

#define PACK_COLLECT_SEND_BUFFER(input_type, output_type) \
  do { \
    for (iCount = 0; iCount < Count; ++iCount) { \
      ((output_type *)CollectSendBuffers[iSend][iCount])[iSendPoint] = (output_type)( \
      (const input_type *)GridData[iCount])[iGridPoint]; \
    } \
  } while (false)

void PRIVATE(ExchangeCollect)(const ovk_exchange *Exchange, ovk_data_type DataType, int Count,
  ovk_collect_op CollectOp, const void **GridData, ovk_array_layout GridDataLayout,
  void **DonorData) {

  int iDim, iSend, iRecv, iCount;
  size_t iDonor, iGridPoint, iBuffer;
  int i, j, k;

  int NumDims = Exchange->num_dims;
  const ovk_connectivity *Connectivity = Exchange->connectivity;

  const ovk_connectivity_d *Donors;
  ovkGetConnectivityDonorSide(Connectivity, &Donors);

  const ovk_connectivity_d_properties *DonorsProperties;
  ovkGetConnectivityDonorSideProperties(Donors, &DonorsProperties);

  size_t NumDonors;
  ovkGetConnectivityDonorSidePropertyCount(DonorsProperties, &NumDonors);

  int MaxSize;
  ovkGetConnectivityDonorSidePropertyMaxSize(DonorsProperties, &MaxSize);

  const ovk_grid *DonorGrid;
  ovkGetConnectivityDonorSideGrid(Donors, &DonorGrid);

  const ovk_grid_properties *DonorGridProperties;
  ovkGetGridProperties(DonorGrid, &DonorGridProperties);

  MPI_Comm DonorGridComm;
  ovkGetGridPropertyComm(DonorGridProperties, &DonorGridComm);

  ovk_range LocalRange;
  ovkGetGridPropertyLocalRange(DonorGridProperties, &LocalRange);

  int NumCollectSends = Exchange->num_collect_sends;
  int NumCollectRecvs = Exchange->num_collect_recvs;

  int DataSize = DataTypeSize(DataType);

  MPI_Datatype MPIDataType;
  int MPIDataSize;
  DataTypeToMPI(DataType, &MPIDataType, &MPIDataSize);

  void ***CollectSendBuffers = malloc(NumCollectSends*sizeof(void **));
  for (iSend = 0; iSend < NumCollectSends; ++iSend) {
    CollectSendBuffers[iSend] = malloc(Count*sizeof(void *));
    for (iCount = 0; iCount < Count; ++iCount) {
      CollectSendBuffers[iSend][iCount] = malloc(Exchange->num_collect_send_points[iSend] *
        MPIDataSize);
    }
  }

  void ***CollectRecvBuffers = malloc(NumCollectRecvs*sizeof(void **));
  for (iRecv = 0; iRecv < NumCollectRecvs; ++iRecv) {
    CollectRecvBuffers[iRecv] = malloc(Count*sizeof(void *));
    for (iCount = 0; iCount < Count; ++iCount) {
      CollectRecvBuffers[iRecv][iCount] = malloc(Exchange->num_collect_recv_points[iRecv] *
        MPIDataSize);
    }
  }

  MPI_Request *Requests = malloc(Count*(NumCollectSends+NumCollectRecvs)*sizeof(MPI_Request));

  int iRequest = 0;

  for (iRecv = 0; iRecv < NumCollectRecvs; ++iRecv) {
    for (iCount = 0; iCount < Count; ++iCount) {
      MPI_Irecv(CollectRecvBuffers[iRecv][iCount], Exchange->num_collect_recv_points[iRecv],
        MPIDataType, Exchange->collect_recv_source_ranks[iRecv], iCount, DonorGridComm,
        &Requests[iRequest]);
      ++iRequest;
    }
  }

  for (iSend = 0; iSend < NumCollectSends; ++iSend) {
    size_t NumSendPoints = Exchange->num_collect_send_points[iSend];
    size_t iSendPoint;
    for (iSendPoint = 0; iSendPoint < NumSendPoints; ++iSendPoint) {
      int Point[MAX_DIMS] = {
        Exchange->collect_send_points[iSend][0][iSendPoint],
        Exchange->collect_send_points[iSend][1][iSendPoint],
        Exchange->collect_send_points[iSend][2][iSendPoint]
      };
      ovkRangeTupleToIndex(&LocalRange, GridDataLayout, Point, &iGridPoint);
      switch (DataType) {
      case OVK_BOOL:
        PACK_COLLECT_SEND_BUFFER(bool, unsigned char);
        break;
      case OVK_DOUBLE:
        PACK_COLLECT_SEND_BUFFER(double, double);
        break;
      default:
        OVK_DEBUG_ASSERT(false, "Data type not yet implemented.");
        break;
      }
    }
  }

  for (iSend = 0; iSend < NumCollectSends; ++iSend) {
    for (iCount = 0; iCount < Count; ++iCount) {
      MPI_Isend(CollectSendBuffers[iSend][iCount], Exchange->num_collect_send_points[iSend],
        MPIDataType, Exchange->collect_send_dest_ranks[iSend], iCount, DonorGridComm,
        &Requests[iRequest]);
      ++iRequest;
    }
  }

  MPI_Waitall(Count*(NumCollectSends+NumCollectRecvs), Requests, MPI_STATUSES_IGNORE);

  for (iSend = 0; iSend < NumCollectSends; ++iSend) {
    for (iCount = 0; iCount < Count; ++iCount) {
      free(CollectSendBuffers[iSend][iCount]);
    }
    free(CollectSendBuffers[iSend]);
  }
  free(CollectSendBuffers);

  free(Requests);

  int MaxPointsInCell = 1;
  for (iDim = 0; iDim < NumDims; ++iDim) {
    MaxPointsInCell *= MaxSize;
  }

  void **DonorPointData = malloc(Count*sizeof(void *));
  for (iCount = 0; iCount < Count; ++iCount) {
    DonorPointData[iCount] = malloc(MaxPointsInCell*DataSize);
  }

  for (iDonor = 0; iDonor < NumDonors; ++iDonor) {

    int DonorBegin[MAX_DIMS] = {
      Donors->extents[0][0][iDonor],
      Donors->extents[0][1][iDonor],
      Donors->extents[0][2][iDonor]
    };
    int DonorEnd[MAX_DIMS] = {
      Donors->extents[1][0][iDonor],
      Donors->extents[1][1][iDonor],
      Donors->extents[1][2][iDonor]
    };

    ovk_range DonorRange;
    ovkSetRange(&DonorRange, NumDims, DonorBegin, DonorEnd);

    ovk_range LocalDonorRange;
    ovkRangeIntersect(&LocalRange, &DonorRange, &LocalDonorRange);

    int NumRemotePoints;
    int iDonorPoint, iRemotePoint;

    switch (DataType) {
    case OVK_DOUBLE:
      // Fill in the local data
      for (k = LocalDonorRange.b[2]; k < LocalDonorRange.e[2]; ++k) {
        for (j = LocalDonorRange.b[1]; j < LocalDonorRange.e[1]; ++j) {
          for (i = LocalDonorRange.b[0]; i < LocalDonorRange.e[0]; ++i) {
            int Point[MAX_DIMS] = {i, j, k};
            ovkRangeTupleToIndexSmall(&DonorRange, OVK_COLUMN_MAJOR, Point, &iDonorPoint);
            ovkRangeTupleToIndex(&LocalRange, GridDataLayout, Point, &iGridPoint);
            for (iCount = 0; iCount < Count; ++iCount) {
              ((double *)DonorPointData[iCount])[iDonorPoint] =
                ((const double *)GridData[iCount])[iGridPoint];
            }
          }
        }
      }
      // Fill in the remote data from recv buffers
      NumRemotePoints = Exchange->num_remote_donor_points[iDonor];
      for (iRemotePoint = 0; iRemotePoint < NumRemotePoints; ++iRemotePoint) {
        iDonorPoint = Exchange->remote_donor_points[iDonor][iRemotePoint];
        iRecv = Exchange->remote_donor_point_collect_recv_indices[iDonor][iRemotePoint];
        iBuffer = Exchange->remote_donor_point_collect_recv_buffer_offsets[iDonor][iRemotePoint];
        for (iCount = 0; iCount < Count; ++iCount) {
          ((double *)DonorPointData[iCount])[iDonorPoint] =
            ((const double *)CollectRecvBuffers[iRecv][iCount])[iBuffer];
        }
      }
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Data type not yet implemented.");
      break;
    }

    switch (CollectOp) {
    case OVK_COLLECT_INTERPOLATE:
      switch (DataType) {
      case OVK_DOUBLE:
        iDonorPoint = 0;
        for (iCount = 0; iCount < Count; ++iCount) {
          ((double *)DonorData[iCount])[iDonor] = 0.;
        }
        for (k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
          for (j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
            for (i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
              int Point[MAX_DIMS] = {i, j, k};
              double Coef = 1.;
              for (iDim = 0; iDim < NumDims; ++iDim) {
                Coef *= Donors->interp_coefs[iDim][Point[iDim]-DonorBegin[iDim]][iDonor];
              }
              for (iCount = 0; iCount < Count; ++iCount) {
                ((double *)DonorData[iCount])[iDonor] += Coef *
                  ((double *)DonorPointData[iCount])[iDonorPoint];
              }
              ++iDonorPoint;
            }
          }
        }
        break;
      default:
        OVK_DEBUG_ASSERT(false, "Data type not yet implemented.");
        break;
      }
      break;
    default:
      OVK_DEBUG_ASSERT(false, "Collect operation type not yet implemented.");
      break;
    }

  }

  for (iRecv = 0; iRecv < NumCollectRecvs; ++iRecv) {
    for (iCount = 0; iCount < Count; ++iCount) {
      free(CollectRecvBuffers[iRecv][iCount]);
    }
    free(CollectRecvBuffers[iRecv]);
  }
  free(CollectRecvBuffers);

  for (iCount = 0; iCount < Count; ++iCount) {
    free(DonorPointData[iCount]);
  }
  free(DonorPointData);

}

void PRIVATE(ExchangeSend)(const ovk_exchange *Exchange, ovk_data_type DataType, int Count,
  const void **DonorData, int Tag, ovk_request **Request) {

  int iSend, iCount;
  size_t iDonor, iDonorOrder;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  const ovk_connectivity_d *Donors;
  ovkGetConnectivityDonorSide(Connectivity, &Donors);

  const ovk_connectivity_d_properties *DonorsProperties;
  ovkGetConnectivityDonorSideProperties(Donors, &DonorsProperties);

  size_t NumDonors;
  ovkGetConnectivityDonorSidePropertyCount(DonorsProperties, &NumDonors);

  int NumSends = Exchange->num_sends;

  MPI_Datatype MPIDataType;
  int MPIDataSize;
  DataTypeToMPI(DataType, &MPIDataType, &MPIDataSize);

  CreateSendRequest(Exchange, MPIDataSize, Count, Request);
  t_send_request_data *RequestData = (*Request)->data;

  size_t *NextBufferEntry = malloc(NumSends*sizeof(size_t));
  for (iSend = 0; iSend < NumSends; ++iSend) {
    NextBufferEntry[iSend] = 0;
  }

  switch (DataType) {
  case OVK_BOOL:
    for (iDonorOrder = 0; iDonorOrder < NumDonors; ++iDonorOrder) {
      iDonor = Exchange->donors_sorted[iDonorOrder];
      iSend = Exchange->donor_send_indices[iDonor];
      if (iSend >= 0) {
        unsigned char *BufferData = (unsigned char *)RequestData->buffers[iSend] +
          Count*NextBufferEntry[iSend];
        for (iCount = 0; iCount < Count; ++iCount) {
          BufferData[iCount] = (unsigned char)((const bool *)DonorData[iCount])[iDonor];
        }
        ++NextBufferEntry[iSend];
      }
    }
    break;
  case OVK_DOUBLE:
    for (iDonorOrder = 0; iDonorOrder < NumDonors; ++iDonorOrder) {
      iDonor = Exchange->donors_sorted[iDonorOrder];
      iSend = Exchange->donor_send_indices[iDonor];
      if (iSend >= 0) {
        double *BufferData = (double *)RequestData->buffers[iSend] + Count*NextBufferEntry[iSend];
        for (iCount = 0; iCount < Count; ++iCount) {
          BufferData[iCount] = ((const double *)DonorData[iCount])[iDonor];
        }
        ++NextBufferEntry[iSend];
      }
    }
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Data type not yet implemented.");
    break;
  }

  free(NextBufferEntry);

  for (iSend = 0; iSend < NumSends; ++iSend) {
    MPI_Isend(RequestData->buffers[iSend], Count*Exchange->send_counts[iSend], MPIDataType,
      Exchange->send_ranks[iSend], Tag, Exchange->comm, &RequestData->mpi_requests[iSend]);
  }

}

void PRIVATE(ExchangeReceive)(const ovk_exchange *Exchange, ovk_data_type DataType, int Count,
  void **ReceiverData, int Tag, ovk_request **Request) {

  int iRecv;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  const ovk_connectivity_r *Receivers;
  ovkGetConnectivityReceiverSide(Connectivity, &Receivers);

  const ovk_connectivity_r_properties *ReceiversProperties;
  ovkGetConnectivityReceiverSideProperties(Receivers, &ReceiversProperties);

  size_t NumReceivers;
  ovkGetConnectivityReceiverSidePropertyCount(ReceiversProperties, &NumReceivers);

  int NumRecvs = Exchange->num_recvs;

  MPI_Datatype MPIDataType;
  int MPIDataSize;
  DataTypeToMPI(DataType, &MPIDataType, &MPIDataSize);

  CreateReceiveRequest(Exchange, DataType, MPIDataSize, Count, ReceiverData, Request);
  t_recv_request_data *RequestData = (*Request)->data;

  for (iRecv = 0; iRecv < NumRecvs; ++iRecv) {
    MPI_Irecv(RequestData->buffers[iRecv], Count*Exchange->recv_counts[iRecv], MPIDataType,
      Exchange->recv_ranks[iRecv], Tag, Exchange->comm, &RequestData->mpi_requests[iRecv]);
  }

}

static void CompleteReceive(t_recv_request_data *RequestData) {

  int iRecv, iCount;
  size_t iReceiver, iReceiverOrder;

  const ovk_exchange *Exchange = RequestData->exchange;
  const ovk_connectivity *Connectivity = Exchange->connectivity;

  const ovk_connectivity_properties *ConnectivityProperties;
  ovkGetConnectivityProperties(Connectivity, &ConnectivityProperties);

  int DonorGridID;
  ovkGetConnectivityPropertyDonorGridID(ConnectivityProperties, &DonorGridID);

  int ReceiverGridID;
  ovkGetConnectivityPropertyReceiverGridID(ConnectivityProperties, &ReceiverGridID);

  ovk_data_type DataType = RequestData->data_type;
  int Count = RequestData->count;
  void **ReceiverData = RequestData->receiver_data;

  const ovk_connectivity_r *Receivers;
  ovkGetConnectivityReceiverSide(Connectivity, &Receivers);

  const ovk_connectivity_r_properties *ReceiversProperties;
  ovkGetConnectivityReceiverSideProperties(Receivers, &ReceiversProperties);

  size_t NumReceivers;
  ovkGetConnectivityReceiverSidePropertyCount(ReceiversProperties, &NumReceivers);

  int NumRecvs = Exchange->num_recvs;

  size_t *NextBufferEntry = malloc(NumRecvs*sizeof(size_t));
  for (iRecv = 0; iRecv < NumRecvs; ++iRecv) {
    NextBufferEntry[iRecv] = 0;
  }

  switch (DataType) {
  case OVK_BOOL:
    for (iReceiverOrder = 0; iReceiverOrder < NumReceivers; ++iReceiverOrder) {
      iReceiver = Exchange->receivers_sorted[iReceiverOrder];
      iRecv = Exchange->receiver_recv_indices[iReceiver];
      if (iRecv >= 0) {
        unsigned char *BufferData = (unsigned char *)RequestData->buffers[iRecv] +
          Count*NextBufferEntry[iRecv];
        for (iCount = 0; iCount < Count; ++iCount) {
          ((bool *)ReceiverData[iCount])[iReceiver] = (bool)BufferData[iCount];
        }
        ++NextBufferEntry[iRecv];
      }
    }
    break;
  case OVK_DOUBLE:
    for (iReceiverOrder = 0; iReceiverOrder < NumReceivers; ++iReceiverOrder) {
      iReceiver = Exchange->receivers_sorted[iReceiverOrder];
      iRecv = Exchange->receiver_recv_indices[iReceiver];
      if (iRecv >= 0) {
        double *BufferData = (double *)RequestData->buffers[iRecv] + Count*NextBufferEntry[iRecv];
        for (iCount = 0; iCount < Count; ++iCount) {
          ((double *)ReceiverData[iCount])[iReceiver] = BufferData[iCount];
        }
        ++NextBufferEntry[iRecv];
      }
    }
    break;
  default:
    OVK_DEBUG_ASSERT(false, "Data type not yet implemented.");
    break;
  }

  free(NextBufferEntry);

}

void PRIVATE(ExchangeWaitAll)(int NumRequests, ovk_request **Requests) {

  int iRequest, iSend, iRecv;

  int NumMPIRequests = 0;
  for (iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests[iRequest]) {
      switch (Requests[iRequest]->type) {
      case SEND_REQUEST:
        NumMPIRequests += ((t_send_request_data *)(Requests[iRequest]->data))->exchange->num_sends;
        break;
      case RECV_REQUEST:
        NumMPIRequests += ((t_recv_request_data *)(Requests[iRequest]->data))->exchange->num_recvs;
        break;
      }
    }
  }

  MPI_Request *MPIRequests = malloc(NumMPIRequests*sizeof(MPI_Request));

  int iMPIRequest = 0;
  for (iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests[iRequest]) {
      switch (Requests[iRequest]->type) {
      case SEND_REQUEST:
        {
          t_send_request_data *RequestData = Requests[iRequest]->data;
          const ovk_exchange *Exchange = RequestData->exchange;
          for (iSend = 0; iSend < Exchange->num_sends; ++iSend) {
            MPIRequests[iMPIRequest] = RequestData->mpi_requests[iSend];
            ++iMPIRequest;
          }
        }
        break;
      case RECV_REQUEST:
        {
          t_recv_request_data *RequestData = Requests[iRequest]->data;
          const ovk_exchange *Exchange = RequestData->exchange;
          for (iRecv = 0; iRecv < Exchange->num_recvs; ++iRecv) {
            MPIRequests[iMPIRequest] = RequestData->mpi_requests[iRecv];
            ++iMPIRequest;
          }
        }
        break;
      }
    }
  }

  MPI_Waitall(NumMPIRequests, MPIRequests, MPI_STATUSES_IGNORE);

  free(MPIRequests);

  for (iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests[iRequest]) {
      switch (Requests[iRequest]->type) {
      case SEND_REQUEST:
        DestroySendRequest(&Requests[iRequest]);
        break;
      case RECV_REQUEST:
        CompleteReceive((t_recv_request_data *)(Requests[iRequest]->data));
        DestroyReceiveRequest(&Requests[iRequest]);
        break;
      }
    }
  }

}

void PRIVATE(ExchangeWaitAny)(int NumRequests, ovk_request **Requests, int *Index) {

  OVK_DEBUG_ASSERT(false, "ovkWaitAny not yet implemented.");

}

void PRIVATE(ExchangeDisperse)(const ovk_exchange *Exchange, ovk_data_type DataType, int Count,
  ovk_disperse_op DisperseOp, const void **ReceiverData, void **GridData,
  ovk_array_layout GridDataLayout) {

  int iCount;
  size_t iReceiver, iPoint;

  const ovk_connectivity *Connectivity = Exchange->connectivity;

  const ovk_connectivity_r *Receivers;
  ovkGetConnectivityReceiverSide(Connectivity, &Receivers);

  const ovk_connectivity_r_properties *ReceiversProperties;
  ovkGetConnectivityReceiverSideProperties(Receivers, &ReceiversProperties);

  size_t NumReceivers;
  ovkGetConnectivityReceiverSidePropertyCount(ReceiversProperties, &NumReceivers);

  const ovk_grid *ReceiverGrid;
  ovkGetConnectivityReceiverSideGrid(Receivers, &ReceiverGrid);

  const ovk_grid_properties *ReceiverGridProperties;
  ovkGetGridProperties(ReceiverGrid, &ReceiverGridProperties);

  ovk_range LocalRange;
  ovkGetGridPropertyLocalRange(ReceiverGridProperties, &LocalRange);

  switch (DisperseOp) {
  case OVK_DISPERSE_OVERWRITE:
    for (iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      int Point[MAX_DIMS] = {
        Receivers->points[0][iReceiver],
        Receivers->points[1][iReceiver],
        Receivers->points[2][iReceiver]
      };
      ovkRangeTupleToIndex(&LocalRange, GridDataLayout, Point, &iPoint);
      switch (DataType) {
      case OVK_BOOL:
        for (iCount = 0; iCount < Count; ++iCount) {
          ((bool *)GridData[iCount])[iPoint] = ((const bool *)ReceiverData[iCount])[iReceiver];
        }
        break;
      case OVK_DOUBLE:
        for (iCount = 0; iCount < Count; ++iCount) {
          ((double *)GridData[iCount])[iPoint] = ((const double *)ReceiverData[iCount])[iReceiver];
        }
        break;
      default:
        OVK_DEBUG_ASSERT(false, "Data type not yet implemented.");
        break;
      }
    }
    break;
  }

}

static void CreateSendRequest(const ovk_exchange *Exchange, int DataSize, int Count,
  ovk_request **Request) {

  int iSend;

  *Request = malloc(sizeof(ovk_request));

  (*Request)->type = SEND_REQUEST;
  (*Request)->data = malloc(sizeof(t_send_request_data));

  t_send_request_data *RequestData = (*Request)->data;

  RequestData->exchange = Exchange;

  int NumSends = Exchange->num_sends;

  RequestData->buffers = malloc(NumSends*sizeof(void *));
  for (iSend = 0; iSend < NumSends; ++iSend) {
    RequestData->buffers[iSend] = malloc(DataSize*Count*Exchange->send_counts[iSend]);
  }

  RequestData->mpi_requests = malloc(NumSends*sizeof(MPI_Request));
  for (iSend = 0; iSend < NumSends; ++iSend) {
    RequestData->mpi_requests[iSend] = MPI_REQUEST_NULL;
  }

}

static void DestroySendRequest(ovk_request **Request) {

  int iSend;

  t_send_request_data *RequestData = (*Request)->data;

  const ovk_exchange *Exchange = RequestData->exchange;

  int NumSends = Exchange->num_sends;

  for (iSend = 0; iSend < NumSends; ++iSend) {
    free(RequestData->buffers[iSend]);
  }
  free(RequestData->buffers);

  free(RequestData->mpi_requests);

  free(RequestData);

  free_null(Request);

}

static void CreateReceiveRequest(const ovk_exchange *Exchange, ovk_data_type DataType, int DataSize,
  int Count, void **ReceiverData, ovk_request **Request) {

  int iRecv, iCount;

  *Request = malloc(sizeof(ovk_request));

  (*Request)->type = RECV_REQUEST;
  (*Request)->data = malloc(sizeof(t_recv_request_data));

  t_recv_request_data *RequestData = (*Request)->data;

  RequestData->exchange = Exchange;

  int NumRecvs = Exchange->num_recvs;

  RequestData->buffers = malloc(NumRecvs*sizeof(void *));
  for (iRecv = 0; iRecv < NumRecvs; ++iRecv) {
    RequestData->buffers[iRecv] = malloc(DataSize*Count*Exchange->recv_counts[iRecv]);
  }

  RequestData->mpi_requests = malloc(NumRecvs*sizeof(MPI_Request));
  for (iRecv = 0; iRecv < NumRecvs; ++iRecv) {
    RequestData->mpi_requests[iRecv] = MPI_REQUEST_NULL;
  }

  RequestData->data_type = DataType;
  RequestData->count = Count;

  RequestData->receiver_data = malloc(Count*sizeof(void *));
  for (iCount = 0; iCount < Count; ++iCount) {
    RequestData->receiver_data[iCount] = ReceiverData[iCount];
  }

}

static void DestroyReceiveRequest(ovk_request **Request) {

  int iRecv;

  t_recv_request_data *RequestData = (*Request)->data;

  const ovk_exchange *Exchange = RequestData->exchange;

  int NumRecvs = Exchange->num_recvs;

  for (iRecv = 0; iRecv < NumRecvs; ++iRecv) {
    free(RequestData->buffers[iRecv]);
  }
  free(RequestData->buffers);

  free(RequestData->mpi_requests);

  free(RequestData->receiver_data);

  free(RequestData);

  free_null(Request);

}
