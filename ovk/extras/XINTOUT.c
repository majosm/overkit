// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

// =======
//  Notes
// =======
// * Current implementation will not play well with large datasets on small numbers of ranks
//   (MPI I/O and send/recv routines will fail if data has count larger than INT_MAX). Can be fixed
//   by splitting into multiple calls of size INT_MAX.

#include "ovk/extras/XINTOUT.h"

#include "ovk/extras/Global.h"
#include "ovk/core/Cart.h"
#include "ovk/core/Connectivity.h"
#include "ovk/core/Domain.h"
#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Logger.h"
#include "ovk/core/MPIUtils.h"
#include "ovk/core/PartitionHash.h"
#include "ovk/core/ProfileUtils.h"
#include "ovk/core/Range.h"
#include "ovk/core/TextUtils.h"

#include <limits.h>
#include <stdio.h>
#include <string.h>

typedef struct {
  size_t count;
  int max_size;
  int *extents[2][MAX_DIMS];
  double *coords[MAX_DIMS];
  double **interp_coefs[MAX_DIMS];
  int *destination_grid_ids;
  int *destination_points[MAX_DIMS];
} t_donor_data;

typedef struct {
  size_t count;
  int *points[MAX_DIMS];
  int *source_grid_ids;
  int *source_cells[MAX_DIMS];
} t_receiver_data;

typedef struct {
  size_t count;
  size_t *connection_ids;
  int *donor_grid_ids;
  int *donor_cells[MAX_DIMS];
  int *receiver_grid_ids;
  int *receiver_points[MAX_DIMS];
} t_connection_data;

typedef struct {
  size_t begin;
  size_t end;
  t_donor_data *data;
  size_t starting_connection_id;
} t_xintout_donor_chunk;

typedef struct {
  size_t begin;
  size_t end;
  t_receiver_data *data;
  size_t *connection_ids;
} t_xintout_receiver_chunk;

typedef struct {
  size_t count;
  size_t chunk_size;
  bool has_chunk;
  t_xintout_donor_chunk *chunk;
} t_xintout_donors;

typedef struct {
  size_t count;
  size_t chunk_size;
  bool has_chunk;
  t_xintout_receiver_chunk *chunk;
} t_xintout_receivers;

typedef struct {
  int id;
  char name[OVK_NAME_LENGTH];
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  t_logger *logger;
  t_error_handler *error_handler;
  int global_size[MAX_DIMS];
  t_xintout_donors donors;
  t_xintout_receivers receivers;
} t_xintout_grid;

typedef struct {
  size_t begin;
  size_t end;
  t_connection_data *data;
} t_xintout_connection_bin;

typedef struct {
  size_t count;
  size_t bin_size;
  bool has_bin;
  t_xintout_connection_bin *bin;
} t_xintout_connections;

typedef struct {
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  t_logger *logger;
  t_error_handler *error_handler;
  int num_grids;
  int num_local_grids;
  t_xintout_grid **grids;
  t_xintout_connections connections;
} t_xintout;

static void CreateXINTOUT(t_xintout **XINTOUT, int NumDims, MPI_Comm Comm, int NumGrids,
  int NumLocalGrids, const int *LocalGridIDs, const char (*LocalGridNames)[OVK_NAME_LENGTH],
  const MPI_Comm *LocalGridComms, const int **LocalGridGlobalSizes, t_logger *Logger,
  t_error_handler *ErrorHandler);
static void DestroyXINTOUT(t_xintout **XINTOUT);

static void CreateXINTOUTGrid(t_xintout_grid **XINTOUTGrid, int ID, const char *Name, int NumDims,
  MPI_Comm Comm, const int *GlobalSize, t_logger *Logger, t_error_handler *ErrorHandler);
static void DestroyXINTOUTGrid(t_xintout_grid **XINTOUTGrid);

static ovk_error ReadXINTOUT(t_xintout *XINTOUT, const char *HOPath, const char *XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo, t_profiler *Profiler);

static ovk_error ReadGlobalInfo(const t_xintout *XINTOUT, const char *HOPath,
  ovk_ext_endian *Endian, ovk_ext_xintout_format *Format, bool *WithIBlank, t_profiler *Profiler);

static bool DetectFormat(MPI_File HOFile, ovk_ext_endian *Endian, ovk_ext_xintout_format *Format,
  t_profiler *Profiler);

static ovk_error ReadGridInfo(const t_xintout_grid *XINTOUTGrid, const char *HOPath,
  const char *XPath, size_t *NumDonors, size_t *NumReceivers, size_t *StartingConnectionID,
  MPI_Offset *HODonorCellsOffset, MPI_Offset *HODonorCoordsOffset,
  MPI_Offset *HOReceiverPointsOffset, MPI_Offset *HOReceiverConnectionIDsOffset,
  MPI_Offset *XDonorSizesOffset, MPI_Offset *XDonorInterpCoefsOffset, ovk_ext_endian Endian,
  ovk_ext_xintout_format Format, bool WithIBlank, t_profiler *Profiler);

static ovk_error ReadDonors(t_xintout_grid *XINTOUTGrid, const char *HOPath, const char *XPath,
  size_t NumDonors, size_t StartingConnectionID, MPI_Offset HOCellsOffset, MPI_Offset HOCoordsOffset,
  MPI_Offset XSizesOffset, MPI_Offset XInterpCoefsOffset, ovk_ext_endian Endian,
  int ReadGranularityAdjust, MPI_Info MPIInfo, t_profiler *Profiler);

static ovk_error ReadReceivers(t_xintout_grid *XINTOUTGrid, const char *HOPath, size_t NumReceivers,
  MPI_Offset HOPointsOffset, MPI_Offset HOConnectionIDsOffset, ovk_ext_endian Endian,
  ovk_ext_xintout_format Format, int ReadGranularityAdjust, MPI_Info MPIInfo, t_profiler *Profiler);

static void MatchDonorsAndReceivers(t_xintout *XINTOUT, t_profiler *Profiler);

static void DistributeConnectivityData(const t_xintout *XINTOUT, const ovk_grid **LocalGrids,
  t_donor_data **LocalDonorData, t_receiver_data **LocalReceiverData, t_profiler *Profiler);

static void DistributeGridConnectivityData(const t_xintout_grid *XINTOUTGrid, const ovk_grid *Grid,
  t_donor_data **DonorData, t_receiver_data **ReceiverData, t_profiler *Profiler);

static void ImportConnectivityData(int NumGrids, int NumLocalGrids, int *LocalGridIDs,
  const t_donor_data **LocalDonorData, const t_receiver_data **LocalReceiverData,
  ovk_domain *Domain);

static void ImportDonors(const t_donor_data *DonorData, int NumReceiverGrids,
  ovk_connectivity_d **Donors);
static void ImportReceivers(const t_receiver_data *ReceiverData, int NumDonorGrids,
  ovk_connectivity_r **Receivers);

static void CreateDonorData(t_donor_data **Data, size_t Count, int MaxSize);
static void DestroyDonorData(t_donor_data **Data);

static void CreateReceiverData(t_receiver_data **Data, size_t Count);
static void DestroyReceiverData(t_receiver_data **Data);

static void CreateConnectionData(t_connection_data **Data, size_t Count);
static void DestroyConnectionData(t_connection_data **Data);

static size_t BinDivide(size_t N, int NumChunks);
static void Chunkify(size_t Count, int MaxChunks, size_t TargetChunkSize, int Adjust,
  int *ChunkInterval, int *NumChunks, size_t *ChunkSize);

static int File_read_all_endian(MPI_File File, void *Buffer, int Count, MPI_Datatype DataType,
  ovk_ext_endian Endian, MPI_Status *Status, t_profiler *Profiler, MPI_Comm Comm);
static int File_read_at_endian(MPI_File File, MPI_Offset Offset, void *Buffer, int Count,
  MPI_Datatype DataType, ovk_ext_endian Endian, MPI_Status *Status, t_profiler *Profiler);

static ovk_ext_endian MachineEndian();
static void SwapEndian(void *Data, int ElementSize, int NumElements);

ovk_error ovkEXTImportXINTOUT(ovk_domain *Domain, const char *HOPath, const char *XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(HOPath, "Invalid HO path pointer.");
  OVK_DEBUG_ASSERT(XPath, "Invalid X path pointer.");

  int iGrid, iLocalGrid;
  int iDim;

  t_logger *Logger;
  GetDomainLogger(Domain, &Logger);

  t_error_handler *ErrorHandler;
  GetDomainErrorHandler(Domain, &ErrorHandler);

  OVK_EH_INIT(ErrorHandler);
  ovk_error Error = OVK_NO_ERROR;

  const ovk_domain_properties *DomainProperties;
  ovkGetDomainProperties(Domain, &DomainProperties);

  int NumDims;
  int NumGrids;
  MPI_Comm Comm;
  int CommRank;
  ovkGetDomainPropertyDimension(DomainProperties, &NumDims);
  ovkGetDomainPropertyGridCount(DomainProperties, &NumGrids);
  ovkGetDomainPropertyComm(DomainProperties, &Comm);
  ovkGetDomainPropertyCommRank(DomainProperties, &CommRank);

  t_profiler *Profiler;
  CreateProfiler(&Profiler, Comm);
  if (OVK_PROFILE) EnableProfiler(Profiler);
  AddProfilerTimer(Profiler, "XINTOUT::Create");
  AddProfilerTimer(Profiler, "XINTOUT::Destroy");
  AddProfilerTimer(Profiler, "XINTOUT::Read");
  AddProfilerTimer(Profiler, "XINTOUT::Read::MPI-IO::Open");
  AddProfilerTimer(Profiler, "XINTOUT::Read::MPI-IO::Close");
  AddProfilerTimer(Profiler, "XINTOUT::Read::MPI-IO::Read");
  AddProfilerTimer(Profiler, "XINTOUT::Read::MPI-IO::Other");
  AddProfilerTimer(Profiler, "XINTOUT::Match");
  AddProfilerTimer(Profiler, "XINTOUT::Match::MapToBins");
  AddProfilerTimer(Profiler, "XINTOUT::Match::Handshake");
  AddProfilerTimer(Profiler, "XINTOUT::Match::SendToBins");
  AddProfilerTimer(Profiler, "XINTOUT::Match::FillConnectionData");
  AddProfilerTimer(Profiler, "XINTOUT::Match::RecvFromBins");
  AddProfilerTimer(Profiler, "XINTOUT::Match::Unpack");
  AddProfilerTimer(Profiler, "XINTOUT::Distribute");
  AddProfilerTimer(Profiler, "XINTOUT::Distribute::MapToBins");
  AddProfilerTimer(Profiler, "XINTOUT::Distribute::RetrieveBins");
  AddProfilerTimer(Profiler, "XINTOUT::Distribute::FindRanks");
  AddProfilerTimer(Profiler, "XINTOUT::Distribute::Handshake");
  AddProfilerTimer(Profiler, "XINTOUT::Distribute::SendData");
  AddProfilerTimer(Profiler, "XINTOUT::Import");

  MPI_Barrier(Comm);

  if (NumGrids > 0) {

    int NumLocalGrids = 0;
    for (iGrid = 0; iGrid < NumGrids; ++iGrid) {
      int GridID = iGrid+1;
      if (ovkRankHasGrid(Domain, GridID)) {
        ++NumLocalGrids;
      }
    }

    const ovk_grid **LocalGrids = malloc(NumLocalGrids*sizeof(const ovk_grid *));
    int *LocalGridIDs = malloc(NumLocalGrids*sizeof(int));
    char (*LocalGridNames)[OVK_NAME_LENGTH] = malloc(NumLocalGrids*sizeof(char[OVK_NAME_LENGTH]));
    MPI_Comm *LocalGridComms = malloc(NumLocalGrids*sizeof(MPI_Comm));
    int *LocalGridGlobalSizes[MAX_DIMS];
    LocalGridGlobalSizes[0] = malloc(MAX_DIMS*NumLocalGrids*sizeof(int));
    LocalGridGlobalSizes[1] = LocalGridGlobalSizes[0] + NumLocalGrids;
    LocalGridGlobalSizes[2] = LocalGridGlobalSizes[1] + NumLocalGrids;
    iLocalGrid = 0;
    for (iGrid = 0; iGrid < NumGrids; ++iGrid) {
      int GridID = iGrid+1;
      if (ovkRankHasGrid(Domain, GridID)) {
        const ovk_grid *Grid;
        ovkGetGrid(Domain, GridID, &Grid);
        const ovk_grid_properties *GridProperties;
        ovkGetGridProperties(Grid, &GridProperties);
        char Name[OVK_NAME_LENGTH];
        int GlobalSize[MAX_DIMS];
        MPI_Comm GridComm;
        ovkGetGridPropertyName(GridProperties, Name);
        ovkGetGridPropertySize(GridProperties, GlobalSize);
        ovkGetGridPropertyComm(GridProperties, &GridComm);
        LocalGrids[iLocalGrid] = Grid;
        LocalGridIDs[iLocalGrid] = GridID;
        strncpy(LocalGridNames[iLocalGrid], Name, OVK_NAME_LENGTH);
        LocalGridComms[iLocalGrid] = GridComm;
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          LocalGridGlobalSizes[iDim][iLocalGrid] = GlobalSize[iDim];
        }
        ++iLocalGrid;
      }
    }

    int CreateTime = GetProfilerTimerID(Profiler, "XINTOUT::Create");
    int DestroyTime = GetProfilerTimerID(Profiler, "XINTOUT::Destroy");
    int ReadTime = GetProfilerTimerID(Profiler, "XINTOUT::Read");
    int MatchTime = GetProfilerTimerID(Profiler, "XINTOUT::Match");
    int DistributeTime = GetProfilerTimerID(Profiler, "XINTOUT::Distribute");
    int ImportTime = GetProfilerTimerID(Profiler, "XINTOUT::Import");

    StartProfileSync(Profiler, CreateTime, Comm);

    t_xintout *XINTOUT;
    CreateXINTOUT(&XINTOUT, NumDims, Comm, NumGrids, NumLocalGrids, LocalGridIDs,
      (const char (*)[])LocalGridNames, LocalGridComms, (const int **)LocalGridGlobalSizes, Logger,
      ErrorHandler);

    EndProfile(Profiler, CreateTime);
    StartProfileSync(Profiler, ReadTime, Comm);

    Error = ReadXINTOUT(XINTOUT, HOPath, XPath, ReadGranularityAdjust, MPIInfo, Profiler);

    EndProfile(Profiler, ReadTime);
    OVK_EH_CHECK_GOTO(ErrorHandler, Error, destroy_xintout);
    StartProfileSync(Profiler, MatchTime, Comm);

    MatchDonorsAndReceivers(XINTOUT, Profiler);

    EndProfile(Profiler, MatchTime);
    StartProfileSync(Profiler, DistributeTime, Comm);

    t_donor_data **LocalDonors = malloc(NumLocalGrids*sizeof(t_donor_data *));
    t_receiver_data **LocalReceivers = malloc(NumLocalGrids*sizeof(t_receiver_data *));

    DistributeConnectivityData(XINTOUT, LocalGrids, LocalDonors, LocalReceivers, Profiler);

    EndProfile(Profiler, DistributeTime);
    StartProfileSync(Profiler, ImportTime, Comm);

    ImportConnectivityData(NumGrids, NumLocalGrids, LocalGridIDs, (const t_donor_data **)LocalDonors,
      (const t_receiver_data **)LocalReceivers, Domain);

    EndProfile(Profiler, ImportTime);
    StartProfileSync(Profiler, DestroyTime, Comm);

    for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
      DestroyDonorData(LocalDonors+iLocalGrid);
      DestroyReceiverData(LocalReceivers+iLocalGrid);
    }
    free(LocalDonors);
    free(LocalReceivers);
    free(LocalGridIDs);
    free(LocalGridNames);
    free(LocalGridComms);
    free(LocalGridGlobalSizes[0]);
    free(LocalGrids);

    EndProfile(Profiler, DestroyTime);

    destroy_xintout:
      StartProfileSync(Profiler, DestroyTime, Comm);
      DestroyXINTOUT(&XINTOUT);
      EndProfile(Profiler, DestroyTime);
      OVK_EH_CHECK_GOTO(ErrorHandler, Error, destroy_profiler);

  }

  WriteProfileTimes(Profiler, stdout);

  destroy_profiler:
    DestroyProfiler(&Profiler);

  OVK_EH_FINALIZE(ErrorHandler);

  return OVK_NO_ERROR;

}

ovk_error ovkEXTExportXINTOUT(const ovk_domain *Domain, const char *HOPath, const char *XPath,
  ovk_ext_xintout_format Format, ovk_ext_endian Endian, int WriteGranularityAdjust,
  MPI_Info MPIInfo) {

  OVK_DEBUG_ASSERT(false, "ovkEXTExportXINTOUT is not yet implemented.");

  return OVK_NO_ERROR;

}

static void CreateXINTOUT(t_xintout **XINTOUT_, int NumDims, MPI_Comm Comm, int NumGrids,
  int NumLocalGrids, const int *LocalGridIDs, const char (*LocalGridNames)[OVK_NAME_LENGTH],
  const MPI_Comm *LocalGridComms, const int **LocalGridGlobalSizes, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  int iLocalGrid;

  *XINTOUT_ = malloc(sizeof(t_xintout));
  t_xintout *XINTOUT = *XINTOUT_;

  MPI_Barrier(Comm);

  int CommSize, CommRank;
  MPI_Comm_size(Comm, &CommSize);
  MPI_Comm_rank(Comm, &CommRank);

  XINTOUT->num_dims = NumDims;
  XINTOUT->comm = Comm;
  XINTOUT->comm_size = CommSize;
  XINTOUT->comm_rank = CommRank;

  XINTOUT->logger = Logger;
  XINTOUT->error_handler = ErrorHandler;

  XINTOUT->num_grids = NumGrids;
  XINTOUT->num_local_grids = NumLocalGrids;

  XINTOUT->grids = malloc(NumLocalGrids*sizeof(t_xintout_grid *));
  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    int GlobalSize[MAX_DIMS] = {
      LocalGridGlobalSizes[0][iLocalGrid],
      LocalGridGlobalSizes[1][iLocalGrid],
      LocalGridGlobalSizes[2][iLocalGrid]
    };
    CreateXINTOUTGrid(XINTOUT->grids+iLocalGrid, LocalGridIDs[iLocalGrid],
      LocalGridNames[iLocalGrid], NumDims, LocalGridComms[iLocalGrid], GlobalSize, Logger,
      ErrorHandler);
  }

  XINTOUT->connections.count = 0;
  XINTOUT->connections.bin_size = 0;
  XINTOUT->connections.has_bin = false;
  XINTOUT->connections.bin = NULL;

  MPI_Barrier(Comm);

}

static void DestroyXINTOUT(t_xintout **XINTOUT_) {

  int iLocalGrid;

  t_xintout *XINTOUT = *XINTOUT_;

  MPI_Comm Comm = XINTOUT->comm;

  MPI_Barrier(Comm);

  if (XINTOUT->connections.has_bin) {
    t_xintout_connection_bin *Bin = XINTOUT->connections.bin;
    DestroyConnectionData(&Bin->data);
    free(Bin);
  }

  for (iLocalGrid = 0; iLocalGrid < XINTOUT->num_local_grids; ++iLocalGrid) {
    DestroyXINTOUTGrid(XINTOUT->grids+iLocalGrid);
  }
  free(XINTOUT->grids);

  free_null(XINTOUT_);

  MPI_Barrier(Comm);

}

static void CreateXINTOUTGrid(t_xintout_grid **XINTOUTGrid_, int ID, const char *Name, int NumDims,
  MPI_Comm Comm, const int *GlobalSize, t_logger *Logger, t_error_handler *ErrorHandler) {

  int iDim;

  *XINTOUTGrid_ = malloc(sizeof(t_xintout_grid));
  t_xintout_grid *XINTOUTGrid = *XINTOUTGrid_;

  MPI_Barrier(Comm);

  int CommSize, CommRank;
  MPI_Comm_size(Comm, &CommSize);
  MPI_Comm_rank(Comm, &CommRank);

  XINTOUTGrid->id = ID;

  strncpy(XINTOUTGrid->name, Name, OVK_NAME_LENGTH);

  XINTOUTGrid->num_dims = NumDims;
  XINTOUTGrid->comm = Comm;
  XINTOUTGrid->comm_size = CommSize;
  XINTOUTGrid->comm_rank = CommRank;

  XINTOUTGrid->logger = Logger;
  XINTOUTGrid->error_handler = ErrorHandler;

  for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
    XINTOUTGrid->global_size[iDim] = GlobalSize[iDim];
  }

  XINTOUTGrid->donors.count = 0;
  XINTOUTGrid->donors.chunk_size = 0;
  XINTOUTGrid->donors.has_chunk = false;
  XINTOUTGrid->donors.chunk = NULL;

  XINTOUTGrid->receivers.count = 0;
  XINTOUTGrid->receivers.chunk_size = 0;
  XINTOUTGrid->receivers.has_chunk = false;
  XINTOUTGrid->receivers.chunk = NULL;

}

static void DestroyXINTOUTGrid(t_xintout_grid **XINTOUTGrid_) {

  t_xintout_grid *XINTOUTGrid = *XINTOUTGrid_;

  MPI_Comm Comm = XINTOUTGrid->comm;

  MPI_Barrier(Comm);

  if (XINTOUTGrid->donors.has_chunk) {
    t_xintout_donor_chunk *Chunk = XINTOUTGrid->donors.chunk;
    DestroyDonorData(&Chunk->data);
    free(Chunk);
  }

  if (XINTOUTGrid->receivers.has_chunk) {
    t_xintout_receiver_chunk *Chunk = XINTOUTGrid->receivers.chunk;
    DestroyReceiverData(&Chunk->data);
    free(Chunk->connection_ids);
    free(Chunk);
  }

  free_null(XINTOUTGrid_);

  MPI_Barrier(Comm);

}

static ovk_error ReadXINTOUT(t_xintout *XINTOUT, const char *HOPath, const char *XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo, t_profiler *Profiler) {

  int iLocalGrid;

  MPI_Comm Comm = XINTOUT->comm;
  int CommRank = XINTOUT->comm_rank;

  MPI_Barrier(Comm);

  t_error_handler *ErrorHandler = XINTOUT->error_handler;
  t_logger *Logger = XINTOUT->logger;

  int NumLocalGrids = XINTOUT->num_local_grids;

  OVK_EH_INIT(ErrorHandler);
  ovk_error Error = OVK_NO_ERROR;

  LogStatus(Logger, CommRank == 0, 0, "Reading XINTOUT files '%s' and '%s'...", HOPath, XPath);

  ovk_ext_endian Endian;
  ovk_ext_xintout_format Format;
  bool WithIBlank;
  Error = ReadGlobalInfo(XINTOUT, HOPath, &Endian, &Format, &WithIBlank, Profiler);
  OVK_EH_CHECK(ErrorHandler, Error);

  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {

    t_xintout_grid *XINTOUTGrid = XINTOUT->grids[iLocalGrid];

    size_t NumDonors, NumReceivers;
    size_t DonorStartingConnectionID;
    MPI_Offset HODonorCellsOffset, HODonorCoordsOffset, XDonorSizesOffset, XDonorInterpCoefsOffset;
    MPI_Offset HOReceiverPointsOffset, HOReceiverConnectionIDsOffset;
    Error = ReadGridInfo(XINTOUTGrid, HOPath, XPath, &NumDonors, &NumReceivers,
      &DonorStartingConnectionID, &HODonorCellsOffset, &HODonorCoordsOffset, &HOReceiverPointsOffset,
      &HOReceiverConnectionIDsOffset, &XDonorSizesOffset, &XDonorInterpCoefsOffset, Endian,
      Format, WithIBlank, Profiler);
    OVK_EH_CHECK_GOTO(ErrorHandler, Error, done_reading);

    Error = ReadDonors(XINTOUTGrid, HOPath, XPath, NumDonors, DonorStartingConnectionID,
      HODonorCellsOffset, HODonorCoordsOffset, XDonorSizesOffset, XDonorInterpCoefsOffset, Endian,
      ReadGranularityAdjust, MPIInfo, Profiler);
    OVK_EH_CHECK_GOTO(ErrorHandler, Error, done_reading);

    Error = ReadReceivers(XINTOUTGrid, HOPath, NumReceivers, HOReceiverPointsOffset,
      HOReceiverConnectionIDsOffset, Endian, Format, ReadGranularityAdjust, MPIInfo, Profiler);
    OVK_EH_CHECK_GOTO(ErrorHandler, Error, done_reading);

  }

  done_reading:
    OVK_EH_CHECK_ALL(ErrorHandler, Error, Comm);

  MPI_Barrier(Comm);

  LogStatus(Logger, CommRank == 0, 0, "Finished reading XINTOUT files.");

  OVK_EH_FINALIZE(ErrorHandler);

  return OVK_NO_ERROR;

}

static ovk_error ReadGlobalInfo(const t_xintout *XINTOUT, const char *HOPath,
  ovk_ext_endian *Endian_, ovk_ext_xintout_format *Format_, bool *WithIBlank_,
  t_profiler *Profiler) {

  int iDim;

  MPI_Comm Comm = XINTOUT->comm;
  int CommRank = XINTOUT->comm_rank;

  int NumGrids = XINTOUT->num_grids;

  t_logger *Logger = (t_logger *)XINTOUT->logger;
  t_error_handler *ErrorHandler = (t_error_handler *)XINTOUT->error_handler;

  OVK_EH_INIT(ErrorHandler);
  ovk_error Error = OVK_NO_ERROR;

  ovk_ext_endian Endian;
  ovk_ext_xintout_format Format;
  bool WithIBlank;

  int MPIIOOpenTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Open");
  int MPIIOCloseTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Close");

  if (CommRank == 0) {

    MPI_File HOFile;
    MPI_Status Status;
    int MPIError;
    int ReadSize;

    StartProfile(Profiler, MPIIOOpenTime);
    MPIError = MPI_File_open(MPI_COMM_SELF, (char *)HOPath, MPI_MODE_RDONLY, MPI_INFO_NULL,
      &HOFile);
    EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      LogError(Logger, true, "Unable to open file '%s'.", HOPath);
      Error = OVK_ERROR_FILE_OPEN;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, done_reading);
    }

    bool Success = DetectFormat(HOFile, &Endian, &Format, Profiler);
    if (!Success) {
      LogError(Logger, true, "Unable to detect format of XINTOUT file '%s'.", HOPath);
      Error = OVK_ERROR_FILE_READ;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
    }

    int RecordWrapperSize = Format == OVK_EXT_XINTOUT_STANDARD ? sizeof(int) : 0;

    // Both formats have a record wrapper at the beginning
    MPI_Offset HOHeaderOffset = sizeof(int);
    int NumGridsInFile;
    File_read_at_endian(HOFile, HOHeaderOffset, &NumGridsInFile, 1, MPI_INT, Endian, &Status,
      Profiler);
    MPI_Get_count(&Status, MPI_INT, &ReadSize);
    if (ReadSize < 1) {
      LogError(Logger, true, "Unable to read header of XINTOUT file '%s'.", HOPath);
      Error = OVK_ERROR_FILE_READ;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
    }
    if (NumGridsInFile != NumGrids) {
      LogError(Logger, CommRank == 0, "XINTOUT file '%s' has incorrect number of grids.", HOPath);
      Error = OVK_ERROR_FILE_READ;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
    }

    MPI_Offset HOGridOffset = 0;

    HOGridOffset += RecordWrapperSize;
    if (Format == OVK_EXT_XINTOUT_STANDARD) {
      HOGridOffset += 5*sizeof(int);
    } else {
      HOGridOffset += sizeof(int) + 4*sizeof(long long);
    }
    HOGridOffset += RecordWrapperSize;

    HOGridOffset += RecordWrapperSize;
    size_t NumDonors;
    size_t NumReceivers;
    int GridSize[MAX_DIMS];
    if (Format == OVK_EXT_XINTOUT_STANDARD) {
      int Data[7];
      File_read_at_endian(HOFile, HOGridOffset, Data, 7, MPI_INT, Endian, &Status, Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < 7) {
        LogError(Logger, true, "Unable to read grid 1 header of XINTOUT file '%s'.", HOPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
      }
      NumDonors = Data[1];
      NumReceivers = Data[0];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        GridSize[iDim] = Data[4+iDim];
      }
      HOGridOffset += 7*sizeof(int);
    } else {
      long long LongLongData[4];
      int IntData[3];
      File_read_at_endian(HOFile, HOGridOffset, LongLongData, 4, MPI_LONG_LONG, Endian, &Status,
        Profiler);
      MPI_Get_count(&Status, MPI_LONG_LONG, &ReadSize);
      if (ReadSize < 4) {
        LogError(Logger, true, "Unable to read grid 1 header of XINTOUT file '%s'.", HOPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
      }
      HOGridOffset += 4*sizeof(long long);
      File_read_at_endian(HOFile, HOGridOffset, IntData, 3, MPI_INT, Endian, &Status, Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < 3) {
        LogError(Logger, true, "Unable to read grid 1 header of XINTOUT file '%s'.", HOPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
      }
      HOGridOffset += 3*sizeof(int);
      NumDonors = LongLongData[1];
      NumReceivers = LongLongData[0];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        GridSize[iDim] = IntData[iDim];
      }
      HOGridOffset += sizeof(long long);
    }
    HOGridOffset += RecordWrapperSize;

    HOGridOffset += RecordWrapperSize;
    HOGridOffset += MAX_DIMS*NumDonors*sizeof(int);
    HOGridOffset += MAX_DIMS*NumDonors*sizeof(double);
    HOGridOffset += RecordWrapperSize;

    int ConnectionIDSize = Format == OVK_EXT_XINTOUT_STANDARD ? sizeof(int) : sizeof(long long);
    HOGridOffset += RecordWrapperSize;
    HOGridOffset += MAX_DIMS*NumReceivers*sizeof(int);
    HOGridOffset += NumReceivers*ConnectionIDSize;
    HOGridOffset += RecordWrapperSize;

    if (Format == OVK_EXT_XINTOUT_STANDARD) {
      int RecordWrapper;
      File_read_at_endian(HOFile, HOGridOffset, &RecordWrapper, 1, MPI_INT, Endian, &Status,
        Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize == 1) {
        size_t NumPoints = (size_t)GridSize[0]*(size_t)GridSize[1]*(size_t)GridSize[2];
        if (RecordWrapper == NumPoints*sizeof(int)) {
          WithIBlank = true;
        } else {
          WithIBlank = false;
        }
      } else {
        if (NumGrids == 1) {
          WithIBlank = false;
        } else {
          LogError(Logger, true, "Unable to detect whether XINTOUT file '%s' contains IBlank.",
            HOPath);
          Error = OVK_ERROR_FILE_READ;
          OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
        }
      }
    } else {
      WithIBlank = false;
    }

    close_ho:
      StartProfile(Profiler, MPIIOCloseTime);
      MPI_File_close(&HOFile);
      EndProfile(Profiler, MPIIOCloseTime);

  }

  done_reading:
    OVK_EH_CHECK_ALL(ErrorHandler, Error, Comm);

  int EndianInt;
  if (CommRank == 0) EndianInt = Endian == OVK_EXT_BIG_ENDIAN ? 1 : 0;
  MPI_Bcast(&EndianInt, 1, MPI_INT, 0, Comm);
  Endian = EndianInt == 1 ? OVK_EXT_BIG_ENDIAN : OVK_EXT_LITTLE_ENDIAN;

  int FormatInt;
  if (CommRank == 0) FormatInt = Format == OVK_EXT_XINTOUT_EXTENDED ? 1 : 0;
  MPI_Bcast(&FormatInt, 1, MPI_INT, 0, Comm);
  Format = FormatInt == 1 ? OVK_EXT_XINTOUT_EXTENDED : OVK_EXT_XINTOUT_STANDARD;

  int WithIBlankInt;
  if (CommRank == 0) WithIBlankInt = (int)WithIBlank;
  MPI_Bcast(&WithIBlankInt, 1, MPI_INT, 0, Comm);
  WithIBlank = WithIBlankInt != 0;

  *Endian_ = Endian;
  *Format_ = Format;
  *WithIBlank_ = WithIBlank;

  OVK_EH_FINALIZE(ErrorHandler);

  return OVK_NO_ERROR;

}

static bool DetectFormat(MPI_File HOFile, ovk_ext_endian *Endian_, ovk_ext_xintout_format *Format_,
  t_profiler *Profiler) {

  ovk_ext_endian Endian;
  ovk_ext_xintout_format Format;

  int MPIIOReadTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Read");

  unsigned char InitialBytes[sizeof(int)];
  MPI_Status Status;
  StartProfile(Profiler, MPIIOReadTime);
  int MPIError = MPI_File_read_at(HOFile, 0, InitialBytes, sizeof(int), MPI_BYTE, &Status);
  EndProfile(Profiler, MPIIOReadTime);
  if (MPIError != MPI_SUCCESS) return false;

  // If little endian, the first byte will be the size of the file header data
  if (InitialBytes[0] != 0) {
    Endian = OVK_EXT_LITTLE_ENDIAN;
  } else {
    Endian = OVK_EXT_BIG_ENDIAN;
  }

  int HeaderSize;
  memcpy(&HeaderSize, InitialBytes, sizeof(int));
  if (Endian != MachineEndian()) {
    SwapEndian(&HeaderSize, sizeof(int), 1);
  }

  if (HeaderSize == 5*sizeof(int)) {
    Format = OVK_EXT_XINTOUT_STANDARD;
  } else if (HeaderSize == sizeof(int) + 4*sizeof(long long)) {
    Format = OVK_EXT_XINTOUT_EXTENDED;
  } else {
    return false;
  }

  *Endian_ = Endian;
  *Format_ = Format;

  return true;

}

static ovk_error ReadGridInfo(const t_xintout_grid *XINTOUTGrid, const char *HOPath,
  const char *XPath, size_t *NumDonors_, size_t *NumReceivers_, size_t *StartingConnectionID_,
  MPI_Offset *HODonorCellsOffset_, MPI_Offset *HODonorCoordsOffset_,
  MPI_Offset *HOReceiverPointsOffset_, MPI_Offset *HOReceiverConnectionIDsOffset_,
  MPI_Offset *XDonorSizesOffset_, MPI_Offset *XDonorInterpCoefsOffset_, ovk_ext_endian Endian,
  ovk_ext_xintout_format Format, bool WithIBlank, t_profiler *Profiler) {

  int iDim;

  int GridID = XINTOUTGrid->id;
  int iGrid = GridID-1;

  MPI_Comm Comm = XINTOUTGrid->comm;
  int CommRank = XINTOUTGrid->comm_rank;

  t_logger *Logger = (t_logger *)XINTOUTGrid->logger;
  t_error_handler *ErrorHandler = (t_error_handler *)XINTOUTGrid->error_handler;

  OVK_EH_INIT(ErrorHandler);
  ovk_error Error = OVK_NO_ERROR;

  size_t NumDonors, NumReceivers;
  size_t StartingConnectionID;
  MPI_Offset HODonorCellsOffset, HODonorCoordsOffset, HOReceiverPointsOffset,
    HOReceiverConnectionIDsOffset, XDonorSizesOffset, XDonorInterpCoefsOffset;

  int MPIIOOpenTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Open");
  int MPIIOCloseTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Close");

  if (CommRank == 0) {

    int iOtherGrid;

    MPI_File HOFile, XFile;
    MPI_Status Status;
    int MPIError;
    int ReadSize;
    int GridSize[MAX_DIMS];
    size_t NumInterpCoefs;

    StartProfile(Profiler, MPIIOOpenTime);
    MPIError = MPI_File_open(MPI_COMM_SELF, (char *)HOPath, MPI_MODE_RDONLY, MPI_INFO_NULL,
      &HOFile);
    EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      LogError(Logger, true, "Unable to open file '%s'.", HOPath);
      Error = OVK_ERROR_FILE_OPEN;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, done_reading);
    }

    StartProfile(Profiler, MPIIOOpenTime);
    MPIError = MPI_File_open(MPI_COMM_SELF, (char *)XPath, MPI_MODE_RDONLY, MPI_INFO_NULL,
      &XFile);
    EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      LogError(Logger, true, "Unable to open file '%s'.", XPath);
      Error = OVK_ERROR_FILE_OPEN;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
    }

    int RecordWrapperSize = Format == OVK_EXT_XINTOUT_STANDARD ? sizeof(int) : 0;

    MPI_Offset HOGridOffset = 0;
    MPI_Offset XGridOffset = 0;

    int HeaderSize;
    if (Format == OVK_EXT_XINTOUT_STANDARD) {
      HeaderSize = 5*sizeof(int);
    } else {
      HeaderSize = sizeof(int) + 4*sizeof(long long);
    }
    // Both formats have a record wrapper at the beginning
    HOGridOffset += sizeof(int);
    HOGridOffset += HeaderSize;
    HOGridOffset += RecordWrapperSize;
    // Both formats have a record wrapper at the beginning
    XGridOffset += sizeof(int);
    XGridOffset += HeaderSize;
    XGridOffset += RecordWrapperSize;

    for (iOtherGrid = 0; iOtherGrid < iGrid; ++iOtherGrid) {

      int OtherGridID = iOtherGrid+1;

      HOGridOffset += RecordWrapperSize;
      if (Format == OVK_EXT_XINTOUT_STANDARD) {
        int Data[7];
        File_read_at_endian(HOFile, HOGridOffset, Data, 7, MPI_INT, Endian, &Status, Profiler);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        if (ReadSize < 7) {
          LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
            HOPath);
          Error = OVK_ERROR_FILE_READ;
          OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
        }
        NumDonors = Data[1];
        NumReceivers = Data[0];
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          GridSize[iDim] = Data[4+iDim];
        }
        HOGridOffset += 7*sizeof(int);
      } else {
        long long LongLongData[5];
        int IntData[3];
        File_read_at_endian(HOFile, HOGridOffset, LongLongData, 4, MPI_LONG_LONG, Endian, &Status,
          Profiler);
        MPI_Get_count(&Status, MPI_LONG_LONG, &ReadSize);
        if (ReadSize < 4) {
          LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
            HOPath);
          Error = OVK_ERROR_FILE_READ;
          OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
        }
        HOGridOffset += 4*sizeof(long long);
        File_read_at_endian(HOFile, HOGridOffset, IntData, 3, MPI_INT, Endian, &Status, Profiler);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        if (ReadSize < 3) {
          LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
            HOPath);
          Error = OVK_ERROR_FILE_READ;
          OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
        }
        HOGridOffset += 3*sizeof(int);
        File_read_at_endian(HOFile, HOGridOffset, LongLongData+4, 1, MPI_LONG_LONG, Endian, &Status,
          Profiler);
        MPI_Get_count(&Status, MPI_LONG_LONG, &ReadSize);
        if (ReadSize < 1) {
          LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
            HOPath);
          Error = OVK_ERROR_FILE_READ;
          OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
        }
        HOGridOffset += sizeof(long long);
        NumDonors = LongLongData[1];
        NumReceivers = LongLongData[0];
        NumInterpCoefs = LongLongData[4];
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          GridSize[iDim] = IntData[iDim];
        }
      }
      HOGridOffset += RecordWrapperSize;

      HOGridOffset += RecordWrapperSize;
      HOGridOffset += MAX_DIMS*NumDonors*sizeof(int);
      HOGridOffset += MAX_DIMS*NumDonors*sizeof(double);
      HOGridOffset += RecordWrapperSize;

      int ConnectionIDSize = Format == OVK_EXT_XINTOUT_STANDARD ? sizeof(int) : sizeof(long long);
      HOGridOffset += RecordWrapperSize;
      HOGridOffset += MAX_DIMS*NumReceivers*sizeof(int);
      HOGridOffset += NumReceivers*ConnectionIDSize;
      HOGridOffset += RecordWrapperSize;

      if (WithIBlank) {
        HOGridOffset += RecordWrapperSize;
        HOGridOffset += (size_t)GridSize[0]*(size_t)GridSize[1]*(size_t)GridSize[2]*sizeof(int);
        HOGridOffset += RecordWrapperSize;
      }

      XGridOffset += RecordWrapperSize;
      XGridOffset += MAX_DIMS*NumDonors*sizeof(int);
      XGridOffset += RecordWrapperSize;

      if (Format == OVK_EXT_XINTOUT_STANDARD) {
        // Figure out size of interp coef data
        int RecordWrapper;
        File_read_at_endian(XFile, XGridOffset, &RecordWrapper, 1, MPI_INT, Endian, &Status,
          Profiler);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        if (ReadSize < 1) {
          LogError(Logger, true, "Unable to read grid %i data from XINTOUT file '%s'.",
            OtherGridID, XPath);
          Error = OVK_ERROR_FILE_READ;
          OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_x);
        }
        NumInterpCoefs = RecordWrapper/sizeof(double);
      }

      XGridOffset += RecordWrapperSize;
      XGridOffset += NumInterpCoefs*sizeof(double);
      XGridOffset += RecordWrapperSize;

      XGridOffset += RecordWrapperSize;
      XGridOffset += 3*sizeof(int);
      XGridOffset += RecordWrapperSize;

    }

    HOGridOffset += RecordWrapperSize;
    if (Format == OVK_EXT_XINTOUT_STANDARD) {
      int Data[7];
      File_read_at_endian(HOFile, HOGridOffset, Data, 7, MPI_INT, Endian, &Status, Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < 7) {
        LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
          HOPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
      }
      NumDonors = Data[1];
      NumReceivers = Data[0];
      // Convert to zero-based indexing
      StartingConnectionID = Data[3]-1;
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        GridSize[iDim] = Data[4+iDim];
      }
      HOGridOffset += 7*sizeof(int);
    } else {
      long long LongLongData[5];
      int IntData[3];
      File_read_at_endian(HOFile, HOGridOffset, LongLongData, 4, MPI_LONG_LONG, Endian, &Status,
        Profiler);
      MPI_Get_count(&Status, MPI_LONG_LONG, &ReadSize);
      if (ReadSize < 4) {
        LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
          HOPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
      }
      HOGridOffset += 4*sizeof(long long);
      File_read_at_endian(HOFile, HOGridOffset, IntData, 3, MPI_INT, Endian, &Status, Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < 3) {
        LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
          HOPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
      }
      HOGridOffset += 3*sizeof(int);
      File_read_at_endian(HOFile, HOGridOffset, LongLongData+4, 1, MPI_LONG_LONG, Endian, &Status,
        Profiler);
      MPI_Get_count(&Status, MPI_LONG_LONG, &ReadSize);
      if (ReadSize < 1) {
        LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
          HOPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
      }
      HOGridOffset += sizeof(long long);
      NumDonors = LongLongData[1];
      NumReceivers = LongLongData[0];
      // Convert to zero-based indexing
      StartingConnectionID = LongLongData[3]-1;
      NumInterpCoefs = LongLongData[4];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        GridSize[iDim] = IntData[iDim];
      }
    }
    HOGridOffset += RecordWrapperSize;

    if (GridSize[0] != XINTOUTGrid->global_size[0] || GridSize[1] != XINTOUTGrid->global_size[1] ||
      GridSize[2] != XINTOUTGrid->global_size[2]) {
      LogError(Logger, true, "Grid %i of XINTOUT file '%s' has incorrect size.", GridID, HOPath);
      Error = OVK_ERROR_FILE_READ;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_x);
    }

    HOGridOffset += RecordWrapperSize;
    HODonorCellsOffset = HOGridOffset;
    HOGridOffset += MAX_DIMS*NumDonors*sizeof(int);
    HODonorCoordsOffset = HOGridOffset;
    HOGridOffset += MAX_DIMS*NumDonors*sizeof(double);
    HOGridOffset += RecordWrapperSize;

    HOGridOffset += RecordWrapperSize;
    HOReceiverPointsOffset = HOGridOffset;
    HOGridOffset += MAX_DIMS*NumReceivers*sizeof(int);
    HOReceiverConnectionIDsOffset = HOGridOffset;

    XGridOffset += RecordWrapperSize;
    XDonorSizesOffset = XGridOffset;
    XGridOffset += MAX_DIMS*NumDonors*sizeof(int);
    XGridOffset += RecordWrapperSize;
    XGridOffset += RecordWrapperSize;
    XDonorInterpCoefsOffset = XGridOffset;

    close_x:
      StartProfile(Profiler, MPIIOCloseTime);
      MPI_File_close(&XFile);
      EndProfile(Profiler, MPIIOCloseTime);

    close_ho:
      StartProfile(Profiler, MPIIOCloseTime);
      MPI_File_close(&HOFile);
      EndProfile(Profiler, MPIIOCloseTime);

  }

  done_reading:
    OVK_EH_CHECK_ALL(ErrorHandler, Error, Comm);

  MPI_Bcast(&NumDonors, 1, KMPI_UNSIGNED_SIZE, 0, Comm);
  MPI_Bcast(&NumReceivers, 1, KMPI_UNSIGNED_SIZE, 0, Comm);
  MPI_Bcast(&StartingConnectionID, 1, KMPI_UNSIGNED_SIZE, 0, Comm);
  MPI_Bcast(&HODonorCellsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&HODonorCoordsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&HOReceiverPointsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&HOReceiverConnectionIDsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&XDonorSizesOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&XDonorInterpCoefsOffset, 1, MPI_OFFSET, 0, Comm);

  *NumDonors_ = NumDonors;
  *NumReceivers_ = NumReceivers;
  *StartingConnectionID_ = StartingConnectionID;
  *HODonorCellsOffset_ = HODonorCellsOffset;
  *HODonorCoordsOffset_ = HODonorCoordsOffset;
  *HOReceiverPointsOffset_ = HOReceiverPointsOffset;
  *HOReceiverConnectionIDsOffset_ = HOReceiverConnectionIDsOffset;
  *XDonorSizesOffset_ = XDonorSizesOffset;
  *XDonorInterpCoefsOffset_ = XDonorInterpCoefsOffset;

  OVK_EH_FINALIZE(ErrorHandler);

  return OVK_NO_ERROR;

}

static ovk_error ReadDonors(t_xintout_grid *XINTOUTGrid, const char *HOPath, const char *XPath,
  size_t NumDonors, size_t StartingConnectionID, MPI_Offset HOCellsOffset, MPI_Offset HOCoordsOffset,
  MPI_Offset XSizesOffset, MPI_Offset XInterpCoefsOffset, ovk_ext_endian Endian,
  int ReadGranularityAdjust, MPI_Info MPIInfo, t_profiler *Profiler) {

  int iDim;
  int iPoint;
  size_t iDonor;

  int GridID = XINTOUTGrid->id;

  MPI_Comm Comm = XINTOUTGrid->comm;
  int CommSize = XINTOUTGrid->comm_size;
  int CommRank = XINTOUTGrid->comm_rank;

  t_logger *Logger = XINTOUTGrid->logger;
  t_error_handler *ErrorHandler = XINTOUTGrid->error_handler;

  OVK_EH_INIT(ErrorHandler);
  ovk_error Error = OVK_NO_ERROR;

  t_xintout_donors *XINTOUTDonors = &XINTOUTGrid->donors;

  XINTOUTDonors->count = NumDonors;

  // Read chunk on every Nth rank, where N is a power of 2 (or CommSize) chosen such that the number
  // of donors per chunk is roughly the average number of grid points per rank (subject to user
  // adjustment)
  size_t NumPoints =
    (size_t)XINTOUTGrid->global_size[0]*
    (size_t)XINTOUTGrid->global_size[1]*
    (size_t)XINTOUTGrid->global_size[2];
  size_t AvgPointsPerRank = BinDivide(NumPoints, CommSize);
  int ChunkRankInterval, NumChunks;
  size_t ChunkSize;
  Chunkify(NumDonors, CommSize, AvgPointsPerRank, ReadGranularityAdjust, &ChunkRankInterval,
    &NumChunks, &ChunkSize);
  bool HasChunk = CommRank % ChunkRankInterval == 0;

  if (LoggingStatus(Logger)) {
    char NumDonorsString[NUMBER_STRING_LENGTH+6];
    char NumRanksString[NUMBER_STRING_LENGTH+9];
    PluralizeLabel(NumDonors, "donors", "donor", NumDonorsString);
    PluralizeLabel(NumChunks, "I/O ranks", "I/O rank", NumRanksString);
    LogStatus(Logger, CommRank == 0, 1, "Grid %s has %s; using %s.",
      XINTOUTGrid->name, NumDonorsString, NumRanksString);
  }

  XINTOUTDonors->chunk_size = ChunkSize;
  XINTOUTDonors->has_chunk = HasChunk;

  MPI_Comm ChunkComm;
  MPI_Comm_split(Comm, HasChunk ? 0 : MPI_UNDEFINED, CommRank, &ChunkComm);

  t_xintout_donor_chunk *Chunk = NULL;

  if (HasChunk) {

    Chunk = malloc(sizeof(t_xintout_donor_chunk));
    XINTOUTDonors->chunk = Chunk;

    int ChunkCommRank;
    MPI_Comm_rank(ChunkComm, &ChunkCommRank);

    size_t LocalBegin = ChunkSize*ChunkCommRank;
    size_t LocalEnd = min(ChunkSize*(ChunkCommRank+1), NumDonors);
    size_t NumLocalDonors = LocalEnd - LocalBegin;

    if (NumLocalDonors*sizeof(double) > INT_MAX) {
      LogError(Logger, true, "Donor chunk size too big; increase number of processes or read "
        "granularity.");
      Error = OVK_ERROR_FILE_READ;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, free_comm);
    }

    Chunk->begin = LocalBegin;
    Chunk->end = LocalEnd;
    Chunk->starting_connection_id = StartingConnectionID + LocalBegin;

    MPI_File HOFile, XFile;
    MPI_Status Status;
    int MPIError;
    MPI_Offset DatasetOffset, ReadOffset;
    int ReadSize;

    int MPIIOOpenTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Open");
    int MPIIOCloseTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Close");
    int MPIIOOtherTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Other");

    StartProfileSync(Profiler, MPIIOOpenTime, ChunkComm);
    MPIError = MPI_File_open(ChunkComm, (char *)HOPath, MPI_MODE_RDONLY, MPIInfo, &HOFile);
    EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      LogError(Logger, true, "Unable to open file '%s'.", HOPath);
      Error = OVK_ERROR_FILE_OPEN;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, free_comm);
    }

    StartProfileSync(Profiler, MPIIOOpenTime, ChunkComm);
    MPIError = MPI_File_open(ChunkComm, (char *)XPath, MPI_MODE_RDONLY, MPIInfo, &XFile);
    EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      LogError(Logger, true, "Unable to open file '%s'.", XPath);
      Error = OVK_ERROR_FILE_OPEN;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
    }

    int *Sizes[MAX_DIMS];
    Sizes[0] = malloc(MAX_DIMS*NumLocalDonors*sizeof(int));
    Sizes[1] = Sizes[0] + NumLocalDonors;
    Sizes[2] = Sizes[1] + NumLocalDonors;

//     MPI_Datatype XDonorSizeType;
//     MPI_Type_vector(MAX_DIMS, NumLocalDonors, NumDonors, MPI_INT, &XDonorSizeType);
//     MPI_Type_commit(&XDonorSizeType);
//     MPI_File_set_view(XFile, XSizesOffset+LocalBegin*sizeof(int), MPI_INT, XDonorSizeType, "native",
//       MPIInfo);
//     int DataSize = MAX_DIMS*NumLocalDonors;
//     File_read_all_endian(XFile, Sizes[0], DataSize, MPI_INT, Endian, &Status);
//     MPI_Type_free(&XDonorSizeType);
//     MPI_Get_count(&Status, MPI_INT, &ReadSize);
//     if (ReadSize < DataSize) {
//       LogError(Logger, true, "Unable to read grid %i donor sizes from XINTOUT file '%s'.", GridID,
//         XPath);
//       Error = OVK_ERROR_FILE_READ;
//       OVK_EH_HANDLE_GOTO(ErrorHandler, Error, free_sizes);
//     }

    // The above version causes "Conditional jump or move depends on uninitialized value(s)"
    // errors in valgrind with MPICH 3.2... not sure why
    DatasetOffset = XSizesOffset;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + LocalBegin*sizeof(int);
      StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(XFile, ReadOffset, MPI_INT, MPI_INT, "native", MPIInfo);
      EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(XFile, Sizes[iDim], (int)NumLocalDonors, MPI_INT, Endian,
        &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < NumLocalDonors) {
        LogError(Logger, true, "Unable to read grid %i donor sizes from XINTOUT file '%s'.", GridID,
          XPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, free_sizes);
      }
      DatasetOffset += NumDonors*sizeof(int);
    }

    int MaxSize = 0;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        MaxSize = max(MaxSize, Sizes[iDim][iDonor]);
      }
    }

    CreateDonorData(&Chunk->data, NumLocalDonors, MaxSize);
    t_donor_data *Data = Chunk->data;

//     MPI_Datatype HODonorCellType;
//     MPI_Type_vector(MAX_DIMS, NumLocalDonors, NumDonors, MPI_INT, &HODonorCellType);
//     MPI_Type_commit(&HODonorCellType);
//     MPI_File_set_view(HOFile, HOCellsOffset+LocalBegin*sizeof(int), MPI_INT, HODonorCellType,
//       "native", MPIInfo);
//     DataSize = MAX_DIMS*NumLocalDonors;
//     File_read_all_endian(HOFile, Data->extents[0][0], DataSize, MPI_INT, Endian, &Status);
//     MPI_Type_free(&HODonorCellType);
//     MPI_Get_count(&Status, MPI_INT, &ReadSize);
//     if (ReadSize < DataSize) {
//       LogError(Logger, true, "Unable to read grid %i donor cells from XINTOUT file '%s'.", GridID,
//         HOPath);
//       Error = OVK_ERROR_FILE_READ;
//       OVK_EH_HANDLE_GOTO(ErrorHandler, Error, free_sizes);
//     }

    // The above version causes "Conditional jump or move depends on uninitialized value(s)"
    // errors in valgrind with MPICH 3.2... not sure why
    DatasetOffset = HOCellsOffset;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + LocalBegin*sizeof(int);
      StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(HOFile, ReadOffset, MPI_INT, MPI_INT, "native", MPIInfo);
      EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(HOFile, Data->extents[0][iDim], (int)NumLocalDonors, MPI_INT, Endian,
        &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < NumLocalDonors) {
        LogError(Logger, true, "Unable to read grid %i donor cells from XINTOUT file '%s'.", GridID,
          HOPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, free_sizes);
      }
      DatasetOffset += NumDonors*sizeof(int);
    }

    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        // Convert to zero-based indexing
        --Data->extents[0][iDim][iDonor];
        Data->extents[1][iDim][iDonor] = Data->extents[0][iDim][iDonor] + Sizes[iDim][iDonor];
      }
    }

    free_sizes:
      free(Sizes[0]);
      OVK_EH_CHECK_GOTO(ErrorHandler, Error, close_x);

//     MPI_Datatype HODonorCoordsType;
//     MPI_Type_vector(MAX_DIMS, NumLocalDonors, NumDonors, MPI_DOUBLE, &HODonorCoordsType);
//     MPI_Type_commit(&HODonorCoordsType);
//     MPI_File_set_view(HOFile, HOCoordsOffset+LocalBegin*sizeof(double), MPI_DOUBLE,
//       HODonorCoordsType, "native", MPIInfo);
//     DataSize = MAX_DIMS*NumLocalDonors;
//     File_read_all_endian(HOFile, Data->coords[0], DataSize, MPI_DOUBLE, Endian, &Status);
//     MPI_Type_free(&HODonorCoordsType);
//     MPI_Get_count(&Status, MPI_DOUBLE, &ReadSize);
//     if (ReadSize < DataSize) {
//       LogError(Logger, true, "Unable to read grid %i donor coords from XINTOUT file '%s'.", GridID,
//         HOPath);
//       Error = OVK_ERROR_FILE_READ;
//       OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_x);
//     }

    // The above version causes "Conditional jump or move depends on uninitialized value(s)"
    // errors in valgrind with MPICH 3.2... not sure why
    DatasetOffset = HOCoordsOffset;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + LocalBegin*sizeof(double);
      StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(HOFile, ReadOffset, MPI_DOUBLE, MPI_DOUBLE, "native", MPIInfo);
      EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(HOFile, Data->coords[iDim], (int)NumLocalDonors, MPI_DOUBLE, Endian,
        &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_DOUBLE, &ReadSize);
      if (ReadSize < NumLocalDonors) {
        LogError(Logger, true, "Unable to read grid %i donor coords from XINTOUT file '%s'.",
          GridID, HOPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_x);
      }
      DatasetOffset += NumDonors*sizeof(double);
    }

    size_t NumLocalInterpCoefs[MAX_DIMS] = {0, 0, 0};
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        int Size = Data->extents[1][iDim][iDonor] - Data->extents[0][iDim][iDonor];
        NumLocalInterpCoefs[iDim] += Size;
      }
    }

    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      if (NumLocalInterpCoefs[iDim]*sizeof(double) > INT_MAX) {
        LogError(Logger, true, "Donor chunk size too big; increase number of processes or read "
          "granularity.");
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_x);
      }
    }

    size_t NumInterpCoefs[MAX_DIMS];
    size_t NumInterpCoefsBeforeChunk[MAX_DIMS];
    MPI_Allreduce(&NumLocalInterpCoefs, &NumInterpCoefs, MAX_DIMS, KMPI_UNSIGNED_SIZE, MPI_SUM,
      ChunkComm);
    MPI_Scan(&NumLocalInterpCoefs, &NumInterpCoefsBeforeChunk, MAX_DIMS, KMPI_UNSIGNED_SIZE, MPI_SUM,
      ChunkComm);
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      NumInterpCoefsBeforeChunk[iDim] -= NumLocalInterpCoefs[iDim];
    }

    double *InterpCoefs = malloc((NumLocalInterpCoefs[0]+NumLocalInterpCoefs[1]+
      NumLocalInterpCoefs[2])*sizeof(double));
    double *Buffer = InterpCoefs;
    DatasetOffset = XInterpCoefsOffset;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + NumInterpCoefsBeforeChunk[iDim]*sizeof(double);
      StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(XFile, ReadOffset, MPI_DOUBLE, MPI_DOUBLE, "native", MPIInfo);
      EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(XFile, Buffer, (int)NumLocalInterpCoefs[iDim], MPI_DOUBLE, Endian,
        &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_DOUBLE, &ReadSize);
      if (ReadSize < NumLocalInterpCoefs[iDim]) {
        LogError(Logger, true, "Unable to read grid %i donor interpolation coefficients from "
          "XINTOUT file '%s'.", GridID, XPath);
        Error = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_x);
      }
      Buffer += NumLocalInterpCoefs[iDim];
      DatasetOffset += NumInterpCoefs[iDim]*sizeof(double);
    }

    size_t iNextCoef[MAX_DIMS];
    iNextCoef[0] = 0;
    iNextCoef[1] = iNextCoef[0] + NumLocalInterpCoefs[0];
    iNextCoef[2] = iNextCoef[1] + NumLocalInterpCoefs[1];
    for (iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        int Size = Data->extents[1][iDim][iDonor] - Data->extents[0][iDim][iDonor];
        for (iPoint = 0; iPoint < Size; ++iPoint) {
          Data->interp_coefs[iDim][iPoint][iDonor] = InterpCoefs[iNextCoef[iDim]];
          ++iNextCoef[iDim];
        }
      }
    }

    free(InterpCoefs);

    close_x:
      StartProfileSync(Profiler, MPIIOCloseTime, ChunkComm);
      MPI_File_close(&XFile);
      EndProfile(Profiler, MPIIOCloseTime);

    close_ho:
      StartProfileSync(Profiler, MPIIOCloseTime, ChunkComm);
      MPI_File_close(&HOFile);
      EndProfile(Profiler, MPIIOCloseTime);

    free_comm:
      MPI_Comm_free(&ChunkComm);

  }

  OVK_EH_CHECK_ALL(ErrorHandler, Error, Comm);

  OVK_EH_FINALIZE(ErrorHandler);

  return OVK_NO_ERROR;

}

static ovk_error ReadReceivers(t_xintout_grid *XINTOUTGrid, const char *HOPath, size_t NumReceivers,
  MPI_Offset HOPointsOffset, MPI_Offset HOConnectionIDsOffset, ovk_ext_endian Endian,
  ovk_ext_xintout_format Format, int ReadGranularityAdjust, MPI_Info MPIInfo,
  t_profiler *Profiler) {

  int iDim;
  size_t iReceiver;

  int GridID = XINTOUTGrid->id;

  MPI_Comm Comm = XINTOUTGrid->comm;
  int CommSize = XINTOUTGrid->comm_size;
  int CommRank = XINTOUTGrid->comm_rank;

  t_logger *Logger = XINTOUTGrid->logger;
  t_error_handler *ErrorHandler = XINTOUTGrid->error_handler;

  OVK_EH_INIT(ErrorHandler);
  ovk_error Error = OVK_NO_ERROR;

  t_xintout_receivers *XINTOUTReceivers = &XINTOUTGrid->receivers;

  XINTOUTReceivers->count = NumReceivers;

  // Read chunk on every Nth rank, where N is a power of 2 (or CommSize) chosen such that the number
  // of receivers per chunk is roughly the average number of grid points per rank (subject to user
  // adjustment)
  size_t NumPoints =
    (size_t)XINTOUTGrid->global_size[0]*
    (size_t)XINTOUTGrid->global_size[1]*
    (size_t)XINTOUTGrid->global_size[2];
  size_t AvgPointsPerRank = BinDivide(NumPoints, CommSize);
  int ChunkRankInterval, NumChunks;
  size_t ChunkSize;
  Chunkify(NumReceivers, CommSize, AvgPointsPerRank, ReadGranularityAdjust, &ChunkRankInterval,
    &NumChunks, &ChunkSize);
  bool HasChunk = CommRank % ChunkRankInterval == 0;

  if (LoggingStatus(Logger)) {
    char NumReceiversString[NUMBER_STRING_LENGTH+9];
    char NumRanksString[NUMBER_STRING_LENGTH+9];
    PluralizeLabel(NumReceivers, "receivers", "receiver", NumReceiversString);
    PluralizeLabel(NumChunks, "I/O ranks", "I/O rank", NumRanksString);
    LogStatus(Logger, CommRank == 0, 1, "Grid %s has %s; using %s.",
      XINTOUTGrid->name, NumReceiversString, NumRanksString);
  }

  XINTOUTReceivers->chunk_size = ChunkSize;
  XINTOUTReceivers->has_chunk = HasChunk;

  MPI_Comm ChunkComm;
  MPI_Comm_split(Comm, HasChunk ? 0 : MPI_UNDEFINED, CommRank, &ChunkComm);

  t_xintout_receiver_chunk *Chunk = NULL;

  if (HasChunk) {

    Chunk = malloc(sizeof(t_xintout_receiver_chunk));
    XINTOUTReceivers->chunk = Chunk;

    int ChunkCommRank;
    MPI_Comm_rank(ChunkComm, &ChunkCommRank);

    size_t LocalBegin = ChunkSize*ChunkCommRank;
    size_t LocalEnd = min(ChunkSize*(ChunkCommRank+1), NumReceivers);
    size_t NumLocalReceivers = LocalEnd - LocalBegin;

    if (NumLocalReceivers*sizeof(long long) > INT_MAX) {
      LogError(Logger, true, "Receiver chunk size too big; increase number of processes or read "
        "granularity.");
      Error = OVK_ERROR_FILE_READ;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, free_comm);
    }

    Chunk->begin = LocalBegin;
    Chunk->end = LocalEnd;

    CreateReceiverData(&Chunk->data, NumLocalReceivers);
    t_receiver_data *Data = Chunk->data;

    Chunk->connection_ids = malloc(NumLocalReceivers*sizeof(size_t));

    MPI_File HOFile;
    MPI_Status Status;
    int MPIError;
    MPI_Offset DatasetOffset, ReadOffset;
    int ReadSize;

    int MPIIOOpenTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Open");
    int MPIIOCloseTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Close");
    int MPIIOOtherTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Other");

    StartProfileSync(Profiler, MPIIOOpenTime, ChunkComm);
    MPIError = MPI_File_open(ChunkComm, (char *)HOPath, MPI_MODE_RDONLY, MPIInfo, &HOFile);
    EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      LogError(Logger, true, "Unable to open file '%s'.", HOPath);
      Error = OVK_ERROR_FILE_OPEN;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, free_comm);
    }

//     MPI_Datatype HOReceiverPointType;
//     MPI_Type_vector(MAX_DIMS, NumLocalReceivers, NumReceivers, MPI_INT, &HOReceiverPointType);
//     MPI_Type_commit(&HOReceiverPointType);
//     MPI_File_set_view(HOFile, HOPointsOffset+LocalBegin*sizeof(int), MPI_INT, HOReceiverPointType,
//       "native", MPIInfo);
//     DataSize = MAX_DIMS*NumLocalReceivers;
//     File_read_all_endian(HOFile, Data->points[0], DataSize, MPI_INT, Endian, &Status);
//     MPI_Type_free(&HOReceiverPointType);
//     MPI_Get_count(&Status, MPI_INT, &ReadSize);
//     if (ReadSize < DataSize) {
//       LogError(Logger, true, "Unable to read grid %i receiver points from XINTOUT file '%s'.",
//         GridID, HOPath);
//       Error  = OVK_ERROR_FILE_READ;
//       OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
//     }

    // The above version causes "Conditional jump or move depends on uninitialized value(s)"
    // errors in valgrind with MPICH 3.2... not sure why
    DatasetOffset = HOPointsOffset;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + LocalBegin*sizeof(int);
      StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(HOFile, ReadOffset, MPI_INT, MPI_INT, "native", MPIInfo);
      EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(HOFile, Data->points[iDim], (int)NumLocalReceivers, MPI_INT, Endian,
        &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < NumLocalReceivers) {
        LogError(Logger, true, "Unable to read grid %i receiver points from XINTOUT file '%s'.",
          GridID, HOPath);
        Error  = OVK_ERROR_FILE_READ;
        OVK_EH_HANDLE_GOTO(ErrorHandler, Error, close_ho);
      }
      DatasetOffset += NumReceivers*sizeof(int);
    }

    // Convert to zero-based indexing
    for (iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        --Data->points[iDim][iReceiver];
      }
    }

    int ConnectionIDSize;
    MPI_Datatype ConnectionIDType;
    if (Format == OVK_EXT_XINTOUT_STANDARD) {
      ConnectionIDSize = sizeof(int);
      ConnectionIDType = MPI_INT;
    } else {
      ConnectionIDSize = sizeof(long long);
      ConnectionIDType = MPI_LONG_LONG;
    }

    void *ConnectionIDs = malloc(NumLocalReceivers*ConnectionIDSize);

    DatasetOffset = HOConnectionIDsOffset;
    ReadOffset = DatasetOffset+LocalBegin*ConnectionIDSize;
    StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
    MPI_File_set_view(HOFile, ReadOffset, ConnectionIDType, ConnectionIDType, "native", MPIInfo);
    EndProfile(Profiler, MPIIOOtherTime);
    File_read_all_endian(HOFile, ConnectionIDs, (int)NumLocalReceivers, ConnectionIDType, Endian,
      &Status, Profiler, ChunkComm);
    MPI_Get_count(&Status, ConnectionIDType, &ReadSize);
    if (ReadSize < NumLocalReceivers) {
      LogError(Logger, true, "Unable to read grid %i receiver connection IDs from XINTOUT file "
        "'%s'.", GridID, HOPath);
      Error  = OVK_ERROR_FILE_READ;
      OVK_EH_HANDLE_GOTO(ErrorHandler, Error, free_connection_ids);
    }

    if (Format == OVK_EXT_XINTOUT_STANDARD) {
      for (iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        // Convert to zero-based indexing
        Chunk->connection_ids[iReceiver] = ((int *)ConnectionIDs)[iReceiver]-1;
      }
    } else {
      for (iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        // Convert to zero-based indexing
        Chunk->connection_ids[iReceiver] = ((long long *)ConnectionIDs)[iReceiver]-1;
      }
    }

    free_connection_ids:
      free(ConnectionIDs);
      OVK_EH_CHECK_GOTO(ErrorHandler, Error, close_ho);

    close_ho:
      StartProfileSync(Profiler, MPIIOCloseTime, ChunkComm);
      MPI_File_close(&HOFile);
      EndProfile(Profiler, MPIIOCloseTime);

    free_comm:
      MPI_Comm_free(&ChunkComm);

  }

  OVK_EH_CHECK_ALL(ErrorHandler, Error, Comm);

  OVK_EH_FINALIZE(ErrorHandler);

  return OVK_NO_ERROR;

}

static void MatchDonorsAndReceivers(t_xintout *XINTOUT, t_profiler *Profiler) {

  int iGrid, iLocalGrid;
  int iDim;
  int iBin;
  size_t iDonor, iReceiver;
  size_t iConnection;
  t_ordered_map_entry *Entry;

  MPI_Comm Comm = XINTOUT->comm;
  int CommSize = XINTOUT->comm_size;
  int CommRank = XINTOUT->comm_rank;

  MPI_Barrier(Comm);

  t_logger *Logger = XINTOUT->logger;

  LogStatus(Logger, CommRank == 0, 0, "Matching donors and receivers...");

  int MapToBinsTime = GetProfilerTimerID(Profiler, "XINTOUT::Match::MapToBins");
  int HandshakeTime = GetProfilerTimerID(Profiler, "XINTOUT::Match::Handshake");
  int SendToBinsTime = GetProfilerTimerID(Profiler, "XINTOUT::Match::SendToBins");
  int FillConnectionDataTime = GetProfilerTimerID(Profiler, "XINTOUT::Match::FillConnectionData");
  int RecvFromBinsTime = GetProfilerTimerID(Profiler, "XINTOUT::Match::RecvFromBins");
  int UnpackTime = GetProfilerTimerID(Profiler, "XINTOUT::Match::Unpack");

  int NumGrids = XINTOUT->num_grids;
  int NumLocalGrids = XINTOUT->num_local_grids;

  size_t *NumPointsOnGrid = malloc(NumGrids*sizeof(size_t));
  size_t *NumReceiversOnGrid = malloc(NumGrids*sizeof(size_t));

  for (iGrid = 0; iGrid < NumGrids; ++iGrid) {
    NumPointsOnGrid[iGrid] = 0;
    NumReceiversOnGrid[iGrid] = 0;
  }

  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    t_xintout_grid *XINTOUTGrid = XINTOUT->grids[iLocalGrid];
    int GridID = XINTOUTGrid->id;
    iGrid = GridID-1;
    NumPointsOnGrid[iGrid] = 1;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      NumPointsOnGrid[iGrid] *= XINTOUTGrid->global_size[iDim];
    }
    NumReceiversOnGrid[iGrid] = XINTOUTGrid->receivers.count;
  }

  MPI_Allreduce(MPI_IN_PLACE, NumPointsOnGrid, NumGrids, KMPI_UNSIGNED_SIZE, MPI_MAX, Comm);
  MPI_Allreduce(MPI_IN_PLACE, NumReceiversOnGrid, NumGrids, KMPI_UNSIGNED_SIZE, MPI_MAX, Comm);

  size_t NumPoints = 0;
  size_t NumConnections = 0;
  for (iGrid = 0; iGrid < NumGrids; ++iGrid) {
    NumPoints += NumPointsOnGrid[iGrid];
    NumConnections += NumReceiversOnGrid[iGrid];
  }

  free(NumPointsOnGrid);
  free(NumReceiversOnGrid);

  t_xintout_connections *XINTOUTConnections = &XINTOUT->connections;

  size_t BinSize = min(BinDivide(NumPoints, CommSize), NumConnections);
  bool HasBin = BinSize*CommRank < NumConnections;

  XINTOUTConnections->bin_size = BinSize;
  XINTOUTConnections->has_bin = HasBin;

  t_xintout_connection_bin *Bin = NULL;
  t_connection_data *BinData = NULL;
  if (HasBin) {
    Bin = malloc(sizeof(t_xintout_connection_bin));
    Bin->begin = BinSize*CommRank;
    Bin->end = min(BinSize*(CommRank+1), NumConnections);
    size_t NumLocalConnections = Bin->end - Bin->begin;
    CreateConnectionData(&BinData, NumLocalConnections);
    Bin->data = BinData;
    XINTOUTConnections->bin = Bin;
  }

  typedef struct {
    size_t count;
    t_connection_data *data;
  } t_send_recv;

  StartProfileSync(Profiler, MapToBinsTime, Comm);

  t_ordered_map *DonorSends, *ReceiverSends;
  OMCreate(&DonorSends);
  OMCreate(&ReceiverSends);

  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    t_xintout_grid *XINTOUTGrid = XINTOUT->grids[iLocalGrid];
    t_xintout_donors *XINTOUTDonors = &XINTOUTGrid->donors;
    t_xintout_receivers *XINTOUTReceivers = &XINTOUTGrid->receivers;
    if (XINTOUTDonors->has_chunk) {
      t_xintout_donor_chunk *DonorChunk = XINTOUTDonors->chunk;
      size_t NumLocalDonors = DonorChunk->end - DonorChunk->begin;
      for (iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        size_t ConnectionID = DonorChunk->starting_connection_id + iDonor;
        int ConnectionBinIndex = (int)(ConnectionID/BinSize);
        Entry = OMFind(DonorSends, ConnectionBinIndex);
        t_send_recv *Send;
        if (Entry != OMEnd(DonorSends)) {
          Send = OMData(Entry);
        } else {
          Send = malloc(sizeof(t_send_recv));
          Send->count = 0;
          OMInsert(DonorSends, ConnectionBinIndex, Send);
        }
        ++Send->count;
      }
    }
    if (XINTOUTReceivers->has_chunk) {
      t_xintout_receiver_chunk *ReceiverChunk = XINTOUTReceivers->chunk;
      size_t NumLocalReceivers = ReceiverChunk->end - ReceiverChunk->begin;
      for (iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        size_t ConnectionID = ReceiverChunk->connection_ids[iReceiver];
        int ConnectionBinIndex = (int)(ConnectionID/BinSize);
        Entry = OMFind(ReceiverSends, ConnectionBinIndex);
        t_send_recv *Send;
        if (Entry != OMEnd(ReceiverSends)) {
          Send = OMData(Entry);
        } else {
          Send = malloc(sizeof(t_send_recv));
          Send->count = 0;
          OMInsert(ReceiverSends, ConnectionBinIndex, Send);
        }
        ++Send->count;
      }
    }
  }

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    t_send_recv *Send = OMData(Entry);
    size_t Count = Send->count;
    CreateConnectionData(&Send->data, Count);
    // Reset count for filling in data
    Send->count = 0;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    t_send_recv *Send = OMData(Entry);
    size_t Count = Send->count;
    CreateConnectionData(&Send->data, Count);
    // Reset count for filling in data
    Send->count = 0;
    Entry = OMNext(Entry);
  }

  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    t_xintout_grid *XINTOUTGrid = XINTOUT->grids[iLocalGrid];
    t_xintout_donors *XINTOUTDonors = &XINTOUTGrid->donors;
    t_xintout_receivers *XINTOUTReceivers = &XINTOUTGrid->receivers;
    if (XINTOUTDonors->has_chunk) {
      t_xintout_donor_chunk *DonorChunk = XINTOUTDonors->chunk;
      size_t NumLocalDonors = DonorChunk->end - DonorChunk->begin;
      for (iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        size_t ConnectionID = DonorChunk->starting_connection_id + iDonor;
        int ConnectionBinIndex = (int)(ConnectionID/BinSize);
        Entry = OMFind(DonorSends, ConnectionBinIndex);
        t_send_recv *Send = OMData(Entry);
        size_t iNext = Send->count;
        Send->data->connection_ids[iNext] = ConnectionID;
        Send->data->donor_grid_ids[iNext] = XINTOUTGrid->id;
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send->data->donor_cells[iDim][iNext] = DonorChunk->data->extents[0][iDim][iDonor];
        }
        ++Send->count;
      }
    }
    if (XINTOUTReceivers->has_chunk) {
      t_xintout_receiver_chunk *ReceiverChunk = XINTOUTReceivers->chunk;
      size_t NumLocalReceivers = ReceiverChunk->end - ReceiverChunk->begin;
      for (iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        size_t ConnectionID = ReceiverChunk->connection_ids[iReceiver];
        int ConnectionBinIndex = (int)(ConnectionID/BinSize);
        Entry = OMFind(ReceiverSends, ConnectionBinIndex);
        t_send_recv *Send = OMData(Entry);
        size_t iNext = Send->count;
        Send->data->connection_ids[iNext] = ConnectionID;
        Send->data->receiver_grid_ids[iNext] = XINTOUTGrid->id;
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send->data->receiver_points[iDim][iNext] = ReceiverChunk->data->points[iDim][iReceiver];
        }
        ++Send->count;
      }
    }
  }

  int NumDonorSends = OMSize(DonorSends);
  int NumReceiverSends = OMSize(ReceiverSends);
  int *DonorSendToRanks = malloc(NumDonorSends*sizeof(int));
  int *ReceiverSendToRanks = malloc(NumReceiverSends*sizeof(int));

  Entry = OMBegin(DonorSends);
  iBin = 0;
  while (Entry != OMEnd(DonorSends)) {
    DonorSendToRanks[iBin] = OMKey(Entry);
    Entry = OMNext(Entry);
    ++iBin;
  }

  Entry = OMBegin(ReceiverSends);
  iBin = 0;
  while (Entry != OMEnd(ReceiverSends)) {
    ReceiverSendToRanks[iBin] = OMKey(Entry);
    Entry = OMNext(Entry);
    ++iBin;
  }

  EndProfile(Profiler, MapToBinsTime);
  StartProfileSync(Profiler, HandshakeTime, Comm);

  t_ordered_map *DonorRecvs, *ReceiverRecvs;
  OMCreate(&DonorRecvs);
  OMCreate(&ReceiverRecvs);

  DynamicHandshake(Comm, NumDonorSends, DonorSendToRanks, DonorRecvs);
  DynamicHandshake(Comm, NumReceiverSends, ReceiverSendToRanks, ReceiverRecvs);

  free(DonorSendToRanks);
  free(ReceiverSendToRanks);

  EndProfile(Profiler, HandshakeTime);
  StartProfileSync(Profiler, SendToBinsTime, Comm);

  int NumDonorRecvs = OMSize(DonorRecvs);
  int NumReceiverRecvs = OMSize(ReceiverRecvs);
  // 3 requests per send/recv (1 for each of grid IDs, connection IDs, and points/cells)
  MPI_Request *Requests = malloc(3*(NumDonorSends+NumDonorRecvs+NumReceiverSends+NumReceiverRecvs)*
    sizeof(MPI_Request));

  int NumRequests = 0;

  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    int Rank = OMKey(Entry);
    t_send_recv *Recv = malloc(sizeof(t_send_recv));
    OMSetData(Entry, Recv);
    MPI_Irecv(&Recv->count, 1, KMPI_UNSIGNED_SIZE, Rank, 0, Comm, Requests+NumRequests);
    ++NumRequests;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    int Rank = OMKey(Entry);
    t_send_recv *Recv = malloc(sizeof(t_send_recv));
    OMSetData(Entry, Recv);
    MPI_Irecv(&Recv->count, 1, KMPI_UNSIGNED_SIZE, Rank, 1, Comm, Requests+NumRequests);
    ++NumRequests;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    int Rank = OMKey(Entry);
    t_send_recv *Send = OMData(Entry);
    MPI_Isend(&Send->count, 1, KMPI_UNSIGNED_SIZE, Rank, 0, Comm, Requests+NumRequests);
    ++NumRequests;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    int Rank = OMKey(Entry);
    t_send_recv *Send = OMData(Entry);
    MPI_Isend(&Send->count, 1, KMPI_UNSIGNED_SIZE, Rank, 1, Comm, Requests+NumRequests);
    ++NumRequests;
    Entry = OMNext(Entry);
  }

  MPI_Waitall(NumRequests, Requests, MPI_STATUSES_IGNORE);

  NumRequests = 0;

  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    int Rank = OMKey(Entry);
    t_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    CreateConnectionData(&Recv->data, Count);
    MPI_Irecv(Recv->data->connection_ids, Count, KMPI_UNSIGNED_SIZE, Rank, 0, Comm,
      Requests+NumRequests);
    MPI_Irecv(Recv->data->donor_grid_ids, Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+1);
    MPI_Irecv(Recv->data->donor_cells[0], MAX_DIMS*Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+2);
    NumRequests += 3;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    int Rank = OMKey(Entry);
    t_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    CreateConnectionData(&Recv->data, Count);
    MPI_Irecv(Recv->data->connection_ids, Count, KMPI_UNSIGNED_SIZE, Rank, 1, Comm,
      Requests+NumRequests);
    MPI_Irecv(Recv->data->receiver_grid_ids, Count, MPI_INT, Rank, 1, Comm,
      Requests+NumRequests+1);
    MPI_Irecv(Recv->data->receiver_points[0], MAX_DIMS*Count, MPI_INT, Rank, 1, Comm,
      Requests+NumRequests+2);
    NumRequests += 3;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    int Rank = OMKey(Entry);
    t_send_recv *Send = OMData(Entry);
    int Count = (int)Send->count;
    MPI_Isend(Send->data->connection_ids, Count, KMPI_UNSIGNED_SIZE, Rank, 0, Comm,
      Requests+NumRequests);
    MPI_Isend(Send->data->donor_grid_ids, Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+1);
    MPI_Isend(Send->data->donor_cells[0], MAX_DIMS*Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+2);
    NumRequests += 3;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    int Rank = OMKey(Entry);
    t_send_recv *Send = OMData(Entry);
    int Count = (int)Send->count;
    MPI_Isend(Send->data->connection_ids, Count, KMPI_UNSIGNED_SIZE, Rank, 1, Comm,
      Requests+NumRequests);
    MPI_Isend(Send->data->receiver_grid_ids, Count, MPI_INT, Rank, 1, Comm,
      Requests+NumRequests+1);
    MPI_Isend(Send->data->receiver_points[0], MAX_DIMS*Count, MPI_INT, Rank, 1, Comm,
      Requests+NumRequests+2);
    NumRequests += 3;
    Entry = OMNext(Entry);
  }

  MPI_Waitall(NumRequests, Requests, MPI_STATUSES_IGNORE);

  EndProfile(Profiler, SendToBinsTime);
  StartProfileSync(Profiler, RecvFromBinsTime, Comm);

  NumRequests = 0;

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    int Rank = OMKey(Entry);
    t_send_recv *Send = OMData(Entry);
    int Count = (int)Send->count;
    MPI_Irecv(Send->data->receiver_grid_ids, Count, MPI_INT, Rank, 0, Comm, Requests+NumRequests);
    MPI_Irecv(Send->data->receiver_points[0], MAX_DIMS*Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+1);
    NumRequests += 2;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    int Rank = OMKey(Entry);
    t_send_recv *Send = OMData(Entry);
    int Count = (int)Send->count;
    MPI_Irecv(Send->data->donor_grid_ids, Count, MPI_INT, Rank, 1, Comm, Requests+NumRequests);
    MPI_Irecv(Send->data->donor_cells[0], MAX_DIMS*Count, MPI_INT, Rank, 1, Comm,
      Requests+NumRequests+1);
    NumRequests += 2;
    Entry = OMNext(Entry);
  }

  EndProfile(Profiler, RecvFromBinsTime);
  StartProfileSync(Profiler, FillConnectionDataTime, Comm);

  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    t_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    for (iDonor = 0; iDonor < Count; ++iDonor) {
      iConnection = Recv->data->connection_ids[iDonor] - Bin->begin;
      BinData->donor_grid_ids[iConnection] = Recv->data->donor_grid_ids[iDonor];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        BinData->donor_cells[iDim][iConnection] = Recv->data->donor_cells[iDim][iDonor];
      }
    }
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    t_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    for (iReceiver = 0; iReceiver < Count; ++iReceiver) {
      iConnection = Recv->data->connection_ids[iReceiver] - Bin->begin;
      BinData->receiver_grid_ids[iConnection] = Recv->data->receiver_grid_ids[iReceiver];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        BinData->receiver_points[iDim][iConnection] = Recv->data->receiver_points[iDim][iReceiver];
      }
    }
    Entry = OMNext(Entry);
  }

  EndProfile(Profiler, FillConnectionDataTime);
  StartProfileSync(Profiler, RecvFromBinsTime, Comm);

  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    int Rank = OMKey(Entry);
    t_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    for (iDonor = 0; iDonor < Count; ++iDonor) {
      iConnection = Recv->data->connection_ids[iDonor] - Bin->begin;
      Recv->data->receiver_grid_ids[iDonor] = BinData->receiver_grid_ids[iConnection];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Recv->data->receiver_points[iDim][iDonor] = BinData->receiver_points[iDim][iConnection];
      }
    }
    MPI_Isend(Recv->data->receiver_grid_ids, Count, MPI_INT, Rank, 0, Comm, Requests+NumRequests);
    MPI_Isend(Recv->data->receiver_points[0], MAX_DIMS*Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+1);
    NumRequests += 2;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    int Rank = OMKey(Entry);
    t_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    for (iReceiver = 0; iReceiver < Count; ++iReceiver) {
      iConnection = Recv->data->connection_ids[iReceiver] - Bin->begin;
      Recv->data->donor_grid_ids[iReceiver] = BinData->donor_grid_ids[iConnection];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Recv->data->donor_cells[iDim][iReceiver] = BinData->donor_cells[iDim][iConnection];
      }
    }
    MPI_Isend(Recv->data->donor_grid_ids, Count, MPI_INT, Rank, 1, Comm, Requests+NumRequests);
    MPI_Isend(Recv->data->donor_cells[0], MAX_DIMS*Count, MPI_INT, Rank, 1, Comm,
      Requests+NumRequests+1);
    NumRequests += 2;
    Entry = OMNext(Entry);
  }

  MPI_Waitall(NumRequests, Requests, MPI_STATUSES_IGNORE);

  free(Requests);

  EndProfile(Profiler, RecvFromBinsTime);
  StartProfileSync(Profiler, UnpackTime, Comm);

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    t_send_recv *Send = OMData(Entry);
    // Reset count to 0 for unpacking
    Send->count = 0;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    t_send_recv *Send = OMData(Entry);
    // Reset count to 0 for unpacking
    Send->count = 0;
    Entry = OMNext(Entry);
  }

  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    t_xintout_grid *XINTOUTGrid = XINTOUT->grids[iLocalGrid];
    t_xintout_donors *XINTOUTDonors = &XINTOUTGrid->donors;
    t_xintout_receivers *XINTOUTReceivers = &XINTOUTGrid->receivers;
    if (XINTOUTDonors->has_chunk) {
      t_xintout_donor_chunk *DonorChunk = XINTOUTDonors->chunk;
      size_t NumLocalDonors = DonorChunk->end - DonorChunk->begin;
      for (iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        size_t ConnectionID = DonorChunk->starting_connection_id + iDonor;
        int ConnectionBinIndex = (int)(ConnectionID/BinSize);
        Entry = OMFind(DonorSends, ConnectionBinIndex);
        t_send_recv *Send = OMData(Entry);
        size_t iNext = Send->count;
        DonorChunk->data->destination_grid_ids[iDonor] = Send->data->receiver_grid_ids[iNext];
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          DonorChunk->data->destination_points[iDim][iDonor] =
            Send->data->receiver_points[iDim][iNext];
        }
        ++Send->count;
      }
    }
    if (XINTOUTReceivers->has_chunk) {
      t_xintout_receiver_chunk *ReceiverChunk = XINTOUTReceivers->chunk;
      size_t NumLocalReceivers = ReceiverChunk->end - ReceiverChunk->begin;
      for (iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        size_t ConnectionID = ReceiverChunk->connection_ids[iReceiver];
        int ConnectionBinIndex = (int)(ConnectionID/BinSize);
        Entry = OMFind(ReceiverSends, ConnectionBinIndex);
        t_send_recv *Send = OMData(Entry);
        size_t iNext = Send->count;
        ReceiverChunk->data->source_grid_ids[iReceiver] = Send->data->donor_grid_ids[iNext];
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          ReceiverChunk->data->source_cells[iDim][iReceiver] = Send->data->donor_cells[iDim][iNext];
        }
        ++Send->count;
      }
    }
  }

  EndProfile(Profiler, UnpackTime);

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    t_send_recv *Send = OMData(Entry);
    DestroyConnectionData(&Send->data);
    free(Send);
    Entry = OMNext(Entry);
  }
  OMDestroy(&DonorSends);

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    t_send_recv *Send = OMData(Entry);
    DestroyConnectionData(&Send->data);
    free(Send);
    Entry = OMNext(Entry);
  }
  OMDestroy(&ReceiverSends);

  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    t_send_recv *Recv = OMData(Entry);
    DestroyConnectionData(&Recv->data);
    free(Recv);
    Entry = OMNext(Entry);
  }
  OMDestroy(&DonorRecvs);

  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    t_send_recv *Recv = OMData(Entry);
    DestroyConnectionData(&Recv->data);
    free(Recv);
    Entry = OMNext(Entry);
  }
  OMDestroy(&ReceiverRecvs);

  MPI_Barrier(Comm);

  LogStatus(Logger, CommRank == 0, 0, "Finished matching donors and receivers.");

}

static void DistributeConnectivityData(const t_xintout *XINTOUT, const ovk_grid **LocalGrids,
  t_donor_data **LocalDonorData, t_receiver_data **LocalReceiverData, t_profiler *Profiler) {

  int iLocalGrid;

  MPI_Comm Comm = XINTOUT->comm;
  int CommRank = XINTOUT->comm_rank;

  MPI_Barrier(Comm);

  t_logger *Logger = XINTOUT->logger;

  LogStatus(Logger, CommRank == 0, 0, "Distributing connectivity data to ranks...");

  for (iLocalGrid = 0; iLocalGrid < XINTOUT->num_local_grids; ++iLocalGrid) {
    const t_xintout_grid *XINTOUTGrid = XINTOUT->grids[iLocalGrid];
    const ovk_grid *Grid = LocalGrids[iLocalGrid];
    DistributeGridConnectivityData(XINTOUTGrid, Grid, LocalDonorData+iLocalGrid,
      LocalReceiverData+iLocalGrid, Profiler);
  }

  MPI_Barrier(Comm);

  LogStatus(Logger, CommRank == 0, 0, "Finished distributing connectivity data to ranks.");

}

static void DistributeGridConnectivityData(const t_xintout_grid *XINTOUTGrid, const ovk_grid *Grid,
  t_donor_data **DonorData_, t_receiver_data **ReceiverData_, t_profiler *Profiler) {

  int i, j, k;
  int iDim;
  int iPoint;
  int iPointInCell;
  int iRank;
  int iSend;
  size_t iDonor, iReceiver;
  size_t iDonorPoint;
  size_t iNext;
  t_ordered_map_entry *Entry;

  int NumDims = XINTOUTGrid->num_dims;

  MPI_Comm Comm = XINTOUTGrid->comm;

  int MapToBinsTime = GetProfilerTimerID(Profiler, "XINTOUT::Distribute::MapToBins");
  int RetrieveBinsTime = GetProfilerTimerID(Profiler, "XINTOUT::Distribute::RetrieveBins");
  int FindRanksTime = GetProfilerTimerID(Profiler, "XINTOUT::Distribute::FindRanks");
  int HandshakeTime = GetProfilerTimerID(Profiler, "XINTOUT::Distribute::Handshake");
  int SendDataTime = GetProfilerTimerID(Profiler, "XINTOUT::Distribute::SendData");

  const t_xintout_donors *XINTOUTDonors = &XINTOUTGrid->donors;
  const t_xintout_receivers *XINTOUTReceivers = &XINTOUTGrid->receivers;

  const ovk_grid_properties *Properties;
  ovkGetGridProperties(Grid, &Properties);

  ovk_range GlobalRange;
  ovkGetGridPropertyGlobalRange(Properties, &GlobalRange);

  ovk_cart Cart;
  ovkGetGridCart(Grid, &Cart);

  const t_partition_hash *Hash;
  GetGridPartitionHash(Grid, &Hash);

  StartProfileSync(Profiler, MapToBinsTime, Comm);

  t_ordered_map *Bins;
  OMCreate(&Bins);

  t_xintout_donor_chunk *DonorChunk = NULL;
  size_t NumChunkDonors = 0;
  size_t NumChunkDonorPoints = 0;
  int *ChunkDonorPoints[MAX_DIMS] = {NULL, NULL, NULL};
  int *ChunkDonorPointBinIndices = NULL;
  if (XINTOUTDonors->has_chunk) {
    DonorChunk = XINTOUTDonors->chunk;
    NumChunkDonors = DonorChunk->end - DonorChunk->begin;
    for (iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      int NumPointsInCell = 1;
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        NumPointsInCell *= DonorChunk->data->extents[1][iDim][iDonor] -
          DonorChunk->data->extents[0][iDim][iDonor];
      }
      NumChunkDonorPoints += NumPointsInCell;
    }
    ChunkDonorPoints[0] = malloc(MAX_DIMS*NumChunkDonorPoints*sizeof(int));
    ChunkDonorPoints[1] = ChunkDonorPoints[0] + NumChunkDonorPoints;
    ChunkDonorPoints[2] = ChunkDonorPoints[1] + NumChunkDonorPoints;
    iDonorPoint = 0;
    for (iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      ovk_range DonorRange;
      ovkDefaultRange(&DonorRange, NumDims);
      for (iDim = 0; iDim < NumDims; ++iDim) {
        DonorRange.b[iDim] = DonorChunk->data->extents[0][iDim][iDonor];
        DonorRange.e[iDim] = DonorChunk->data->extents[1][iDim][iDonor];
      }
      bool AwayFromEdge = ovkRangeIncludes(&GlobalRange, &DonorRange);
      if (AwayFromEdge) {
        for (k = DonorRange.b[2]; k < DonorRange.e[2]; ++k) {
          for (j = DonorRange.b[1]; j < DonorRange.e[1]; ++j) {
            for (i = DonorRange.b[0]; i < DonorRange.e[0]; ++i) {
              int Point[MAX_DIMS] = {i, j, k};
              for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
                ChunkDonorPoints[iDim][iDonorPoint] = Point[iDim];
              }
              ++iDonorPoint;
            }
          }
        }
      } else {
        for (k = DonorRange.b[2]; k < DonorRange.e[2]; ++k) {
          for (j = DonorRange.b[1]; j < DonorRange.e[1]; ++j) {
            for (i = DonorRange.b[0]; i < DonorRange.e[0]; ++i) {
              int Point[MAX_DIMS] = {i, j, k};
              ovkCartPeriodicAdjust(&Cart, Point, Point);
              for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
                ChunkDonorPoints[iDim][iDonorPoint] = Point[iDim];
              }
              ++iDonorPoint;
            }
          }
        }
      }
    }
    ChunkDonorPointBinIndices = malloc(NumChunkDonorPoints*sizeof(int));
    MapToPartitionBins(Hash, NumChunkDonorPoints, (const int **)ChunkDonorPoints,
      ChunkDonorPointBinIndices);
    for (iDonorPoint = 0; iDonorPoint < NumChunkDonorPoints; ++iDonorPoint) {
      int BinIndex = ChunkDonorPointBinIndices[iDonorPoint];
      if (!OMExists(Bins, BinIndex)) {
        OMInsert(Bins, BinIndex, NULL);
      }
    }
  }

  t_xintout_receiver_chunk *ReceiverChunk = NULL;
  size_t NumChunkReceivers = 0;
  int *ChunkReceiverBinIndices = NULL;
  if (XINTOUTReceivers->has_chunk) {
    ReceiverChunk = XINTOUTReceivers->chunk;
    NumChunkReceivers = ReceiverChunk->end - ReceiverChunk->begin;
    ChunkReceiverBinIndices = malloc(NumChunkReceivers*sizeof(int));
    MapToPartitionBins(Hash, NumChunkReceivers, (const int **)ReceiverChunk->data->points,
      ChunkReceiverBinIndices);
    for (iReceiver = 0; iReceiver < NumChunkReceivers; ++iReceiver) {
      int BinIndex = ChunkReceiverBinIndices[iReceiver];
      if (!OMExists(Bins, BinIndex)) {
        OMInsert(Bins, BinIndex, NULL);
      }
    }
  }

  EndProfile(Profiler, MapToBinsTime);
  StartProfileSync(Profiler, RetrieveBinsTime, Comm);

  RetrievePartitionBins(Hash, Bins);

  EndProfile(Profiler, RetrieveBinsTime);
  StartProfileSync(Profiler, FindRanksTime, Comm);

  int *NumChunkDonorRanks;
  int **ChunkDonorRanks;
  if (XINTOUTDonors->has_chunk) {
    NumChunkDonorRanks = malloc(NumChunkDonors*sizeof(int));
    ChunkDonorRanks = malloc(NumChunkDonors*sizeof(int *));
    ChunkDonorRanks[0] = malloc(NumChunkDonorPoints*sizeof(int));
    FindPartitions(Hash, Bins, NumChunkDonorPoints, (const int **)ChunkDonorPoints,
      ChunkDonorPointBinIndices, ChunkDonorRanks[0]);
    int MaxSize = DonorChunk->data->max_size;
    int MaxPointsInCell = 1;
    for (iDim = 0; iDim < NumDims; ++iDim) {
      MaxPointsInCell *= MaxSize;
    }
    int *UniqueRanks = malloc(MaxPointsInCell*sizeof(int));
    size_t iDonorPoint = 0;
    for (iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      NumChunkDonorRanks[iDonor] = 0;
      ChunkDonorRanks[iDonor] = ChunkDonorRanks[0] + iDonorPoint;
      int NumPointsInCell = 1;
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        NumPointsInCell *= DonorChunk->data->extents[1][iDim][iDonor] -
          DonorChunk->data->extents[0][iDim][iDonor];
      }
      for (iPointInCell = 0; iPointInCell < NumPointsInCell; ++iPointInCell) {
        int Rank = ChunkDonorRanks[iDonor][iPointInCell];
        iRank = 0;
        while (iRank < NumChunkDonorRanks[iDonor] && UniqueRanks[iRank] != Rank) {
          ++iRank;
        }
        if (iRank == NumChunkDonorRanks[iDonor]) {
          UniqueRanks[iRank] = Rank;
          ++NumChunkDonorRanks[iDonor];
        }
      }
      for (iRank = 0; iRank < NumChunkDonorRanks[iDonor]; ++iRank) {
        ChunkDonorRanks[iDonor][iRank] = UniqueRanks[iRank];
      }
      iDonorPoint += NumPointsInCell;
    }
    free(UniqueRanks);
    free(ChunkDonorPoints[0]);
    free(ChunkDonorPointBinIndices);
  }

  int *ChunkReceiverRanks;
  if (XINTOUTReceivers->has_chunk) {
    ChunkReceiverRanks = malloc(NumChunkReceivers*sizeof(int));
    FindPartitions(Hash, Bins, NumChunkReceivers, (const int **)ReceiverChunk->data->points,
      ChunkReceiverBinIndices, ChunkReceiverRanks);
    free(ChunkReceiverBinIndices);
  }

  ClearPartitionBins(Bins);
  OMDestroy(&Bins);

  typedef struct {
    size_t count;
    int max_size;
    t_donor_data *data;
  } t_donor_send_recv;

  typedef struct {
    size_t count;
    t_receiver_data *data;
  } t_receiver_send_recv;

  t_ordered_map *DonorSends, *ReceiverSends;
  OMCreate(&DonorSends);
  OMCreate(&ReceiverSends);

  if (XINTOUTDonors->has_chunk) {
    for (iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      for (iRank = 0; iRank < NumChunkDonorRanks[iDonor]; ++iRank) {
        int Rank = ChunkDonorRanks[iDonor][iRank];
        Entry = OMFind(DonorSends, Rank);
        t_donor_send_recv *Send;
        if (Entry != OMEnd(DonorSends)) {
          Send = OMData(Entry);
        } else {
          Send = malloc(sizeof(t_donor_send_recv));
          Send->count = 0;
          Send->max_size = 0;
          Send->data = NULL;
          OMInsert(DonorSends, Rank, Send);
        }
        ++Send->count;
        Send->max_size = max(Send->max_size, DonorChunk->data->max_size);
      }
    }
  }

  if (XINTOUTReceivers->has_chunk) {
    for (iReceiver = 0; iReceiver < NumChunkReceivers; ++iReceiver) {
      int Rank = ChunkReceiverRanks[iReceiver];
      Entry = OMFind(ReceiverSends, Rank);
      t_receiver_send_recv *Send;
      if (Entry != OMEnd(ReceiverSends)) {
        Send = OMData(Entry);
      } else {
        Send = malloc(sizeof(t_receiver_send_recv));
        Send->count = 0;
        Send->data = NULL;
        OMInsert(ReceiverSends, Rank, Send);
      }
      ++Send->count;
    }
  }

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    t_donor_send_recv *Send = OMData(Entry);
    CreateDonorData(&Send->data, Send->count, Send->max_size);
    // Reset count for filling in data
    Send->count = 0;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    t_receiver_send_recv *Send = OMData(Entry);
    CreateReceiverData(&Send->data, Send->count);
    // Reset count for filling in data
    Send->count = 0;
    Entry = OMNext(Entry);
  }

  if (XINTOUTDonors->has_chunk) {
    for (iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      for (iRank = 0; iRank < NumChunkDonorRanks[iDonor]; ++iRank) {
        int Rank = ChunkDonorRanks[iDonor][iRank];
        t_donor_send_recv *Send = OMData(OMFind(DonorSends, Rank));
        iNext = Send->count;
        int MaxSize = Send->max_size;
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send->data->extents[0][iDim][iNext] = DonorChunk->data->extents[0][iDim][iDonor];
          Send->data->extents[1][iDim][iNext] = DonorChunk->data->extents[1][iDim][iDonor];
        }
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send->data->coords[iDim][iNext] = DonorChunk->data->coords[iDim][iDonor];
        }
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          int Size = DonorChunk->data->extents[1][iDim][iDonor] -
            DonorChunk->data->extents[0][iDim][iDonor];
          for (iPoint = 0; iPoint < Size; ++iPoint) {
            Send->data->interp_coefs[iDim][iPoint][iNext] =
              DonorChunk->data->interp_coefs[iDim][iPoint][iDonor];
          }
          // Initialize the rest with zeros since we're technically touching all of the data
          // when sending it
          for (iPoint = Size; iPoint < MaxSize; ++iPoint) {
            Send->data->interp_coefs[iDim][iPoint][iNext] = 0.;
          }
        }
        Send->data->destination_grid_ids[iNext] = DonorChunk->data->destination_grid_ids[iDonor];
        for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send->data->destination_points[iDim][iNext] =
            DonorChunk->data->destination_points[iDim][iDonor];
        }
        ++Send->count;
      }
    }
    free(NumChunkDonorRanks);
    free(ChunkDonorRanks[0]);
    free(ChunkDonorRanks);
  }

  if (XINTOUTReceivers->has_chunk) {
    for (iReceiver = 0; iReceiver < NumChunkReceivers; ++iReceiver) {
      int Rank = ChunkReceiverRanks[iReceiver];
      t_receiver_send_recv *Send = OMData(OMFind(ReceiverSends, Rank));
      iNext = Send->count;
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Send->data->points[iDim][iNext] = ReceiverChunk->data->points[iDim][iReceiver];
      }
      Send->data->source_grid_ids[iNext] = ReceiverChunk->data->source_grid_ids[iReceiver];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Send->data->source_cells[iDim][iNext] = ReceiverChunk->data->source_cells[iDim][iReceiver];
      }
      ++Send->count;
    }
    free(ChunkReceiverRanks);
  }

  int NumDonorSends = OMSize(DonorSends);
  int NumReceiverSends = OMSize(ReceiverSends);
  int *DonorSendToRanks = malloc(NumDonorSends*sizeof(int));
  int *ReceiverSendToRanks = malloc(NumReceiverSends*sizeof(int));

  Entry = OMBegin(DonorSends);
  iSend = 0;
  while (Entry != OMEnd(DonorSends)) {
    DonorSendToRanks[iSend] = OMKey(Entry);
    Entry = OMNext(Entry);
    ++iSend;
  }

  Entry = OMBegin(ReceiverSends);
  iSend = 0;
  while (Entry != OMEnd(ReceiverSends)) {
    ReceiverSendToRanks[iSend] = OMKey(Entry);
    Entry = OMNext(Entry);
    ++iSend;
  }

  EndProfile(Profiler, FindRanksTime);
  StartProfileSync(Profiler, HandshakeTime, Comm);

  t_ordered_map *DonorRecvs, *ReceiverRecvs;
  OMCreate(&DonorRecvs);
  OMCreate(&ReceiverRecvs);

  DynamicHandshake(Comm, NumDonorSends, DonorSendToRanks, DonorRecvs);
  DynamicHandshake(Comm, NumReceiverSends, ReceiverSendToRanks, ReceiverRecvs);

  free(DonorSendToRanks);
  free(ReceiverSendToRanks);

  EndProfile(Profiler, HandshakeTime);
  StartProfileSync(Profiler, SendDataTime, Comm);

  int NumDonorRecvs = OMSize(DonorRecvs);
  int NumReceiverRecvs = OMSize(ReceiverRecvs);
  // 5 requests per donor send/recv (1 for each of extents, coords, interp coefs, destination grid
  // IDs, and destination points)
  // 3 requests per receiver send/recv (1 for each of points, source grid IDs, and source cells)
  MPI_Request *Requests = malloc((5*(NumDonorSends+NumDonorRecvs)+3*(NumReceiverSends+
    NumReceiverRecvs))*sizeof(MPI_Request));

  int NumRequests = 0;

  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    int Rank = OMKey(Entry);
    t_donor_send_recv *Recv = malloc(sizeof(t_donor_send_recv));
    OMSetData(Entry, Recv);
    MPI_Irecv(&Recv->count, 1, KMPI_UNSIGNED_SIZE, Rank, 0, Comm, Requests+NumRequests);
    MPI_Irecv(&Recv->max_size, 1, MPI_INT, Rank, 0, Comm, Requests+NumRequests+1);
    NumRequests += 2;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    int Rank = OMKey(Entry);
    t_receiver_send_recv *Recv = malloc(sizeof(t_receiver_send_recv));
    OMSetData(Entry, Recv);
    MPI_Irecv(&Recv->count, 1, KMPI_UNSIGNED_SIZE, Rank, 1, Comm, Requests+NumRequests);
    ++NumRequests;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    int Rank = OMKey(Entry);
    t_donor_send_recv *Send = OMData(Entry);
    MPI_Isend(&Send->count, 1, KMPI_UNSIGNED_SIZE, Rank, 0, Comm, Requests+NumRequests);
    MPI_Isend(&Send->max_size, 1, MPI_INT, Rank, 0, Comm, Requests+NumRequests+1);
    NumRequests += 2;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    int Rank = OMKey(Entry);
    t_receiver_send_recv *Send = OMData(Entry);
    MPI_Isend(&Send->count, 1, KMPI_UNSIGNED_SIZE, Rank, 1, Comm, Requests+NumRequests);
    ++NumRequests;
    Entry = OMNext(Entry);
  }

  MPI_Waitall(NumRequests, Requests, MPI_STATUSES_IGNORE);

  NumRequests = 0;

  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    int Rank = OMKey(Entry);
    t_donor_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    int MaxSize = Recv->max_size;
    CreateDonorData(&Recv->data, Count, MaxSize);
    MPI_Irecv(Recv->data->extents[0][0], 2*MAX_DIMS*Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests);
    MPI_Irecv(Recv->data->coords[0], MAX_DIMS*Count, MPI_DOUBLE, Rank, 0, Comm,
      Requests+NumRequests+1);
    MPI_Irecv(Recv->data->interp_coefs[0][0], MAX_DIMS*MaxSize*Count, MPI_DOUBLE, Rank, 0, Comm,
      Requests+NumRequests+2);
    MPI_Irecv(Recv->data->destination_grid_ids, Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+3);
    MPI_Irecv(Recv->data->destination_points[0], MAX_DIMS*Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+4);
    NumRequests += 5;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    int Rank = OMKey(Entry);
    t_receiver_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    CreateReceiverData(&Recv->data, Count);
    MPI_Irecv(Recv->data->points[0], MAX_DIMS*Count, MPI_INT, Rank, 1, Comm, Requests+NumRequests);
    MPI_Irecv(Recv->data->source_grid_ids, Count, MPI_INT, Rank, 1, Comm, Requests+NumRequests+1);
    MPI_Irecv(Recv->data->source_cells[0], MAX_DIMS*Count, MPI_INT, Rank, 1, Comm,
      Requests+NumRequests+2);
    NumRequests += 3;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    int Rank = OMKey(Entry);
    t_donor_send_recv *Send = OMData(Entry);
    int Count = (int)Send->count;
    int MaxSize = Send->max_size;
    MPI_Isend(Send->data->extents[0][0], 2*MAX_DIMS*Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests);
    MPI_Isend(Send->data->coords[0], MAX_DIMS*Count, MPI_DOUBLE, Rank, 0, Comm,
      Requests+NumRequests+1);
    MPI_Isend(Send->data->interp_coefs[0][0], MAX_DIMS*MaxSize*Count, MPI_DOUBLE, Rank, 0, Comm,
      Requests+NumRequests+2);
    MPI_Isend(Send->data->destination_grid_ids, Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+3);
    MPI_Isend(Send->data->destination_points[0], MAX_DIMS*Count, MPI_INT, Rank, 0, Comm,
      Requests+NumRequests+4);
    NumRequests += 5;
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    int Rank = OMKey(Entry);
    t_receiver_send_recv *Send = OMData(Entry);
    int Count = (int)Send->count;
    MPI_Isend(Send->data->points[0], MAX_DIMS*Count, MPI_INT, Rank, 1, Comm, Requests+NumRequests);
    MPI_Isend(Send->data->source_grid_ids, Count, MPI_INT, Rank, 1, Comm, Requests+NumRequests+1);
    MPI_Isend(Send->data->source_cells[0], MAX_DIMS*Count, MPI_INT, Rank, 1, Comm,
      Requests+NumRequests+2);
    NumRequests += 3;
    Entry = OMNext(Entry);
  }

  MPI_Waitall(NumRequests, Requests, MPI_STATUSES_IGNORE);

  free(Requests);

  Entry = OMBegin(DonorSends);
  while (Entry != OMEnd(DonorSends)) {
    t_donor_send_recv *Send = OMData(Entry);
    DestroyDonorData(&Send->data);
    free(Send);
    Entry = OMNext(Entry);
  }
  OMDestroy(&DonorSends);

  Entry = OMBegin(ReceiverSends);
  while (Entry != OMEnd(ReceiverSends)) {
    t_receiver_send_recv *Send = OMData(Entry);
    DestroyReceiverData(&Send->data);
    free(Send);
    Entry = OMNext(Entry);
  }
  OMDestroy(&ReceiverSends);

  size_t NumLocalDonors = 0;
  int MaxSize = 0;
  size_t NumLocalReceivers = 0;

  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    t_donor_send_recv *Recv = OMData(Entry);
    NumLocalDonors += Recv->count;
    MaxSize = max(MaxSize, Recv->max_size);
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    t_receiver_send_recv *Recv = OMData(Entry);
    NumLocalReceivers += Recv->count;
    Entry = OMNext(Entry);
  }

  CreateDonorData(DonorData_, NumLocalDonors, MaxSize);
  t_donor_data *DonorData = *DonorData_;

  iNext = 0;
  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    t_donor_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    for (iDonor = 0; iDonor < Count; ++iDonor) {
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorData->extents[0][iDim][iNext] = Recv->data->extents[0][iDim][iDonor];
        DonorData->extents[1][iDim][iNext] = Recv->data->extents[1][iDim][iDonor];
      }
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorData->coords[iDim][iNext] = Recv->data->coords[iDim][iDonor];
      }
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        int Size = Recv->data->extents[1][iDim][iDonor] - Recv->data->extents[0][iDim][iDonor];
        for (iPoint = 0; iPoint < Size; ++iPoint) {
          DonorData->interp_coefs[iDim][iPoint][iNext] =
            Recv->data->interp_coefs[iDim][iPoint][iDonor];
        }
      }
      DonorData->destination_grid_ids[iNext] = Recv->data->destination_grid_ids[iDonor];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorData->destination_points[iDim][iNext] = Recv->data->destination_points[iDim][iDonor];
      }
      ++iNext;
    }
    Entry = OMNext(Entry);
  }

  CreateReceiverData(ReceiverData_, NumLocalReceivers);
  t_receiver_data *ReceiverData = *ReceiverData_;

  iNext = 0;
  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    t_receiver_send_recv *Recv = OMData(Entry);
    int Count = (int)Recv->count;
    for (iReceiver = 0; iReceiver < Count; ++iReceiver) {
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        ReceiverData->points[iDim][iNext] = Recv->data->points[iDim][iReceiver];
      }
      ReceiverData->source_grid_ids[iNext] = Recv->data->source_grid_ids[iReceiver];
      for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
        ReceiverData->source_cells[iDim][iNext] = Recv->data->source_cells[iDim][iReceiver];
      }
      ++iNext;
    }
    Entry = OMNext(Entry);
  }

  Entry = OMBegin(DonorRecvs);
  while (Entry != OMEnd(DonorRecvs)) {
    t_donor_send_recv *Recv = OMData(Entry);
    DestroyDonorData(&Recv->data);
    free(Recv);
    Entry = OMNext(Entry);
  }
  OMDestroy(&DonorRecvs);

  Entry = OMBegin(ReceiverRecvs);
  while (Entry != OMEnd(ReceiverRecvs)) {
    t_receiver_send_recv *Recv = OMData(Entry);
    DestroyReceiverData(&Recv->data);
    free(Recv);
    Entry = OMNext(Entry);
  }
  OMDestroy(&ReceiverRecvs);

  EndProfile(Profiler, SendDataTime);

}

static void ImportConnectivityData(int NumGrids, int NumLocalGrids, int *LocalGridIDs,
  const t_donor_data **LocalDonorData, const t_receiver_data **LocalReceiverData,
  ovk_domain *Domain) {

  int iLocalGrid, iGrid, jGrid;
  int iPair;
  size_t iReceiver;

  const ovk_domain_properties *DomainProperties;
  ovkGetDomainProperties(Domain, &DomainProperties);

  char DomainName[OVK_NAME_LENGTH];
  ovkGetDomainPropertyName(DomainProperties, DomainName);

  MPI_Comm Comm;
  int CommRank;
  ovkGetDomainPropertyComm(DomainProperties, &Comm);
  ovkGetDomainPropertyCommRank(DomainProperties, &CommRank);

  MPI_Barrier(Comm);

  t_logger *Logger;
  GetDomainLogger(Domain, &Logger);

  LogStatus(Logger, CommRank == 0, 0, "Importing connectivity data into domain %s...", DomainName);

  size_t **NumConnections = malloc(NumGrids*sizeof(size_t *));
  NumConnections[0] = malloc(NumGrids*NumGrids*sizeof(size_t));
  for (iGrid = 1; iGrid < NumGrids; ++iGrid) {
    NumConnections[iGrid] = NumConnections[iGrid-1] + NumGrids;
  }
  for (iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (jGrid = 0; jGrid < NumGrids; ++jGrid) {
      NumConnections[iGrid][jGrid] = 0;
    }
  }

  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    jGrid = LocalGridIDs[iLocalGrid]-1;
    for (iReceiver = 0; iReceiver < LocalReceiverData[iLocalGrid]->count; ++iReceiver) {
      iGrid = LocalReceiverData[iLocalGrid]->source_grid_ids[iReceiver]-1;
      ++NumConnections[iGrid][jGrid];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, NumConnections[0], NumGrids*NumGrids, KMPI_UNSIGNED_SIZE, MPI_SUM,
    Comm);

  typedef struct {
    ovk_connectivity *connectivity;
    ovk_connectivity_d *donors;
    ovk_connectivity_r *receivers;
  } t_local_connectivity;

  t_ordered_map *LocalConnectivities;
  OMCreate(&LocalConnectivities);

  for (iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections[iGrid][jGrid] > 0) {
        int DonorGridID = iGrid+1;
        int ReceiverGridID = jGrid+1;
        bool DonorGridIsLocal = ovkRankHasGrid(Domain, DonorGridID);
        bool ReceiverGridIsLocal = ovkRankHasGrid(Domain, ReceiverGridID);
        if (DonorGridIsLocal || ReceiverGridIsLocal) {
          iPair = NumGrids*iGrid + jGrid;
          t_local_connectivity *Data = malloc(sizeof(t_local_connectivity));
          ovk_connectivity *Connectivity;
          ovkEditConnectivityLocal(Domain, DonorGridID, ReceiverGridID, &Connectivity);
          Data->connectivity = Connectivity;
          if (DonorGridIsLocal) {
            ovk_connectivity_d *Donors;
            ovkEditConnectivityDonorSideLocal(Connectivity, &Donors);
            Data->donors = Donors;
          } else {
            ovkEditConnectivityDonorSideRemote(Connectivity);
            Data->donors = NULL;
          }
          if (ReceiverGridIsLocal) {
            ovk_connectivity_r *Receivers;
            ovkEditConnectivityReceiverSideLocal(Connectivity, &Receivers);
            Data->receivers = Receivers;
          } else {
            ovkEditConnectivityReceiverSideRemote(Connectivity);
            Data->receivers = NULL;
          }
          OMInsert(LocalConnectivities, iPair, Data);
        } else {
          ovkEditConnectivityRemote(Domain, DonorGridID, ReceiverGridID);
        }
      }
    }
  }

  for (iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {

    int GridID = LocalGridIDs[iLocalGrid];
    const ovk_grid *Grid;
    ovkGetGrid(Domain, GridID, &Grid);
    int iGrid = GridID-1;

    int NumReceiverGrids = 0;
    int NumDonorGrids = 0;
    for (jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections[iGrid][jGrid] > 0) {
        ++NumReceiverGrids;
      }
      if (NumConnections[jGrid][iGrid] > 0) {
        ++NumDonorGrids;
      }
    }

    ovk_connectivity_d **Donors = malloc(NumReceiverGrids*sizeof(ovk_connectivity_d *));
    ovk_connectivity_r **Receivers = malloc(NumDonorGrids*sizeof(ovk_connectivity_r *));

    int iReceiverGrid = 0;
    int iDonorGrid = 0;
    for (jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections[iGrid][jGrid] > 0) {
        iPair = NumGrids*iGrid + jGrid;
        t_local_connectivity *Data = OMData(OMFind(LocalConnectivities, iPair));
        Donors[iReceiverGrid] = Data->donors;
        ++iReceiverGrid;
      }
      if (NumConnections[jGrid][iGrid] > 0) {
        iPair = NumGrids*jGrid + iGrid;
        t_local_connectivity *Data = OMData(OMFind(LocalConnectivities, iPair));
        Receivers[iDonorGrid] = Data->receivers;
        ++iDonorGrid;
      }
    }

    ImportDonors(LocalDonorData[iLocalGrid], NumReceiverGrids, Donors);
    ImportReceivers(LocalReceiverData[iLocalGrid], NumDonorGrids, Receivers);

    free(Donors);
    free(Receivers);

  }

  for (iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections[iGrid][jGrid] > 0) {
        int DonorGridID = iGrid+1;
        int ReceiverGridID = jGrid+1;
        bool DonorGridIsLocal = ovkRankHasGrid(Domain, DonorGridID);
        bool ReceiverGridIsLocal = ovkRankHasGrid(Domain, ReceiverGridID);
        if (DonorGridIsLocal || ReceiverGridIsLocal) {
          iPair = NumGrids*iGrid + jGrid;
          t_local_connectivity *Data = OMData(OMFind(LocalConnectivities, iPair));
          ovk_connectivity *Connectivity = Data->connectivity;
          if (DonorGridIsLocal) {
            ovk_connectivity_d *Donors = Data->donors;
            ovkReleaseConnectivityDonorSideLocal(Connectivity, &Donors);
          } else {
            ovkReleaseConnectivityDonorSideRemote(Connectivity);
          }
          if (ReceiverGridIsLocal) {
            ovk_connectivity_r *Receivers = Data->receivers;
            ovkReleaseConnectivityReceiverSideLocal(Connectivity, &Receivers);
          } else {
            ovkReleaseConnectivityReceiverSideRemote(Connectivity);
          }
          ovkReleaseConnectivityLocal(Domain, DonorGridID, ReceiverGridID, &Connectivity);
          free(Data);
        } else {
          ovkReleaseConnectivityRemote(Domain, DonorGridID, ReceiverGridID);
        }
      }
    }
  }

  OMDestroy(&LocalConnectivities);

  free(NumConnections[0]);
  free(NumConnections);

  MPI_Barrier(Comm);

  LogStatus(Logger, CommRank == 0, 0, "Finished importing connectivity data into domain %s.",
    DomainName);

}

static void ImportDonors(const t_donor_data *DonorData, int NumDestinationGrids,
  ovk_connectivity_d **Donors) {

  int iDestinationGrid;
  int iDim;
  int iPoint;
  size_t iDonor;

  if (NumDestinationGrids == 0) return;

  const ovk_connectivity_d_properties *Properties;
  ovkGetConnectivityDonorSideProperties(Donors[0], &Properties);

  MPI_Comm Comm;
  ovkGetConnectivityDonorSidePropertyComm(Properties, &Comm);

  int *DestinationGridIDs = malloc(NumDestinationGrids*sizeof(int));

  for (iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    ovkGetConnectivityDonorSideProperties(Donors[iDestinationGrid], &Properties);
    ovkGetConnectivityDonorSidePropertyDestinationGridID(Properties, DestinationGridIDs+
      iDestinationGrid);
  }

  t_ordered_map *OverkitDonorData;
  OMCreate(&OverkitDonorData);

  for (iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    t_donor_data *OverkitData = malloc(sizeof(t_donor_data));
    OverkitData->count = 0;
    OverkitData->max_size = 0;
    OMInsert(OverkitDonorData, DestinationGridIDs[iDestinationGrid], OverkitData);
  }

  for (iDonor = 0; iDonor < DonorData->count; ++iDonor) {
    int DestinationGridID = DonorData->destination_grid_ids[iDonor];
    t_donor_data *OverkitData = OMData(OMFind(OverkitDonorData, DestinationGridID));
    ++OverkitData->count;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      int Size = DonorData->extents[1][iDim][iDonor] - DonorData->extents[0][iDim][iDonor];
      OverkitData->max_size = max(OverkitData->max_size, Size);
    }
  }

  t_ordered_map_entry *Entry = OMBegin(OverkitDonorData);
  while (Entry != OMEnd(OverkitDonorData)) {
    t_donor_data *OverkitData = OMData(Entry);
    MPI_Allreduce(MPI_IN_PLACE, &OverkitData->max_size, 1, MPI_INT, MPI_MAX, Comm);
    Entry = OMNext(Entry);
  }

  for (iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    t_donor_data *OverkitData = OMData(OMFind(OverkitDonorData,
      DestinationGridIDs[iDestinationGrid]));
    int MaxSize = OverkitData->max_size;
    ovkResizeDonors(Donors[iDestinationGrid], OverkitData->count, MaxSize);
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkEditDonorExtents(Donors[iDestinationGrid], iDim, &OverkitData->extents[0][iDim],
        &OverkitData->extents[1][iDim]);
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkEditDonorCoords(Donors[iDestinationGrid], iDim, &OverkitData->coords[iDim]);
    }
    OverkitData->interp_coefs[0] = malloc(MAX_DIMS*MaxSize*sizeof(double *));
    OverkitData->interp_coefs[1] = OverkitData->interp_coefs[0] + MaxSize;
    OverkitData->interp_coefs[2] = OverkitData->interp_coefs[1] + MaxSize;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (iPoint = 0; iPoint < MaxSize; ++iPoint) {
        ovkEditDonorInterpCoefs(Donors[iDestinationGrid], iDim, iPoint,
          &OverkitData->interp_coefs[iDim][iPoint]);
      }
    }
    OverkitData->destination_grid_ids = NULL;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkEditDonorDestinations(Donors[iDestinationGrid], iDim,
        &OverkitData->destination_points[iDim]);
    }
    // Reset count to 0 for filling in data
    OverkitData->count = 0;
  }

  for (iDonor = 0; iDonor < DonorData->count; ++iDonor) {
    int DestinationGridID = DonorData->destination_grid_ids[iDonor];
    t_donor_data *OverkitData = OMData(OMFind(OverkitDonorData, DestinationGridID));
    size_t iNext = OverkitData->count;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      OverkitData->extents[0][iDim][iNext] = DonorData->extents[0][iDim][iDonor];
      OverkitData->extents[1][iDim][iNext] = DonorData->extents[1][iDim][iDonor];
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      OverkitData->coords[iDim][iNext] = DonorData->coords[iDim][iDonor];
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      int Size = DonorData->extents[1][iDim][iDonor]-DonorData->extents[0][iDim][iDonor];
      for (iPoint = 0; iPoint < Size; ++iPoint) {
        OverkitData->interp_coefs[iDim][iPoint][iNext] =
          DonorData->interp_coefs[iDim][iPoint][iDonor];
      }
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      OverkitData->destination_points[iDim][iNext] = DonorData->destination_points[iDim][iDonor];
    }
    ++OverkitData->count;
  }

  for (iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    t_donor_data *OverkitData = OMData(OMFind(OverkitDonorData,
      DestinationGridIDs[iDestinationGrid]));
    int MaxSize = OverkitData->max_size;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkReleaseDonorExtents(Donors[iDestinationGrid], iDim, &OverkitData->extents[0][iDim],
        &OverkitData->extents[1][iDim]);
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkReleaseDonorCoords(Donors[iDestinationGrid], iDim, &OverkitData->coords[iDim]);
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (iPoint = 0; iPoint < MaxSize; ++iPoint) {
        ovkReleaseDonorInterpCoefs(Donors[iDestinationGrid], iDim, iPoint,
          &OverkitData->interp_coefs[iDim][iPoint]);
      }
    }
    free(OverkitData->interp_coefs[0]);
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkReleaseDonorDestinations(Donors[iDestinationGrid], iDim,
        &OverkitData->destination_points[iDim]);
    }
    free(OverkitData);
  }

  OMDestroy(&OverkitDonorData);

  free(DestinationGridIDs);

}

static void ImportReceivers(const t_receiver_data *ReceiverData, int NumSourceGrids,
  ovk_connectivity_r **Receivers) {

  int iSourceGrid;
  int iDim;
  size_t iReceiver;

  if (NumSourceGrids == 0) return;

  const ovk_connectivity_r_properties *Properties;
  ovkGetConnectivityReceiverSideProperties(Receivers[0], &Properties);

  MPI_Comm Comm;
  ovkGetConnectivityReceiverSidePropertyComm(Properties, &Comm);

  int *SourceGridIDs = malloc(NumSourceGrids*sizeof(int));

  for (iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    ovkGetConnectivityReceiverSideProperties(Receivers[iSourceGrid], &Properties);
    ovkGetConnectivityReceiverSidePropertySourceGridID(Properties, SourceGridIDs+iSourceGrid);
  }

  t_ordered_map *OverkitReceiverData;
  OMCreate(&OverkitReceiverData);

  for (iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    t_receiver_data *OverkitData = malloc(sizeof(t_receiver_data));
    OverkitData->count = 0;
    OMInsert(OverkitReceiverData, SourceGridIDs[iSourceGrid], OverkitData);
  }

  for (iReceiver = 0; iReceiver < ReceiverData->count; ++iReceiver) {
    int SourceGridID = ReceiverData->source_grid_ids[iReceiver];
    t_receiver_data *OverkitData = OMData(OMFind(OverkitReceiverData, SourceGridID));
    ++OverkitData->count;
  }

  for (iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    t_receiver_data *OverkitData = OMData(OMFind(OverkitReceiverData, SourceGridIDs[iSourceGrid]));
    ovkResizeReceivers(Receivers[iSourceGrid], OverkitData->count);
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkEditReceiverPoints(Receivers[iSourceGrid], iDim, &OverkitData->points[iDim]);
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkEditReceiverSources(Receivers[iSourceGrid], iDim, &OverkitData->source_cells[iDim]);
    }
    // Reset count to 0 for filling in data
    OverkitData->count = 0;
  }

  for (iReceiver = 0; iReceiver < ReceiverData->count; ++iReceiver) {
    int SourceGridID = ReceiverData->source_grid_ids[iReceiver];
    t_receiver_data *OverkitData = OMData(OMFind(OverkitReceiverData, SourceGridID));
    int iNext = OverkitData->count;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      OverkitData->points[iDim][iNext] = ReceiverData->points[iDim][iReceiver];
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      OverkitData->source_cells[iDim][iNext] = ReceiverData->source_cells[iDim][iReceiver];
    }
    ++OverkitData->count;
  }

  for (iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    t_receiver_data *OverkitData = OMData(OMFind(OverkitReceiverData, SourceGridIDs[iSourceGrid]));
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkReleaseReceiverPoints(Receivers[iSourceGrid], iDim, &OverkitData->points[iDim]);
    }
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ovkReleaseReceiverSources(Receivers[iSourceGrid], iDim, &OverkitData->source_cells[iDim]);
    }
    free(OverkitData);
  }

  OMDestroy(&OverkitReceiverData);

  free(SourceGridIDs);

}

static void CreateDonorData(t_donor_data **Data_, size_t Count, int MaxSize) {

  int iDim;
  int iPoint;

  *Data_ = malloc(sizeof(t_donor_data));
  t_donor_data *Data = *Data_;

  Data->count = Count;
  Data->max_size = MaxSize;

  Data->extents[0][0] = malloc(2*MAX_DIMS*Count*sizeof(int));
  Data->extents[1][0] = Data->extents[0][0] + MAX_DIMS*Count;
  for (iDim = 1; iDim < MAX_DIMS; ++iDim) {
    Data->extents[0][iDim] = Data->extents[0][iDim-1] + Count;
    Data->extents[1][iDim] = Data->extents[1][iDim-1] + Count;
  }

  Data->coords[0] = malloc(MAX_DIMS*Count*sizeof(double));
  Data->coords[1] = Data->coords[0] + Count;
  Data->coords[2] = Data->coords[1] + Count;

  Data->interp_coefs[0] = malloc(MAX_DIMS*MaxSize*sizeof(double *));
  Data->interp_coefs[1] = Data->interp_coefs[0] + MaxSize;
  Data->interp_coefs[2] = Data->interp_coefs[1] + MaxSize;
  if (MaxSize > 0) {
    Data->interp_coefs[0][0] = malloc(MAX_DIMS*MaxSize*Count*sizeof(double));
    Data->interp_coefs[1][0] = Data->interp_coefs[0][0] + MaxSize*Count;
    Data->interp_coefs[2][0] = Data->interp_coefs[1][0] + MaxSize*Count;
    for (iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (iPoint = 1; iPoint < MaxSize; ++iPoint) {
        Data->interp_coefs[iDim][iPoint] = Data->interp_coefs[iDim][iPoint-1] + Count;
      }
    }
  }

  Data->destination_grid_ids = malloc(Count*sizeof(int));

  Data->destination_points[0] = malloc(MAX_DIMS*Count*sizeof(int));
  Data->destination_points[1] = Data->destination_points[0] + Count;
  Data->destination_points[2] = Data->destination_points[1] + Count;

}

static void DestroyDonorData(t_donor_data **Data_) {

  t_donor_data *Data = *Data_;

  free(Data->extents[0][0]);
  free(Data->coords[0]);
  if(Data->max_size > 0) {
    free(Data->interp_coefs[0][0]);
  }
  free(Data->interp_coefs[0]);
  free(Data->destination_grid_ids);
  free(Data->destination_points[0]);

  free_null(Data_);

}

static void CreateReceiverData(t_receiver_data **Data_, size_t Count) {

  *Data_ = malloc(sizeof(t_receiver_data));
  t_receiver_data *Data = *Data_;

  Data->count = Count;

  Data->points[0] = malloc(MAX_DIMS*Count*sizeof(int));
  Data->points[1] = Data->points[0] + Count;
  Data->points[2] = Data->points[1] + Count;

  Data->source_grid_ids = malloc(Count*sizeof(int));

  Data->source_cells[0] = malloc(MAX_DIMS*Count*sizeof(int));
  Data->source_cells[1] = Data->source_cells[0] + Count;
  Data->source_cells[2] = Data->source_cells[1] + Count;

}

static void DestroyReceiverData(t_receiver_data **Data_) {

  t_receiver_data *Data = *Data_;

  free(Data->points[0]);
  free(Data->source_grid_ids);
  free(Data->source_cells[0]);

  free_null(Data_);

}

static void CreateConnectionData(t_connection_data **Data_, size_t Count) {

  *Data_ = malloc(sizeof(t_connection_data));
  t_connection_data *Data = *Data_;

  Data->count = Count;

  Data->connection_ids = malloc(Count*sizeof(size_t));

  Data->donor_grid_ids = malloc(Count*sizeof(int));
  Data->donor_cells[0] = malloc(MAX_DIMS*Count*sizeof(int));
  Data->donor_cells[1] = Data->donor_cells[0] + Count;
  Data->donor_cells[2] = Data->donor_cells[1] + Count;

  Data->receiver_grid_ids = malloc(Count*sizeof(int));
  Data->receiver_points[0] = malloc(MAX_DIMS*Count*sizeof(int));
  Data->receiver_points[1] = Data->receiver_points[0] + Count;
  Data->receiver_points[2] = Data->receiver_points[1] + Count;

}

static void DestroyConnectionData(t_connection_data **Data_) {

  t_connection_data *Data = *Data_;

  free(Data->connection_ids);
  free(Data->donor_grid_ids);
  free(Data->donor_cells[0]);
  free(Data->receiver_grid_ids);
  free(Data->receiver_points[0]);

  free_null(Data_);

}

static size_t BinDivide(size_t N, int NumBins) {

  return (N + (size_t)NumBins - 1)/(size_t)NumBins;

}

static void Chunkify(size_t Count, int MaxChunks, size_t TargetChunkSize, int Adjust,
  int *ChunkInterval_, int *NumChunks_, size_t *ChunkSize) {

  int ChunkInterval = 1;
  int NumChunks = MaxChunks;

  while (BinDivide(Count, NumChunks) < TargetChunkSize && NumChunks > 1) {
    ChunkInterval = (int)min((size_t)ChunkInterval << 1, MaxChunks);
    NumChunks = BinDivide((size_t)MaxChunks, ChunkInterval);
  }

  if (Adjust > 0) {
    int RemainingAdjustAmount = Adjust;
    while (RemainingAdjustAmount > 0 && NumChunks < MaxChunks) {
      ChunkInterval = max(ChunkInterval >> 1, 1);
      NumChunks = BinDivide((size_t)MaxChunks, ChunkInterval);
      --RemainingAdjustAmount;
    }
  } else if (Adjust < 0) {
    int RemainingAdjustAmount = -Adjust;
    while (RemainingAdjustAmount > 0 && NumChunks > 1) {
      ChunkInterval = (int)min((size_t)ChunkInterval << 1, MaxChunks);
      NumChunks = BinDivide((size_t)MaxChunks, ChunkInterval);
      --RemainingAdjustAmount;
    }
  }

  *ChunkInterval_ = ChunkInterval;
  *NumChunks_ = NumChunks;
  *ChunkSize = BinDivide(Count, NumChunks);

}

static int File_read_all_endian(MPI_File File, void *Buffer, int Count, MPI_Datatype DataType,
  ovk_ext_endian Endian, MPI_Status *Status, t_profiler *Profiler, MPI_Comm Comm) {

  int MPIIOReadTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Read");
  StartProfileSync(Profiler, MPIIOReadTime, Comm);
  int MPIError = MPI_File_read_all(File, Buffer, Count, DataType, Status);
  EndProfile(Profiler, MPIIOReadTime);

  if (MPIError == MPI_SUCCESS) {
    if (Endian != MachineEndian()) {
      int Size;
      MPI_Type_size(DataType, &Size);
      SwapEndian(Buffer, Size, Count);
    }
  }

  return MPIError;

}

static int File_read_at_endian(MPI_File File, MPI_Offset Offset, void *Buffer, int Count,
  MPI_Datatype DataType, ovk_ext_endian Endian, MPI_Status *Status, t_profiler *Profiler) {

  int MPIIOReadTime = GetProfilerTimerID(Profiler, "XINTOUT::Read::MPI-IO::Read");
  StartProfile(Profiler, MPIIOReadTime);
  int MPIError = MPI_File_read_at(File, Offset, Buffer, Count, DataType, Status);
  EndProfile(Profiler, MPIIOReadTime);

  if (MPIError == MPI_SUCCESS) {
    if (Endian != MachineEndian()) {
      int Size;
      MPI_Type_size(DataType, &Size);
      SwapEndian(Buffer, Size, Count);
    }
  }

  return MPIError;

}

static ovk_ext_endian MachineEndian() {

  unsigned char EndianTest[2] = {1, 0};

  if(*(short *)EndianTest == 1) {
    return OVK_EXT_LITTLE_ENDIAN;
  } else {
    return OVK_EXT_BIG_ENDIAN;
  }

}

static void SwapEndian(void *Data, int ElementSize, int NumElements) {

  int i, j;

  char *DataBytes = Data;

  char *Element = malloc(ElementSize);

  for (i = 0; i < NumElements; ++i) {
    for (j = 0; j < ElementSize; ++j) {
      Element[j] = DataBytes[ElementSize*i+j];
    }
    for (j = 0; j < ElementSize; ++j) {
      DataBytes[ElementSize*i+j] = Element[ElementSize-j-1];
    }
  }

  free(Element);

}
