// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

// Note: Current implementation will not play well with large datasets on small numbers of ranks
// (MPI I/O and send/recv routines will fail if data has count larger than INT_MAX). Can be fixed
// by splitting into multiple calls of size INT_MAX.

#include "ovk/extras/XINTOUT.hpp"

#include "ovk/extras/Constants.hpp"
#include "ovk/extras/Global.hpp"
#include "ovk/core/Array.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/ErrorHandler.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/ScopeGuard.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace ovk {

namespace {

struct donor_data {
  long long Count;
  int MaxSize;
  array<int,3> Extents;
  array<double,2> Coords;
  array<double,3> InterpCoefs;
  array<int> DestinationGridIDs;
  array<int,2> DestinationPoints;
};

struct receiver_data {
  long long Count;
  array<int,2> Points;
  array<int> SourceGridIDs;
  array<int,2> SourceCells;
};

struct connection_data {
  long long Count;
  array<long long> ConnectionIDs;
  array<int> DonorGridIDs;
  array<int,2> DonorCells;
  array<int> ReceiverGridIDs;
  array<int,2> ReceiverPoints;
};

struct xintout_donor_chunk {
  long long Begin;
  long long End;
  donor_data Data;
  long long StartingConnectionID;
};

struct xintout_receiver_chunk {
  long long Begin;
  long long End;
  receiver_data Data;
  array<long long> ConnectionIDs;
};

struct xintout_donors {
  long long Count;
  long long ChunkSize;
  bool HasChunk;
  xintout_donor_chunk Chunk;
};

struct xintout_receivers {
  long long Count;
  long long ChunkSize;
  bool HasChunk;
  xintout_receiver_chunk Chunk;
};

struct xintout_grid {
  mutable core::logger *Logger;
  mutable core::error_handler *ErrorHandler;
  int ID;
  std::string Name;
  int NumDims;
  core::comm Comm;
  elem<int,MAX_DIMS> GlobalSize;
  xintout_donors Donors;
  xintout_receivers Receivers;
};

struct xintout_connection_bin {
  long long Begin;
  long long End;
  connection_data Data;
};

struct xintout_connections {
  long long Count;
  long long BinSize;
  bool HasBin;
  xintout_connection_bin Bin;
};

struct xintout {
  mutable core::logger *Logger;
  mutable core::error_handler *ErrorHandler;
  int NumDims;
  core::comm Comm;
  int NumGrids;
  int NumLocalGrids;
  std::vector<xintout_grid> Grids;
  xintout_connections Connections;
};

void CreateXINTOUT(xintout &XINTOUT, int NumDims, core::comm Comm, int NumGrids, int NumLocalGrids,
  const std::vector<int> &LocalGridIDs, const std::vector<std::string> &LocalGridNames,
  const std::vector<core::comm> &LocalGridComms, const std::vector<elem<int,MAX_DIMS>>
  &LocalGridGlobalSizes, core::logger &Logger, core::error_handler &ErrorHandler);

void CreateXINTOUTGrid(xintout_grid &XINTOUTGrid, int ID, const std::string &Name, int NumDims,
  core::comm Comm, const elem<int,MAX_DIMS> &GlobalSize, core::logger &Logger,
  core::error_handler &ErrorHandler);

error ReadXINTOUT(xintout &XINTOUT, const std::string &HOPath, const std::string &XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo, core::profiler &Profiler);

error ReadGlobalInfo(const xintout &XINTOUT, const std::string &HOPath,
  endian &Endian, xintout_format &Format, bool &WithIBlank, core::profiler &Profiler);

bool DetectFormat(MPI_File HOFile, endian &Endian, xintout_format &Format, core::profiler
  &Profiler);

error ReadGridInfo(const xintout_grid &XINTOUTGrid, const std::string &HOPath,
  const std::string &XPath, long long &NumDonors, long long &NumReceivers, long long
  &StartingConnectionID, MPI_Offset &HODonorCellsOffset, MPI_Offset &HODonorCoordsOffset,
  MPI_Offset &HOReceiverPointsOffset, MPI_Offset &HOReceiverConnectionIDsOffset,
  MPI_Offset &XDonorSizesOffset, MPI_Offset &XDonorInterpCoefsOffset, endian Endian,
  xintout_format Format, bool WithIBlank, core::profiler &Profiler);

error ReadDonors(xintout_grid &XINTOUTGrid, const std::string &HOPath, const std::string &XPath,
  long long NumDonors, long long StartingConnectionID, MPI_Offset HOCellsOffset,
  MPI_Offset HOCoordsOffset, MPI_Offset XSizesOffset, MPI_Offset XInterpCoefsOffset,
  endian Endian, int ReadGranularityAdjust, MPI_Info MPIInfo, core::profiler &Profiler);

error ReadReceivers(xintout_grid &XINTOUTGrid, const std::string &HOPath, long long NumReceivers,
  MPI_Offset HOPointsOffset, MPI_Offset HOConnectionIDsOffset, endian Endian, xintout_format Format,
  int ReadGranularityAdjust, MPI_Info MPIInfo, core::profiler &Profiler);

void MatchDonorsAndReceivers(xintout &XINTOUT, core::profiler &Profiler);

void DistributeConnectivityData(const xintout &XINTOUT, const std::vector<const grid *> &LocalGrids,
  std::vector<donor_data> &LocalDonorData, std::vector<receiver_data> &LocalReceiverData,
  core::profiler &Profiler);

void DistributeGridConnectivityData(const xintout_grid &XINTOUTGrid, const grid &Grid,
  donor_data &DonorData, receiver_data &ReceiverData, core::profiler &Profiler);

void ImportConnectivityData(int NumGrids, int NumLocalGrids, const std::vector<int> &LocalGridIDs,
  const std::vector<donor_data> &LocalDonorData, const std::vector<receiver_data>
  &LocalReceiverData, domain &Domain);

void ImportDonors(const donor_data &GridDonors, const core::comm &Comm, int NumReceiverGrids, const
  std::vector<connectivity_d *> &OverkitDonors);
void ImportReceivers(const receiver_data &GridReceivers, const core::comm &Comm, int NumDonorGrids,
  const std::vector<connectivity_r *> &OverkitReceivers);

void CreateDonorData(donor_data &Data, long long Count, int MaxSize);
void CreateReceiverData(receiver_data &Data, long long Count);
void CreateConnectionData(connection_data &Data, long long Count);

long long BinDivide(long long N, int NumChunks);
void Chunkify(long long Count, int MaxChunks, long long TargetChunkSize, int Adjust,
  int &ChunkRankInterval, int &NumChunks, long long &ChunkSize);

int File_read_all_endian(MPI_File File, void *Buffer, int Count, MPI_Datatype DataType,
  endian Endian, MPI_Status *Status, core::profiler &Profiler, const core::comm &Comm);
int File_read_at_endian(MPI_File File, MPI_Offset Offset, void *Buffer, int Count,
  MPI_Datatype DataType, endian Endian, MPI_Status *Status, core::profiler &Profiler);

endian MachineEndian();
void SwapEndian(void *Data, int ElementSize, int NumElements);

}

error ImportXINTOUT(domain &Domain, const std::string &HOPath, const std::string &XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo) {

  core::logger &Logger = core::GetDomainLogger(Domain);
  core::error_handler &ErrorHandler = core::GetDomainErrorHandler(Domain);

  int NumDims;
  int NumGrids;
  GetDomainDimension(Domain, NumDims);
  GetDomainGridCount(Domain, NumGrids);

  const core::comm &Comm = core::GetDomainComm(Domain);

  core::profiler &Profiler = core::GetDomainProfiler(Domain);
  core::AddProfilerTimer(Profiler, "XINTOUT::Import");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Read");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Read::MPI-IO::Open");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Read::MPI-IO::Close");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Read::MPI-IO::Read");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Read::MPI-IO::Other");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Match");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Match::MapToBins");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Match::Handshake");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Match::SendToBins");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Match::FillConnectionData");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Match::RecvFromBins");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Match::Unpack");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Distribute");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Distribute::MapToBins");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Distribute::RetrieveBins");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Distribute::FindRanks");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Distribute::Handshake");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::Distribute::SendData");
  core::AddProfilerTimer(Profiler, "XINTOUT::Import::SetConnectivities");

  MPI_Barrier(Comm);

  if (NumGrids > 0) {

    int NumLocalGrids = 0;
    for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
      int GridID = iGrid+1;
      if (RankHasGrid(Domain, GridID)) {
        ++NumLocalGrids;
      }
    }

    std::vector<const grid *> LocalGrids(NumLocalGrids);
    std::vector<int> LocalGridIDs(NumLocalGrids);
    std::vector<std::string> LocalGridNames(NumLocalGrids);
    std::vector<core::comm> LocalGridComms(NumLocalGrids);
    std::vector<elem<int,MAX_DIMS>> LocalGridGlobalSizes(NumLocalGrids);
    int iLocalGrid = 0;
    for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
      int GridID = iGrid+1;
      if (RankHasGrid(Domain, GridID)) {
        const grid *GridPtr;
        GetGrid(Domain, GridID, GridPtr);
        const grid &Grid = *GridPtr;
        std::string Name;
        elem<int,MAX_DIMS> GlobalSize;
        GetGridName(Grid, Name);
        GetGridSize(Grid, GlobalSize.Data());
        core::comm GridComm = core::GetGridComm(Grid);
        LocalGrids[iLocalGrid] = GridPtr;
        LocalGridIDs[iLocalGrid] = GridID;
        LocalGridNames[iLocalGrid] = Name;
        LocalGridComms[iLocalGrid] = std::move(GridComm);
        LocalGridGlobalSizes[iLocalGrid] = GlobalSize;
        ++iLocalGrid;
      }
    }

    int ImportTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import");
    int ReadTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read");
    int MatchTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Match");
    int DistributeTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Distribute");
    int SetConnectivitiesTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::SetConnectivities");

    core::StartProfileSync(Profiler, ImportTime, Comm);

    xintout XINTOUT;
    CreateXINTOUT(XINTOUT, NumDims, Comm, NumGrids, NumLocalGrids, LocalGridIDs,
      LocalGridNames, LocalGridComms, LocalGridGlobalSizes, Logger, ErrorHandler);

    core::StartProfileSync(Profiler, ReadTime, Comm);
    error Error = ReadXINTOUT(XINTOUT, HOPath, XPath, ReadGranularityAdjust, MPIInfo, Profiler);
    core::EndProfile(Profiler, ReadTime);
    OVK_EH_CHECK(ErrorHandler, Error);

    core::StartProfileSync(Profiler, MatchTime, Comm);
    MatchDonorsAndReceivers(XINTOUT, Profiler);
    core::EndProfile(Profiler, MatchTime);

    core::StartProfileSync(Profiler, DistributeTime, Comm);
    std::vector<donor_data> LocalDonors(NumLocalGrids);
    std::vector<receiver_data> LocalReceivers(NumLocalGrids);
    DistributeConnectivityData(XINTOUT, LocalGrids, LocalDonors, LocalReceivers, Profiler);
    core::EndProfile(Profiler, DistributeTime);

    core::StartProfileSync(Profiler, SetConnectivitiesTime, Comm);
    ImportConnectivityData(NumGrids, NumLocalGrids, LocalGridIDs, LocalDonors, LocalReceivers,
      Domain);
    core::EndProfile(Profiler, SetConnectivitiesTime);

    core::EndProfile(Profiler, ImportTime);

  }

  return error::NONE;

}

error ExportXINTOUT(const domain &Domain, const std::string &HOPath, const std::string &XPath,
  xintout_format Format, endian Endian, int WriteGranularityAdjust, MPI_Info MPIInfo) {

  OVK_DEBUG_ASSERT(false, "ExportXINTOUT is not yet implemented.");

  return error::NONE;

}

namespace {

void CreateXINTOUT(xintout &XINTOUT, int NumDims, core::comm Comm, int NumGrids, int NumLocalGrids,
  const std::vector<int> &LocalGridIDs, const std::vector<std::string> &LocalGridNames,
  const std::vector<core::comm> &LocalGridComms, const std::vector<elem<int,MAX_DIMS>>
  &LocalGridGlobalSizes, core::logger &Logger, core::error_handler &ErrorHandler) {

  XINTOUT.Comm = std::move(Comm);

  MPI_Barrier(XINTOUT.Comm);

  XINTOUT.NumDims = NumDims;

  XINTOUT.Logger = &Logger;
  XINTOUT.ErrorHandler = &ErrorHandler;

  XINTOUT.NumGrids = NumGrids;
  XINTOUT.NumLocalGrids = NumLocalGrids;

  XINTOUT.Grids.resize(NumLocalGrids);
  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    CreateXINTOUTGrid(XINTOUT.Grids[iLocalGrid], LocalGridIDs[iLocalGrid],
      LocalGridNames[iLocalGrid], NumDims, LocalGridComms[iLocalGrid],
      LocalGridGlobalSizes[iLocalGrid], Logger, ErrorHandler);
  }

  XINTOUT.Connections.Count = 0;
  XINTOUT.Connections.BinSize = 0;
  XINTOUT.Connections.HasBin = false;

  MPI_Barrier(XINTOUT.Comm);

}

void CreateXINTOUTGrid(xintout_grid &XINTOUTGrid, int ID, const std::string &Name, int NumDims,
  core::comm Comm, const elem<int,MAX_DIMS> &GlobalSize, core::logger &Logger,
  core::error_handler &ErrorHandler) {

  XINTOUTGrid.Comm = std::move(Comm);

  MPI_Barrier(XINTOUTGrid.Comm);

  XINTOUTGrid.Logger = &Logger;
  XINTOUTGrid.ErrorHandler = &ErrorHandler;

  XINTOUTGrid.ID = ID;
  XINTOUTGrid.Name = Name;
  XINTOUTGrid.NumDims = NumDims;
  XINTOUTGrid.GlobalSize = GlobalSize;

  XINTOUTGrid.Donors.Count = 0;
  XINTOUTGrid.Donors.ChunkSize = 0;
  XINTOUTGrid.Donors.HasChunk = false;

  XINTOUTGrid.Receivers.Count = 0;
  XINTOUTGrid.Receivers.ChunkSize = 0;
  XINTOUTGrid.Receivers.HasChunk = false;

  MPI_Barrier(XINTOUTGrid.Comm);

}

error ReadXINTOUT(xintout &XINTOUT, const std::string &HOPath, const std::string &XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo, core::profiler &Profiler) {

  const core::comm &Comm = XINTOUT.Comm;

  MPI_Barrier(Comm);

  core::error_handler &ErrorHandler = *XINTOUT.ErrorHandler;
  core::logger &Logger = *XINTOUT.Logger;

  int NumLocalGrids = XINTOUT.NumLocalGrids;

  core::LogStatus(Logger, Comm.Rank() == 0, 0, "Reading XINTOUT files '%s' and '%s'...", HOPath,
    XPath);

  error Error = error::NONE;

  endian Endian;
  xintout_format Format;
  bool WithIBlank;
  Error = ReadGlobalInfo(XINTOUT, HOPath, Endian, Format, WithIBlank, Profiler);
  OVK_EH_CHECK(ErrorHandler, Error);

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {

    xintout_grid &XINTOUTGrid = XINTOUT.Grids[iLocalGrid];

    long long NumDonors, NumReceivers;
    long long DonorStartingConnectionID;
    MPI_Offset HODonorCellsOffset, HODonorCoordsOffset, XDonorSizesOffset, XDonorInterpCoefsOffset;
    MPI_Offset HOReceiverPointsOffset, HOReceiverConnectionIDsOffset;
    Error = ReadGridInfo(XINTOUTGrid, HOPath, XPath, NumDonors, NumReceivers,
      DonorStartingConnectionID, HODonorCellsOffset, HODonorCoordsOffset, HOReceiverPointsOffset,
      HOReceiverConnectionIDsOffset, XDonorSizesOffset, XDonorInterpCoefsOffset, Endian,
      Format, WithIBlank, Profiler);
    OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);

    Error = ReadDonors(XINTOUTGrid, HOPath, XPath, NumDonors, DonorStartingConnectionID,
      HODonorCellsOffset, HODonorCoordsOffset, XDonorSizesOffset, XDonorInterpCoefsOffset, Endian,
      ReadGranularityAdjust, MPIInfo, Profiler);
    OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);

    Error = ReadReceivers(XINTOUTGrid, HOPath, NumReceivers, HOReceiverPointsOffset,
      HOReceiverConnectionIDsOffset, Endian, Format, ReadGranularityAdjust, MPIInfo, Profiler);
    OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);

  }

  done_reading:
    OVK_EH_SYNC(ErrorHandler, Error, Comm);
    OVK_EH_CHECK(ErrorHandler, Error);

  MPI_Barrier(Comm);

  core::LogStatus(Logger, Comm.Rank() == 0, 0, "Finished reading XINTOUT files.");

  return error::NONE;

}

error ReadGlobalInfo(const xintout &XINTOUT, const std::string &HOPath, endian &Endian,
  xintout_format &Format, bool &WithIBlank, core::profiler &Profiler) {

  const core::comm &Comm = XINTOUT.Comm;

  int NumGrids = XINTOUT.NumGrids;

  core::logger &Logger = *XINTOUT.Logger;
  core::error_handler &ErrorHandler = *XINTOUT.ErrorHandler;

  int MPIIOOpenTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Open");
  int MPIIOCloseTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Close");

  error Error = error::NONE;

  if (Comm.Rank() == 0) {

    MPI_File HOFile;
    MPI_Status Status;
    int MPIError;
    int ReadSize;

    core::StartProfile(Profiler, MPIIOOpenTime);
    // MPI_File_open missing const qualifier for path string on some platforms
    MPIError = MPI_File_open(MPI_COMM_SELF, const_cast<char *>(HOPath.c_str()), MPI_MODE_RDONLY,
      MPI_INFO_NULL, &HOFile);
    core::EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      Error = error::FILE_OPEN;
      core::LogError(Logger, true, "Unable to open file '%s'.", HOPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }
    auto CloseHO = core::OnScopeExit([&]() {
      core::StartProfile(Profiler, MPIIOCloseTime);
      MPI_File_close(&HOFile);
      core::EndProfile(Profiler, MPIIOCloseTime);
    });

    bool Success = DetectFormat(HOFile, Endian, Format, Profiler);
    if (!Success) {
      Error = error::FILE_READ;
      core::LogError(Logger, true, "Unable to detect format of XINTOUT file '%s'.", HOPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }

    int RecordWrapperSize = Format == xintout_format::STANDARD ? sizeof(int) : 0;

    // Both formats have a record wrapper at the beginning
    MPI_Offset HOHeaderOffset = sizeof(int);
    int NumGridsInFile;
    File_read_at_endian(HOFile, HOHeaderOffset, &NumGridsInFile, 1, MPI_INT, Endian, &Status,
      Profiler);
    MPI_Get_count(&Status, MPI_INT, &ReadSize);
    if (ReadSize < 1) {
      Error = error::FILE_READ;
      core::LogError(Logger, true, "Unable to read header of XINTOUT file '%s'.", HOPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }
    if (NumGridsInFile != NumGrids) {
      Error = error::FILE_READ;
      core::LogError(Logger, Comm.Rank() == 0, "XINTOUT file '%s' has incorrect number of grids.",
        HOPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }

    MPI_Offset HOGridOffset = 0;

    HOGridOffset += RecordWrapperSize;
    if (Format == xintout_format::STANDARD) {
      HOGridOffset += 5*sizeof(int);
    } else {
      HOGridOffset += sizeof(int) + 4*sizeof(long long);
    }
    HOGridOffset += RecordWrapperSize;

    HOGridOffset += RecordWrapperSize;
    long long NumDonors;
    long long NumReceivers;
    elem<int,MAX_DIMS> GridSize;
    if (Format == xintout_format::STANDARD) {
      int Data[7];
      File_read_at_endian(HOFile, HOGridOffset, Data, 7, MPI_INT, Endian, &Status, Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < 7) {
        Error = error::FILE_READ;
        core::LogError(Logger, true, "Unable to read grid 1 header of XINTOUT file '%s'.", HOPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      NumDonors = Data[1];
      NumReceivers = Data[0];
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
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
        core::LogError(Logger, true, "Unable to read grid 1 header of XINTOUT file '%s'.", HOPath);
        Error = error::FILE_READ;
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      HOGridOffset += 4*sizeof(long long);
      File_read_at_endian(HOFile, HOGridOffset, IntData, 3, MPI_INT, Endian, &Status, Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < 3) {
        Error = error::FILE_READ;
        core::LogError(Logger, true, "Unable to read grid 1 header of XINTOUT file '%s'.", HOPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      HOGridOffset += 3*sizeof(int);
      NumDonors = LongLongData[1];
      NumReceivers = LongLongData[0];
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        GridSize[iDim] = IntData[iDim];
      }
      HOGridOffset += sizeof(long long);
    }
    HOGridOffset += RecordWrapperSize;

    HOGridOffset += RecordWrapperSize;
    HOGridOffset += MAX_DIMS*NumDonors*sizeof(int);
    HOGridOffset += MAX_DIMS*NumDonors*sizeof(double);
    HOGridOffset += RecordWrapperSize;

    int ConnectionIDSize = Format == xintout_format::STANDARD ? sizeof(int) : sizeof(long long);
    HOGridOffset += RecordWrapperSize;
    HOGridOffset += MAX_DIMS*NumReceivers*sizeof(int);
    HOGridOffset += NumReceivers*ConnectionIDSize;
    HOGridOffset += RecordWrapperSize;

    if (Format == xintout_format::STANDARD) {
      int RecordWrapper;
      File_read_at_endian(HOFile, HOGridOffset, &RecordWrapper, 1, MPI_INT, Endian, &Status,
        Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize == 1) {
        long long NumPoints =
          (long long)(GridSize[0])*
          (long long)(GridSize[1])*
          (long long)(GridSize[2]);
        if ((long long)(RecordWrapper) == NumPoints*int(sizeof(int))) {
          WithIBlank = true;
        } else {
          WithIBlank = false;
        }
      } else {
        if (NumGrids == 1) {
          WithIBlank = false;
        } else {
          Error = error::FILE_READ;
          core::LogError(Logger, true, "Unable to detect whether XINTOUT file '%s' contains IBlank.",
            HOPath);
          OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
        }
      }
    } else {
      WithIBlank = false;
    }

  }

  done_reading:
    OVK_EH_SYNC(ErrorHandler, Error, Comm);
    OVK_EH_CHECK(ErrorHandler, Error);

  int EndianInt;
  if (Comm.Rank() == 0) EndianInt = Endian == endian::BIG ? 1 : 0;
  MPI_Bcast(&EndianInt, 1, MPI_INT, 0, Comm);
  Endian = EndianInt == 1 ? endian::BIG : endian::LITTLE;

  int FormatInt;
  if (Comm.Rank() == 0) FormatInt = Format == xintout_format::EXTENDED ? 1 : 0;
  MPI_Bcast(&FormatInt, 1, MPI_INT, 0, Comm);
  Format = FormatInt == 1 ? xintout_format::EXTENDED : xintout_format::STANDARD;

  int WithIBlankInt;
  if (Comm.Rank() == 0) WithIBlankInt = int(WithIBlank);
  MPI_Bcast(&WithIBlankInt, 1, MPI_INT, 0, Comm);
  WithIBlank = WithIBlankInt != 0;

  return error::NONE;

}

bool DetectFormat(MPI_File HOFile, endian &Endian, xintout_format &Format, core::profiler &Profiler) {

  int MPIIOReadTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Read");

  unsigned char InitialBytes[sizeof(int)];
  MPI_Status Status;
  core::StartProfile(Profiler, MPIIOReadTime);
  int MPIError = MPI_File_read_at(HOFile, 0, InitialBytes, sizeof(int), MPI_BYTE, &Status);
  core::EndProfile(Profiler, MPIIOReadTime);
  if (MPIError != MPI_SUCCESS) return false;

  // If little endian, the first byte will be the size of the file header data
  if (InitialBytes[0] != 0) {
    Endian = endian::LITTLE;
  } else {
    Endian = endian::BIG;
  }

  int HeaderSize;
  memcpy(&HeaderSize, InitialBytes, sizeof(int));
  if (Endian != MachineEndian()) {
    SwapEndian(&HeaderSize, sizeof(int), 1);
  }

  if (HeaderSize == 5*sizeof(int)) {
    Format = xintout_format::STANDARD;
  } else if (HeaderSize == sizeof(int) + 4*sizeof(long long)) {
    Format = xintout_format::EXTENDED;
  } else {
    return false;
  }

  return true;

}

error ReadGridInfo(const xintout_grid &XINTOUTGrid, const std::string &HOPath,
  const std::string &XPath, long long &NumDonors, long long &NumReceivers, long long
  &StartingConnectionID, MPI_Offset &HODonorCellsOffset, MPI_Offset &HODonorCoordsOffset,
  MPI_Offset &HOReceiverPointsOffset, MPI_Offset &HOReceiverConnectionIDsOffset,
  MPI_Offset &XDonorSizesOffset, MPI_Offset &XDonorInterpCoefsOffset, endian Endian,
  xintout_format Format, bool WithIBlank, core::profiler &Profiler) {

  int GridID = XINTOUTGrid.ID;
  int iGrid = GridID-1;

  const core::comm &Comm = XINTOUTGrid.Comm;

  core::logger &Logger = *XINTOUTGrid.Logger;
  core::error_handler &ErrorHandler = *XINTOUTGrid.ErrorHandler;

  int MPIIOOpenTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Open");
  int MPIIOCloseTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Close");

  error Error = error::NONE;

  if (Comm.Rank() == 0) {

    MPI_File HOFile, XFile;
    MPI_Status Status;
    int MPIError;
    int ReadSize;
    elem<int,MAX_DIMS> GridSize;
    long long NumInterpCoefs = 0;

    core::StartProfile(Profiler, MPIIOOpenTime);
    // MPI_File_open missing const qualifier for path string on some platforms
    MPIError = MPI_File_open(MPI_COMM_SELF, const_cast<char *>(HOPath.c_str()), MPI_MODE_RDONLY,
      MPI_INFO_NULL, &HOFile);
    core::EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      Error = error::FILE_OPEN;
      core::LogError(Logger, true, "Unable to open file '%s'.", HOPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }
    auto CloseHO = core::OnScopeExit([&]() {
      core::StartProfile(Profiler, MPIIOCloseTime);
      MPI_File_close(&HOFile);
      core::EndProfile(Profiler, MPIIOCloseTime);
    });

    core::StartProfile(Profiler, MPIIOOpenTime);
    // MPI_File_open missing const qualifier for path string on some platforms
    MPIError = MPI_File_open(MPI_COMM_SELF, const_cast<char *>(XPath.c_str()), MPI_MODE_RDONLY,
      MPI_INFO_NULL, &XFile);
    core::EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      Error = error::FILE_OPEN;
      core::LogError(Logger, true, "Unable to open file '%s'.", XPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }
    auto CloseX = core::OnScopeExit([&]() {
      core::StartProfile(Profiler, MPIIOCloseTime);
      MPI_File_close(&XFile);
      core::EndProfile(Profiler, MPIIOCloseTime);
    });

    int RecordWrapperSize = Format == xintout_format::STANDARD ? sizeof(int) : 0;

    MPI_Offset HOGridOffset = 0;
    MPI_Offset XGridOffset = 0;

    int HeaderSize;
    if (Format == xintout_format::STANDARD) {
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

    for (int iOtherGrid = 0; iOtherGrid < iGrid; ++iOtherGrid) {

      int OtherGridID = iOtherGrid+1;

      HOGridOffset += RecordWrapperSize;
      if (Format == xintout_format::STANDARD) {
        int Data[7];
        File_read_at_endian(HOFile, HOGridOffset, Data, 7, MPI_INT, Endian, &Status, Profiler);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        if (ReadSize < 7) {
          Error = error::FILE_READ;
          core::LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
            HOPath);
          OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
        }
        NumDonors = Data[1];
        NumReceivers = Data[0];
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
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
          Error = error::FILE_READ;
          core::LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
            HOPath);
          OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
        }
        HOGridOffset += 4*sizeof(long long);
        File_read_at_endian(HOFile, HOGridOffset, IntData, 3, MPI_INT, Endian, &Status, Profiler);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        if (ReadSize < 3) {
          Error = error::FILE_READ;
          core::LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
            HOPath);
          OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
        }
        HOGridOffset += 3*sizeof(int);
        File_read_at_endian(HOFile, HOGridOffset, LongLongData+4, 1, MPI_LONG_LONG, Endian, &Status,
          Profiler);
        MPI_Get_count(&Status, MPI_LONG_LONG, &ReadSize);
        if (ReadSize < 1) {
          Error = error::FILE_READ;
          core::LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
            HOPath);
          OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
        }
        HOGridOffset += sizeof(long long);
        NumDonors = LongLongData[1];
        NumReceivers = LongLongData[0];
        NumInterpCoefs = LongLongData[4];
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          GridSize[iDim] = IntData[iDim];
        }
      }
      HOGridOffset += RecordWrapperSize;

      HOGridOffset += RecordWrapperSize;
      HOGridOffset += MAX_DIMS*NumDonors*sizeof(int);
      HOGridOffset += MAX_DIMS*NumDonors*sizeof(double);
      HOGridOffset += RecordWrapperSize;

      int ConnectionIDSize = Format == xintout_format::STANDARD ? sizeof(int) : sizeof(long long);
      HOGridOffset += RecordWrapperSize;
      HOGridOffset += MAX_DIMS*NumReceivers*sizeof(int);
      HOGridOffset += NumReceivers*ConnectionIDSize;
      HOGridOffset += RecordWrapperSize;

      if (WithIBlank) {
        HOGridOffset += RecordWrapperSize;
        HOGridOffset +=
          (long long)(GridSize[0])*
          (long long)(GridSize[1])*
          (long long)(GridSize[2])*
          sizeof(int);
        HOGridOffset += RecordWrapperSize;
      }

      XGridOffset += RecordWrapperSize;
      XGridOffset += MAX_DIMS*NumDonors*sizeof(int);
      XGridOffset += RecordWrapperSize;

      if (Format == xintout_format::STANDARD) {
        // Figure out size of interp coef data
        int RecordWrapper;
        File_read_at_endian(XFile, XGridOffset, &RecordWrapper, 1, MPI_INT, Endian, &Status,
          Profiler);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        if (ReadSize < 1) {
          Error = error::FILE_READ;
          core::LogError(Logger, true, "Unable to read grid %i data from XINTOUT file '%s'.",
            OtherGridID, XPath);
          OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
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
    if (Format == xintout_format::STANDARD) {
      int Data[7];
      File_read_at_endian(HOFile, HOGridOffset, Data, 7, MPI_INT, Endian, &Status, Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < 7) {
        Error = error::FILE_READ;
        core::LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
          HOPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      NumDonors = Data[1];
      NumReceivers = Data[0];
      // Convert to zero-based indexing
      StartingConnectionID = Data[3]-1;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
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
        Error = error::FILE_READ;
        core::LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
          HOPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      HOGridOffset += 4*sizeof(long long);
      File_read_at_endian(HOFile, HOGridOffset, IntData, 3, MPI_INT, Endian, &Status, Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < 3) {
        Error = error::FILE_READ;
        core::LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
          HOPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      HOGridOffset += 3*sizeof(int);
      File_read_at_endian(HOFile, HOGridOffset, LongLongData+4, 1, MPI_LONG_LONG, Endian, &Status,
        Profiler);
      MPI_Get_count(&Status, MPI_LONG_LONG, &ReadSize);
      if (ReadSize < 1) {
        Error = error::FILE_READ;
        core::LogError(Logger, true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
          HOPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      HOGridOffset += sizeof(long long);
      NumDonors = LongLongData[1];
      NumReceivers = LongLongData[0];
      // Convert to zero-based indexing
      StartingConnectionID = LongLongData[3]-1;
      NumInterpCoefs = LongLongData[4];
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        GridSize[iDim] = IntData[iDim];
      }
    }
    HOGridOffset += RecordWrapperSize;

    if (GridSize[0] != XINTOUTGrid.GlobalSize[0] || GridSize[1] != XINTOUTGrid.GlobalSize[1] ||
      GridSize[2] != XINTOUTGrid.GlobalSize[2]) {
      Error = error::FILE_READ;
      core::LogError(Logger, true, "Grid %i of XINTOUT file '%s' has incorrect size.", GridID,
        HOPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
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

  }

  done_reading:
    OVK_EH_SYNC(ErrorHandler, Error, Comm);
    OVK_EH_CHECK(ErrorHandler, Error);

  MPI_Bcast(&NumDonors, 1, MPI_LONG_LONG, 0, Comm);
  MPI_Bcast(&NumReceivers, 1, MPI_LONG_LONG, 0, Comm);
  MPI_Bcast(&StartingConnectionID, 1, MPI_LONG_LONG, 0, Comm);
  MPI_Bcast(&HODonorCellsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&HODonorCoordsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&HOReceiverPointsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&HOReceiverConnectionIDsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&XDonorSizesOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&XDonorInterpCoefsOffset, 1, MPI_OFFSET, 0, Comm);

  return error::NONE;

}

error ReadDonors(xintout_grid &XINTOUTGrid, const std::string &HOPath, const std::string &XPath,
  long long NumDonors, long long StartingConnectionID, MPI_Offset HOCellsOffset,
  MPI_Offset HOCoordsOffset, MPI_Offset XSizesOffset, MPI_Offset XInterpCoefsOffset,
  endian Endian, int ReadGranularityAdjust, MPI_Info MPIInfo, core::profiler &Profiler) {

  int GridID = XINTOUTGrid.ID;

  const core::comm &Comm = XINTOUTGrid.Comm;

  core::logger &Logger = *XINTOUTGrid.Logger;
  core::error_handler &ErrorHandler = *XINTOUTGrid.ErrorHandler;

  xintout_donors &XINTOUTDonors = XINTOUTGrid.Donors;

  XINTOUTDonors.Count = NumDonors;

  // Read chunk on every Nth rank, where N is a power of 2 (or Comm.Size()) chosen such that the
  // number of donors per chunk is roughly the average number of grid points per rank (subject to
  // user adjustment)
  long long NumPoints =
    (long long)(XINTOUTGrid.GlobalSize[0])*
    (long long)(XINTOUTGrid.GlobalSize[1])*
    (long long)(XINTOUTGrid.GlobalSize[2]);
  long long AvgPointsPerRank = BinDivide(NumPoints, Comm.Size());
  int ChunkRankInterval, NumChunks;
  long long ChunkSize;
  Chunkify(NumDonors, Comm.Size(), AvgPointsPerRank, ReadGranularityAdjust, ChunkRankInterval,
    NumChunks, ChunkSize);
  bool HasChunk = Comm.Rank() % ChunkRankInterval == 0;

  if (LoggingStatus(Logger)) {
    std::string NumDonorsString = core::FormatNumber(NumDonors, "donors", "donor");
    std::string NumRanksString = core::FormatNumber(NumChunks, "I/O ranks", "I/O rank");
    core::LogStatus(Logger, Comm.Rank() == 0, 1, "Grid %s has %s; using %s.",
      XINTOUTGrid.Name, NumDonorsString, NumRanksString);
  }

  XINTOUTDonors.ChunkSize = ChunkSize;
  XINTOUTDonors.HasChunk = HasChunk;

  core::comm ChunkComm = core::CreateSubsetComm(Comm, HasChunk);

  error Error = error::NONE;

  if (HasChunk) {

    xintout_donor_chunk &Chunk = XINTOUTDonors.Chunk;

    long long LocalBegin = ChunkSize*ChunkComm.Rank();
    long long LocalEnd = std::min(ChunkSize*(ChunkComm.Rank()+1), NumDonors);
    long long NumLocalDonors = LocalEnd - LocalBegin;

    if (NumLocalDonors*sizeof(double) > std::numeric_limits<int>::max()) {
      Error = error::FILE_READ;
    }
    OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
    if (Error != error::NONE) {
      core::LogError(Logger, ChunkComm.Rank() == 0, "Donor chunk size too big; increase number of "
        "processes or read granularity.");
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }

    Chunk.Begin = LocalBegin;
    Chunk.End = LocalEnd;
    Chunk.StartingConnectionID = StartingConnectionID + LocalBegin;

    MPI_File HOFile, XFile;
    MPI_Status Status;
    int MPIError;
    MPI_Offset DatasetOffset, ReadOffset;
    int ReadSize;

    int MPIIOOpenTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Open");
    int MPIIOCloseTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Close");
    int MPIIOOtherTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Other");

    core::StartProfileSync(Profiler, MPIIOOpenTime, ChunkComm);
    // MPI_File_open missing const qualifier for path string on some platforms
    MPIError = MPI_File_open(ChunkComm, const_cast<char *>(HOPath.c_str()), MPI_MODE_RDONLY,
      MPIInfo, &HOFile);
    core::EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      Error = error::FILE_OPEN;
      core::LogError(Logger, true, "Unable to open file '%s'.", HOPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }
    auto CloseHO = core::OnScopeExit([&]() {
      core::StartProfileSync(Profiler, MPIIOCloseTime, ChunkComm);
      MPI_File_close(&HOFile);
      core::EndProfile(Profiler, MPIIOCloseTime);
    });

    core::StartProfileSync(Profiler, MPIIOOpenTime, ChunkComm);
    // MPI_File_open missing const qualifier for path string on some platforms
    MPIError = MPI_File_open(ChunkComm, const_cast<char *>(XPath.c_str()), MPI_MODE_RDONLY,
      MPIInfo, &XFile);
    core::EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      Error = error::FILE_OPEN;
      core::LogError(Logger, true, "Unable to open file '%s'.", XPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }
    auto CloseX = core::OnScopeExit([&]() {
      core::StartProfileSync(Profiler, MPIIOCloseTime, ChunkComm);
      MPI_File_close(&XFile);
      core::EndProfile(Profiler, MPIIOCloseTime);
    });

    array<int,2> Sizes({{MAX_DIMS,NumLocalDonors}});

//     MPI_Datatype XDonorSizeType;
//     MPI_Type_vector(MAX_DIMS, NumLocalDonors, NumDonors, MPI_INT, &XDonorSizeType);
//     MPI_Type_commit(&XDonorSizeType);
//     MPI_File_set_view(XFile, XSizesOffset+LocalBegin*sizeof(int), MPI_INT, XDonorSizeType, "native",
//       MPIInfo);
//     int DataSize = MAX_DIMS*NumLocalDonors;
//     File_read_all_endian(XFile, Sizes.Data(), DataSize, MPI_INT, Endian, &Status);
//     MPI_Type_free(&XDonorSizeType);
//     MPI_Get_count(&Status, MPI_INT, &ReadSize);
//     if (ReadSize < DataSize) Error = error::FILE_READ;
//     OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
//     if (Error != error::NONE) {
//       core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i donor sizes from "
//         "XINTOUT file '%s'.", GridID, XPath);
//       OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
//     }

    // The above version causes "Conditional jump or move depends on uninitialized value(s)"
    // errors in valgrind with MPICH 3.2... not sure why
    DatasetOffset = XSizesOffset;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + LocalBegin*sizeof(int);
      core::StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(XFile, ReadOffset, MPI_INT, MPI_INT, "native", MPIInfo);
      core::EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(XFile, Sizes.Data(iDim,0), int(NumLocalDonors), MPI_INT, Endian,
        &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < NumLocalDonors) Error = error::FILE_READ;
      OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
      if (Error != error::NONE) {
        core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i donor sizes from "
          "XINTOUT file '%s'.", GridID, XPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      DatasetOffset += NumDonors*sizeof(int);
    }

    int MaxSize = 0;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (long long iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        MaxSize = std::max(MaxSize, Sizes(iDim,iDonor));
      }
    }

    CreateDonorData(Chunk.Data, NumLocalDonors, MaxSize);
    donor_data &Data = Chunk.Data;

//     MPI_Datatype HODonorCellType;
//     MPI_Type_vector(MAX_DIMS, NumLocalDonors, NumDonors, MPI_INT, &HODonorCellType);
//     MPI_Type_commit(&HODonorCellType);
//     MPI_File_set_view(HOFile, HOCellsOffset+LocalBegin*sizeof(int), MPI_INT, HODonorCellType,
//       "native", MPIInfo);
//     DataSize = MAX_DIMS*NumLocalDonors;
//     File_read_all_endian(HOFile, Data.Extents.Data(), DataSize, MPI_INT, Endian, &Status);
//     MPI_Type_free(&HODonorCellType);
//     MPI_Get_count(&Status, MPI_INT, &ReadSize);
//     if (ReadSize < DataSize) Error = error::FILE_READ;
//     OVK_EH_SYNC(ErrorHandler, Error, ChunkComm); 
//     if (Error != error::NONE) {
//       core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i donor cells from "
//         XINTOUT file '%s'.", GridID, HOPath);
//       OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
//     }

    // The above version causes "Conditional jump or move depends on uninitialized value(s)"
    // errors in valgrind with MPICH 3.2... not sure why
    DatasetOffset = HOCellsOffset;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + LocalBegin*sizeof(int);
      core::StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(HOFile, ReadOffset, MPI_INT, MPI_INT, "native", MPIInfo);
      core::EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(HOFile, Data.Extents.Data(0,iDim,0), int(NumLocalDonors), MPI_INT,
        Endian, &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < NumLocalDonors) Error = error::FILE_READ;
      OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
      if (Error != error::NONE) {
        core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i donor cells from "
          "XINTOUT file '%s'.", GridID, HOPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      DatasetOffset += NumDonors*sizeof(int);
    }

    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (long long iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        // Convert to zero-based indexing
        --Data.Extents(0,iDim,iDonor);
        Data.Extents(1,iDim,iDonor) = Data.Extents(0,iDim,iDonor) + Sizes(iDim,iDonor);
      }
    }

//     MPI_Datatype HODonorCoordsType;
//     MPI_Type_vector(MAX_DIMS, NumLocalDonors, NumDonors, MPI_DOUBLE, &HODonorCoordsType);
//     MPI_Type_commit(&HODonorCoordsType);
//     MPI_File_set_view(HOFile, HOCoordsOffset+LocalBegin*sizeof(double), MPI_DOUBLE,
//       HODonorCoordsType, "native", MPIInfo);
//     DataSize = MAX_DIMS*NumLocalDonors;
//     File_read_all_endian(HOFile, Data.Coords.Data(), DataSize, MPI_DOUBLE, Endian, &Status);
//     MPI_Type_free(&HODonorCoordsType);
//     MPI_Get_count(&Status, MPI_DOUBLE, &ReadSize);
//     if (ReadSize < DataSize) Error = error::FILE_READ;
//     OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
//     if (Error != error::NONE) {
//       core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i donor Coords from "
//         "XINTOUT file '%s'.", GridID, HOPath);
//       OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
//     }

    // The above version causes "Conditional jump or move depends on uninitialized value(s)"
    // errors in valgrind with MPICH 3.2... not sure why
    DatasetOffset = HOCoordsOffset;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + LocalBegin*sizeof(double);
      core::StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(HOFile, ReadOffset, MPI_DOUBLE, MPI_DOUBLE, "native", MPIInfo);
      core::EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(HOFile, Data.Coords.Data(iDim,0), int(NumLocalDonors), MPI_DOUBLE,
        Endian, &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_DOUBLE, &ReadSize);
      if (ReadSize < NumLocalDonors) Error = error::FILE_READ;
      OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
      if (Error != error::NONE) {
        core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i donor Coords from "
          "XINTOUT file '%s'.", GridID, HOPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      DatasetOffset += NumDonors*sizeof(double);
    }

    elem<long long,MAX_DIMS> NumLocalInterpCoefs = {0,0,0};
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (long long iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        int Size = Data.Extents(1,iDim,iDonor) - Data.Extents(0,iDim,iDonor);
        NumLocalInterpCoefs[iDim] += Size;
      }
    }

    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      if (NumLocalInterpCoefs[iDim]*sizeof(double) > std::numeric_limits<int>::max()) {
        Error = error::FILE_READ;
        break;
      }
    }
    OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
    if (Error != error::NONE) {
      core::LogError(Logger, ChunkComm.Rank() == 0, "Donor chunk size too big; increase number of "
        "processes or read granularity.");
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }

    elem<long long,MAX_DIMS> NumInterpCoefs;
    elem<long long,MAX_DIMS> NumInterpCoefsBeforeChunk;
    MPI_Allreduce(&NumLocalInterpCoefs, &NumInterpCoefs, MAX_DIMS, MPI_LONG_LONG, MPI_SUM,
      ChunkComm);
    MPI_Scan(&NumLocalInterpCoefs, &NumInterpCoefsBeforeChunk, MAX_DIMS, MPI_LONG_LONG, MPI_SUM,
      ChunkComm);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      NumInterpCoefsBeforeChunk[iDim] -= NumLocalInterpCoefs[iDim];
    }

    std::vector<double> InterpCoefs(NumLocalInterpCoefs[0]+NumLocalInterpCoefs[1]+
      NumLocalInterpCoefs[2]);
    double *Buffer = InterpCoefs.data();
    DatasetOffset = XInterpCoefsOffset;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + NumInterpCoefsBeforeChunk[iDim]*sizeof(double);
      core::StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(XFile, ReadOffset, MPI_DOUBLE, MPI_DOUBLE, "native", MPIInfo);
      core::EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(XFile, Buffer, int(NumLocalInterpCoefs[iDim]), MPI_DOUBLE, Endian,
        &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_DOUBLE, &ReadSize);
      if (ReadSize < NumLocalInterpCoefs[iDim]) Error = error::FILE_READ;
      OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
      if (Error != error::NONE) {
        core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i donor interpolation "
          "coefficients from XINTOUT file '%s'.", GridID, XPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      Buffer += NumLocalInterpCoefs[iDim];
      DatasetOffset += NumInterpCoefs[iDim]*sizeof(double);
    }

    elem<long long,MAX_DIMS> iNextCoef;
    iNextCoef[0] = 0;
    iNextCoef[1] = iNextCoef[0] + NumLocalInterpCoefs[0];
    iNextCoef[2] = iNextCoef[1] + NumLocalInterpCoefs[1];
    for (long long iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        int Size = Data.Extents(1,iDim,iDonor) - Data.Extents(0,iDim,iDonor);
        for (int iPoint = 0; iPoint < Size; ++iPoint) {
          Data.InterpCoefs(iDim,iPoint,iDonor) = InterpCoefs[iNextCoef[iDim]];
          ++iNextCoef[iDim];
        }
      }
    }

  }

  done_reading:
    OVK_EH_SYNC(ErrorHandler, Error, Comm);
    OVK_EH_CHECK(ErrorHandler, Error);

  return error::NONE;

}

error ReadReceivers(xintout_grid &XINTOUTGrid, const std::string &HOPath,
  long long NumReceivers, MPI_Offset HOPointsOffset, MPI_Offset HOConnectionIDsOffset,
  endian Endian, xintout_format Format, int ReadGranularityAdjust, MPI_Info MPIInfo,
  core::profiler &Profiler) {

  int GridID = XINTOUTGrid.ID;

  const core::comm &Comm = XINTOUTGrid.Comm;

  core::logger &Logger = *XINTOUTGrid.Logger;
  core::error_handler &ErrorHandler = *XINTOUTGrid.ErrorHandler;

  xintout_receivers &XINTOUTReceivers = XINTOUTGrid.Receivers;

  XINTOUTReceivers.Count = NumReceivers;

  // Read chunk on every Nth rank, where N is a power of 2 (or Comm.Size()) chosen such that the
  // number of receivers per chunk is roughly the average number of grid points per rank (subject to
  // user adjustment)
  long long NumPoints =
    (long long)(XINTOUTGrid.GlobalSize[0])*
    (long long)(XINTOUTGrid.GlobalSize[1])*
    (long long)(XINTOUTGrid.GlobalSize[2]);
  long long AvgPointsPerRank = BinDivide(NumPoints, Comm.Size());
  int ChunkRankInterval, NumChunks;
  long long ChunkSize;
  Chunkify(NumReceivers, Comm.Size(), AvgPointsPerRank, ReadGranularityAdjust, ChunkRankInterval,
    NumChunks, ChunkSize);
  bool HasChunk = Comm.Rank() % ChunkRankInterval == 0;

  if (LoggingStatus(Logger)) {
    std::string NumReceiversString = core::FormatNumber(NumReceivers, "receivers", "receiver");
    std::string NumRanksString = core::FormatNumber(NumChunks, "I/O ranks", "I/O rank");
    core::LogStatus(Logger, Comm.Rank() == 0, 1, "Grid %s has %s; using %s.",
      XINTOUTGrid.Name, NumReceiversString, NumRanksString);
  }

  XINTOUTReceivers.ChunkSize = ChunkSize;
  XINTOUTReceivers.HasChunk = HasChunk;

  core::comm ChunkComm = core::CreateSubsetComm(Comm, HasChunk);

  error Error = error::NONE;

  if (HasChunk) {

    xintout_receiver_chunk &Chunk = XINTOUTReceivers.Chunk;

    long long LocalBegin = ChunkSize*ChunkComm.Rank();
    long long LocalEnd = std::min(ChunkSize*(ChunkComm.Rank()+1), NumReceivers);
    long long NumLocalReceivers = LocalEnd - LocalBegin;

    if (NumLocalReceivers*sizeof(long long) > std::numeric_limits<int>::max()) {
      Error = error::FILE_READ;
    }
    OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
    if (Error != error::NONE) {
      core::LogError(Logger, ChunkComm.Rank() == 0, "Receiver chunk size too big; increase number of "
        "processes or read granularity.");
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }

    Chunk.Begin = LocalBegin;
    Chunk.End = LocalEnd;

    CreateReceiverData(Chunk.Data, NumLocalReceivers);
    receiver_data &Data = Chunk.Data;

    Chunk.ConnectionIDs.Resize({NumLocalReceivers});

    MPI_File HOFile;
    MPI_Status Status;
    int MPIError;
    MPI_Offset DatasetOffset, ReadOffset;
    int ReadSize;

    int MPIIOOpenTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Open");
    int MPIIOCloseTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Close");
    int MPIIOOtherTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Other");

    core::StartProfileSync(Profiler, MPIIOOpenTime, ChunkComm);
    // MPI_File_open missing const qualifier for path string on some platforms
    MPIError = MPI_File_open(ChunkComm, const_cast<char *>(HOPath.c_str()), MPI_MODE_RDONLY,
      MPIInfo, &HOFile);
    core::EndProfile(Profiler, MPIIOOpenTime);
    if (MPIError != MPI_SUCCESS) {
      Error = error::FILE_OPEN;
      core::LogError(Logger, true, "Unable to open file '%s'.", HOPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }
    auto CloseHO = core::OnScopeExit([&]() {
      core::StartProfileSync(Profiler, MPIIOCloseTime, ChunkComm);
      MPI_File_close(&HOFile);
      core::EndProfile(Profiler, MPIIOCloseTime);
    });

//     MPI_Datatype HOReceiverPointType;
//     MPI_Type_vector(MAX_DIMS, NumLocalReceivers, NumReceivers, MPI_INT, &HOReceiverPointType);
//     MPI_Type_commit(&HOReceiverPointType);
//     MPI_File_set_view(HOFile, HOPointsOffset+LocalBegin*sizeof(int), MPI_INT, HOReceiverPointType,
//       "native", MPIInfo);
//     DataSize = MAX_DIMS*NumLocalReceivers;
//     File_read_all_endian(HOFile, Data.Points.Data(), DataSize, MPI_INT, Endian, &Status);
//     MPI_Type_free(&HOReceiverPointType);
//     MPI_Get_count(&Status, MPI_INT, &ReadSize);
//     if (ReadSize < DataSize) Error = error::FILE_READ;
//     OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
//     if (Error != error::NONE) {
//       core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i receiver points from "
//         "XINTOUT file '%s'.", GridID, HOPath);
//       OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
//     }

    // The above version causes "Conditional jump or move depends on uninitialized value(s)"
    // errors in valgrind with MPICH 3.2... not sure why
    DatasetOffset = HOPointsOffset;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReadOffset = DatasetOffset + LocalBegin*sizeof(int);
      core::StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
      MPI_File_set_view(HOFile, ReadOffset, MPI_INT, MPI_INT, "native", MPIInfo);
      core::EndProfile(Profiler, MPIIOOtherTime);
      File_read_all_endian(HOFile, Data.Points.Data(iDim,0), int(NumLocalReceivers), MPI_INT,
        Endian, &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < NumLocalReceivers) Error = error::FILE_READ;
      OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
      if (Error != error::NONE) {
        core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i receiver points from "
          "XINTOUT file '%s'.", GridID, HOPath);
        OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
      }
      DatasetOffset += NumReceivers*sizeof(int);
    }

    // Convert to zero-based indexing
    for (long long iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        --Data.Points(iDim,iReceiver);
      }
    }

    int ConnectionIDSize;
    MPI_Datatype ConnectionIDType;
    if (Format == xintout_format::STANDARD) {
      ConnectionIDSize = sizeof(int);
      ConnectionIDType = MPI_INT;
    } else {
      ConnectionIDSize = sizeof(long long);
      ConnectionIDType = MPI_LONG_LONG;
    }

    std::vector<unsigned char> ConnectionIDsBytes(NumLocalReceivers*ConnectionIDSize);

    DatasetOffset = HOConnectionIDsOffset;
    ReadOffset = DatasetOffset+LocalBegin*ConnectionIDSize;
    core::StartProfileSync(Profiler, MPIIOOtherTime, ChunkComm);
    MPI_File_set_view(HOFile, ReadOffset, ConnectionIDType, ConnectionIDType, "native", MPIInfo);
    core::EndProfile(Profiler, MPIIOOtherTime);
    File_read_all_endian(HOFile, ConnectionIDsBytes.data(), int(NumLocalReceivers),
      ConnectionIDType, Endian, &Status, Profiler, ChunkComm);
    MPI_Get_count(&Status, ConnectionIDType, &ReadSize);
    if (ReadSize < NumLocalReceivers) Error = error::FILE_READ;
    OVK_EH_SYNC(ErrorHandler, Error, ChunkComm);
    if (Error != error::NONE) {
      core::LogError(Logger, ChunkComm.Rank() == 0, "Unable to read grid %i receiver connection IDs "
        "from XINTOUT file '%s'.", GridID, HOPath);
      OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
    }

    if (Format == xintout_format::STANDARD) {
      auto ConnectionIDs = reinterpret_cast<int *>(ConnectionIDsBytes.data());
      for (long long iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        // Convert to zero-based indexing
        Chunk.ConnectionIDs(iReceiver) = ConnectionIDs[iReceiver]-1;
      }
    } else {
      auto ConnectionIDs = reinterpret_cast<long long *>(ConnectionIDsBytes.data());
      for (long long iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        // Convert to zero-based indexing
        Chunk.ConnectionIDs(iReceiver) = ConnectionIDs[iReceiver]-1;
      }
    }

  }

  done_reading:
    OVK_EH_SYNC(ErrorHandler, Error, Comm);
    OVK_EH_CHECK(ErrorHandler, Error);

  return error::NONE;

}

void MatchDonorsAndReceivers(xintout &XINTOUT, core::profiler &Profiler) {

  const core::comm &Comm = XINTOUT.Comm;

  MPI_Barrier(Comm);

  core::logger &Logger = *XINTOUT.Logger;

  core::LogStatus(Logger, Comm.Rank() == 0, 0, "Matching donors and receivers...");

  int MapToBinsTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Match::MapToBins");
  int HandshakeTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Match::Handshake");
  int SendToBinsTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Match::SendToBins");
  int FillConnectionDataTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Match::FillConnectionData");
  int RecvFromBinsTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Match::RecvFromBins");
  int UnpackTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Match::Unpack");

  int NumGrids = XINTOUT.NumGrids;
  int NumLocalGrids = XINTOUT.NumLocalGrids;

  std::vector<long long> NumPointsOnGrid(NumGrids, 0);
  std::vector<long long> NumReceiversOnGrid(NumGrids, 0);

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    xintout_grid &XINTOUTGrid = XINTOUT.Grids[iLocalGrid];
    int GridID = XINTOUTGrid.ID;
    int iGrid = GridID-1;
    NumPointsOnGrid[iGrid] = 1;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      NumPointsOnGrid[iGrid] *= XINTOUTGrid.GlobalSize[iDim];
    }
    NumReceiversOnGrid[iGrid] = XINTOUTGrid.Receivers.Count;
  }

  MPI_Allreduce(MPI_IN_PLACE, NumPointsOnGrid.data(), NumGrids, MPI_LONG_LONG, MPI_MAX, Comm);
  MPI_Allreduce(MPI_IN_PLACE, NumReceiversOnGrid.data(), NumGrids, MPI_LONG_LONG, MPI_MAX, Comm);

  long long NumPoints = 0;
  long long NumConnections = 0;
  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    NumPoints += NumPointsOnGrid[iGrid];
    NumConnections += NumReceiversOnGrid[iGrid];
  }

  NumPointsOnGrid.clear();
  NumReceiversOnGrid.clear();

  xintout_connections &XINTOUTConnections = XINTOUT.Connections;

  long long BinSize = std::min(BinDivide(NumPoints, Comm.Size()), NumConnections);
  bool HasBin = BinSize*Comm.Rank() < NumConnections;

  XINTOUTConnections.BinSize = BinSize;
  XINTOUTConnections.HasBin = HasBin;

  xintout_connection_bin &Bin = XINTOUTConnections.Bin;
  connection_data &BinData = Bin.Data;

  if (HasBin) {
    Bin.Begin = BinSize*Comm.Rank();
    Bin.End = std::min(BinSize*(Comm.Rank()+1), NumConnections);
    long long NumLocalConnections = Bin.End - Bin.Begin;
    CreateConnectionData(BinData, NumLocalConnections);
  }

  struct send_recv {
    int Count;
    connection_data Data;
    send_recv():
      Count(0)
    {}
  };

  core::StartProfileSync(Profiler, MapToBinsTime, Comm);

  std::map<int, send_recv> DonorSends, ReceiverSends;

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    xintout_grid &XINTOUTGrid = XINTOUT.Grids[iLocalGrid];
    xintout_donors &XINTOUTDonors = XINTOUTGrid.Donors;
    xintout_receivers &XINTOUTReceivers = XINTOUTGrid.Receivers;
    if (XINTOUTDonors.HasChunk) {
      xintout_donor_chunk &DonorChunk = XINTOUTDonors.Chunk;
      long long NumLocalDonors = DonorChunk.End - DonorChunk.Begin;
      for (long long iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        long long ConnectionID = DonorChunk.StartingConnectionID + iDonor;
        int ConnectionBinIndex = int(ConnectionID/BinSize);
        auto Iter = DonorSends.lower_bound(ConnectionBinIndex);
        if (Iter == DonorSends.end() || Iter->first > ConnectionBinIndex) {
          Iter = DonorSends.emplace_hint(Iter, ConnectionBinIndex, send_recv());
        }
        send_recv &Send = Iter->second;
        ++Send.Count;
      }
    }
    if (XINTOUTReceivers.HasChunk) {
      xintout_receiver_chunk &ReceiverChunk = XINTOUTReceivers.Chunk;
      long long NumLocalReceivers = ReceiverChunk.End - ReceiverChunk.Begin;
      for (long long iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        long long ConnectionID = ReceiverChunk.ConnectionIDs(iReceiver);
        int ConnectionBinIndex = int(ConnectionID/BinSize);
        auto Iter = ReceiverSends.find(ConnectionBinIndex);
        if (Iter == ReceiverSends.end() || Iter->first > ConnectionBinIndex) {
          Iter = ReceiverSends.emplace_hint(Iter, ConnectionBinIndex, send_recv());
        }
        send_recv &Send = Iter->second;
        ++Send.Count;
      }
    }
  }

  for (auto &Pair : DonorSends) {
    send_recv &Send = Pair.second;
    CreateConnectionData(Send.Data, Send.Count);
    // Reset count for filling in data
    Send.Count = 0;
  }

  for (auto &Pair : ReceiverSends) {
    send_recv &Send = Pair.second;
    CreateConnectionData(Send.Data, Send.Count);
    // Reset count for filling in data
    Send.Count = 0;
  }

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    xintout_grid &XINTOUTGrid = XINTOUT.Grids[iLocalGrid];
    xintout_donors &XINTOUTDonors = XINTOUTGrid.Donors;
    xintout_receivers &XINTOUTReceivers = XINTOUTGrid.Receivers;
    if (XINTOUTDonors.HasChunk) {
      xintout_donor_chunk &DonorChunk = XINTOUTDonors.Chunk;
      long long NumLocalDonors = DonorChunk.End - DonorChunk.Begin;
      for (long long iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        long long ConnectionID = DonorChunk.StartingConnectionID + iDonor;
        int ConnectionBinIndex = int(ConnectionID/BinSize);
        send_recv &Send = DonorSends[ConnectionBinIndex];
        long long iNext = Send.Count;
        Send.Data.ConnectionIDs(iNext) = ConnectionID;
        Send.Data.DonorGridIDs(iNext) = XINTOUTGrid.ID;
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send.Data.DonorCells(iDim,iNext) = DonorChunk.Data.Extents(0,iDim,iDonor);
        }
        ++Send.Count;
      }
    }
    if (XINTOUTReceivers.HasChunk) {
      xintout_receiver_chunk &ReceiverChunk = XINTOUTReceivers.Chunk;
      long long NumLocalReceivers = ReceiverChunk.End - ReceiverChunk.Begin;
      for (long long iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        long long ConnectionID = ReceiverChunk.ConnectionIDs(iReceiver);
        int ConnectionBinIndex = int(ConnectionID/BinSize);
        send_recv &Send = ReceiverSends[ConnectionBinIndex];
        long long iNext = Send.Count;
        Send.Data.ConnectionIDs(iNext) = ConnectionID;
        Send.Data.ReceiverGridIDs(iNext) = XINTOUTGrid.ID;
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send.Data.ReceiverPoints(iDim,iNext) = ReceiverChunk.Data.Points(iDim,iReceiver);
        }
        ++Send.Count;
      }
    }
  }

  int NumDonorSends = DonorSends.size();
  int NumReceiverSends = ReceiverSends.size();

  std::vector<int> DonorSendToRanks;
  for (auto &Pair : DonorSends) {
    int Rank = Pair.first;
    DonorSendToRanks.push_back(Rank);
  }

  std::vector<int> ReceiverSendToRanks;
  for (auto &Pair : ReceiverSends) {
    int Rank = Pair.first;
    ReceiverSendToRanks.push_back(Rank);
  }

  core::EndProfile(Profiler, MapToBinsTime);
  core::StartProfileSync(Profiler, HandshakeTime, Comm);

  array<int> DonorRecvFromRanks, ReceiverRecvFromRanks;

  core::DynamicHandshake(Comm, DonorSendToRanks, DonorRecvFromRanks);
  core::DynamicHandshake(Comm, ReceiverSendToRanks, ReceiverRecvFromRanks);

  DonorSendToRanks.clear();
  ReceiverSendToRanks.clear();

  core::EndProfile(Profiler, HandshakeTime);
  core::StartProfileSync(Profiler, SendToBinsTime, Comm);

  std::map<int, send_recv> DonorRecvs;
  for (int Rank : DonorRecvFromRanks) {
    DonorRecvs.emplace(Rank, send_recv());
  }

  std::map<int, send_recv> ReceiverRecvs;
  for (int Rank : ReceiverRecvFromRanks) {
    ReceiverRecvs.emplace(Rank, send_recv());
  }

  DonorRecvFromRanks.Clear();
  ReceiverRecvFromRanks.Clear();

  int NumDonorRecvs = DonorRecvs.size();
  int NumReceiverRecvs = ReceiverRecvs.size();

  std::vector<MPI_Request> Requests;

  // 3 requests per send/recv (1 for each of grid IDs, connection IDs, and points/cells)
  Requests.reserve(3*(NumDonorSends+NumDonorRecvs+NumReceiverSends+NumReceiverRecvs));

  auto Isend = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int DestRank, int Tag,
    MPI_Comm Comm) {
    Requests.emplace_back();
    MPI_Isend(Buffer, Count, DataType, DestRank, Tag, Comm, &Requests.back());

  };
  auto Irecv = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int SourceRank, int Tag,
    MPI_Comm Comm) {
    Requests.emplace_back();
    MPI_Irecv(Buffer, Count, DataType, SourceRank, Tag, Comm, &Requests.back());
  };

  for (auto &Pair : DonorRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    Irecv(&Recv.Count, 1, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    Irecv(&Recv.Count, 1, MPI_LONG_LONG, Rank, 1, Comm);
  }

  for (auto &Pair : DonorSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Isend(&Send.Count, 1, MPI_LONG_LONG, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Isend(&Send.Count, 1, MPI_LONG_LONG, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.size(), Requests.data(), MPI_STATUSES_IGNORE);

  Requests.clear();

  for (auto &Pair : DonorRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    CreateConnectionData(Recv.Data, Recv.Count);
    Irecv(Recv.Data.ConnectionIDs.Data(), Recv.Count, MPI_LONG_LONG, Rank, 0, Comm);
    Irecv(Recv.Data.DonorGridIDs.Data(), Recv.Count, MPI_INT, Rank, 0, Comm);
    Irecv(Recv.Data.DonorCells.Data(), MAX_DIMS*Recv.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    CreateConnectionData(Recv.Data, Recv.Count);
    Irecv(Recv.Data.ConnectionIDs.Data(), Recv.Count, MPI_LONG_LONG, Rank, 1, Comm);
    Irecv(Recv.Data.ReceiverGridIDs.Data(), Recv.Count, MPI_INT, Rank, 1, Comm);
    Irecv(Recv.Data.ReceiverPoints.Data(), MAX_DIMS*Recv.Count, MPI_INT, Rank, 1, Comm);
  }

  for (auto &Pair : DonorSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Isend(Send.Data.ConnectionIDs.Data(), Send.Count, MPI_LONG_LONG, Rank, 0, Comm);
    Isend(Send.Data.DonorGridIDs.Data(), Send.Count, MPI_INT, Rank, 0, Comm);
    Isend(Send.Data.DonorCells.Data(), MAX_DIMS*Send.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Isend(Send.Data.ConnectionIDs.Data(), Send.Count, MPI_LONG_LONG, Rank, 1, Comm);
    Isend(Send.Data.ReceiverGridIDs.Data(), Send.Count, MPI_INT, Rank, 1, Comm);
    Isend(Send.Data.ReceiverPoints.Data(), MAX_DIMS*Send.Count, MPI_INT, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.size(), Requests.data(), MPI_STATUSES_IGNORE);

  Requests.clear();

  core::EndProfile(Profiler, SendToBinsTime);
  core::StartProfileSync(Profiler, RecvFromBinsTime, Comm);

  for (auto &Pair : DonorSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Irecv(Send.Data.ReceiverGridIDs.Data(), Send.Count, MPI_INT, Rank, 0, Comm);
    Irecv(Send.Data.ReceiverPoints.Data(), MAX_DIMS*Send.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverSends) {
    int Rank = Pair.first;
    send_recv &Send = Pair.second;
    Irecv(Send.Data.DonorGridIDs.Data(), Send.Count, MPI_INT, Rank, 1, Comm);
    Irecv(Send.Data.DonorCells.Data(), MAX_DIMS*Send.Count, MPI_INT, Rank, 1, Comm);
  }

  core::EndProfile(Profiler, RecvFromBinsTime);
  core::StartProfileSync(Profiler, FillConnectionDataTime, Comm);

  for (auto &Pair : DonorRecvs) {
    send_recv &Recv = Pair.second;
    for (long long iDonor = 0; iDonor < Recv.Count; ++iDonor) {
      long long iConnection = Recv.Data.ConnectionIDs(iDonor) - Bin.Begin;
      BinData.DonorGridIDs(iConnection) = Recv.Data.DonorGridIDs(iDonor);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        BinData.DonorCells(iDim,iConnection) = Recv.Data.DonorCells(iDim,iDonor);
      }
    }
  }

  for (auto &Pair : ReceiverRecvs) {
    send_recv &Recv = Pair.second;
    for (long long iReceiver = 0; iReceiver < Recv.Count; ++iReceiver) {
      long long iConnection = Recv.Data.ConnectionIDs(iReceiver) - Bin.Begin;
      BinData.ReceiverGridIDs(iConnection) = Recv.Data.ReceiverGridIDs(iReceiver);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        BinData.ReceiverPoints(iDim,iConnection) = Recv.Data.ReceiverPoints(iDim,iReceiver);
      }
    }
  }

  core::EndProfile(Profiler, FillConnectionDataTime);
  core::StartProfileSync(Profiler, RecvFromBinsTime, Comm);

  for (auto &Pair : DonorRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    for (long long iDonor = 0; iDonor < Recv.Count; ++iDonor) {
      long long iConnection = Recv.Data.ConnectionIDs(iDonor) - Bin.Begin;
      Recv.Data.ReceiverGridIDs(iDonor) = BinData.ReceiverGridIDs(iConnection);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Recv.Data.ReceiverPoints(iDim,iDonor) = BinData.ReceiverPoints(iDim,iConnection);
      }
    }
    Isend(Recv.Data.ReceiverGridIDs.Data(), Recv.Count, MPI_INT, Rank, 0, Comm);
    Isend(Recv.Data.ReceiverPoints.Data(), MAX_DIMS*Recv.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverRecvs) {
    int Rank = Pair.first;
    send_recv &Recv = Pair.second;
    for (long long iReceiver = 0; iReceiver < Recv.Count; ++iReceiver) {
      long long iConnection = Recv.Data.ConnectionIDs(iReceiver) - Bin.Begin;
      Recv.Data.DonorGridIDs(iReceiver) = BinData.DonorGridIDs(iConnection);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Recv.Data.DonorCells(iDim,iReceiver) = BinData.DonorCells(iDim,iConnection);
      }
    }
    Isend(Recv.Data.DonorGridIDs.Data(), Recv.Count, MPI_INT, Rank, 1, Comm);
    Isend(Recv.Data.DonorCells.Data(), MAX_DIMS*Recv.Count, MPI_INT, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.size(), Requests.data(), MPI_STATUSES_IGNORE);

  Requests.clear();

  core::EndProfile(Profiler, RecvFromBinsTime);
  core::StartProfileSync(Profiler, UnpackTime, Comm);

  for (auto &Pair : DonorSends) {
    send_recv &Send = Pair.second;
    // Reset count to 0 for unpacking
    Send.Count = 0;
  }

  for (auto &Pair : ReceiverSends) {
    send_recv &Send = Pair.second;
    // Reset count to 0 for unpacking
    Send.Count = 0;
  }

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    xintout_grid &XINTOUTGrid = XINTOUT.Grids[iLocalGrid];
    xintout_donors &XINTOUTDonors = XINTOUTGrid.Donors;
    xintout_receivers &XINTOUTReceivers = XINTOUTGrid.Receivers;
    if (XINTOUTDonors.HasChunk) {
      xintout_donor_chunk &DonorChunk = XINTOUTDonors.Chunk;
      long long NumLocalDonors = DonorChunk.End - DonorChunk.Begin;
      for (long long iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        long long ConnectionID = DonorChunk.StartingConnectionID + iDonor;
        int ConnectionBinIndex = int(ConnectionID/BinSize);
        send_recv &Send = DonorSends[ConnectionBinIndex];
        long long iNext = Send.Count;
        DonorChunk.Data.DestinationGridIDs(iDonor) = Send.Data.ReceiverGridIDs(iNext);
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          DonorChunk.Data.DestinationPoints(iDim,iDonor) = Send.Data.ReceiverPoints(iDim,iNext);
        }
        ++Send.Count;
      }
    }
    if (XINTOUTReceivers.HasChunk) {
      xintout_receiver_chunk &ReceiverChunk = XINTOUTReceivers.Chunk;
      long long NumLocalReceivers = ReceiverChunk.End - ReceiverChunk.Begin;
      for (long long iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
        long long ConnectionID = ReceiverChunk.ConnectionIDs(iReceiver);
        int ConnectionBinIndex = int(ConnectionID/BinSize);
        send_recv &Send = ReceiverSends[ConnectionBinIndex];
        long long iNext = Send.Count;
        ReceiverChunk.Data.SourceGridIDs(iReceiver) = Send.Data.DonorGridIDs(iNext);
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          ReceiverChunk.Data.SourceCells(iDim,iReceiver) = Send.Data.DonorCells(iDim,iNext);
        }
        ++Send.Count;
      }
    }
  }

  core::EndProfile(Profiler, UnpackTime);

  MPI_Barrier(Comm);

  core::LogStatus(Logger, Comm.Rank() == 0, 0, "Finished matching donors and receivers.");

}

void DistributeConnectivityData(const xintout &XINTOUT, const std::vector<const grid *> &LocalGrids,
  std::vector<donor_data> &LocalDonorData, std::vector<receiver_data> &LocalReceiverData,
  core::profiler &Profiler) {

  const core::comm &Comm = XINTOUT.Comm;

  MPI_Barrier(Comm);

  core::logger &Logger = *XINTOUT.Logger;

  core::LogStatus(Logger, Comm.Rank() == 0, 0, "Distributing connectivity data to ranks...");

  for (int iLocalGrid = 0; iLocalGrid < XINTOUT.NumLocalGrids; ++iLocalGrid) {
    const xintout_grid &XINTOUTGrid = XINTOUT.Grids[iLocalGrid];
    const grid &Grid = *LocalGrids[iLocalGrid];
    DistributeGridConnectivityData(XINTOUTGrid, Grid, LocalDonorData[iLocalGrid],
      LocalReceiverData[iLocalGrid], Profiler);
  }

  MPI_Barrier(Comm);

  core::LogStatus(Logger, Comm.Rank() == 0, 0, "Finished distributing connectivity data to ranks.");

}

void DistributeGridConnectivityData(const xintout_grid &XINTOUTGrid, const grid &Grid,
  donor_data &DonorData, receiver_data &ReceiverData, core::profiler &Profiler) {

  int NumDims = XINTOUTGrid.NumDims;
  const core::comm &Comm = XINTOUTGrid.Comm;

  int MapToBinsTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Distribute::MapToBins");
  int RetrieveBinsTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Distribute::RetrieveBins");
  int FindRanksTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Distribute::FindRanks");
  int HandshakeTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Distribute::Handshake");
  int SendDataTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Distribute::SendData");

  const xintout_donors &XINTOUTDonors = XINTOUTGrid.Donors;
  const xintout_receivers &XINTOUTReceivers = XINTOUTGrid.Receivers;

  cart Cart;
  range GlobalRange;
  GetGridCart(Grid, Cart);
  GetGridGlobalRange(Grid, GlobalRange);

  const core::partition_hash &Hash = core::GetGridPartitionHash(Grid);

  core::StartProfileSync(Profiler, MapToBinsTime, Comm);

  std::map<int, core::partition_bin> Bins;

  long long NumChunkDonors = 0;
  long long NumChunkDonorPoints = 0;
  array<int,2> ChunkDonorPoints;
  std::vector<int> ChunkDonorPointBinIndices;
  if (XINTOUTDonors.HasChunk) {
    const xintout_donor_chunk &DonorChunk = XINTOUTDonors.Chunk;
    NumChunkDonors = DonorChunk.End - DonorChunk.Begin;
    for (long long iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      int NumPointsInCell = 1;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        NumPointsInCell *= DonorChunk.Data.Extents(1,iDim,iDonor) -
          DonorChunk.Data.Extents(0,iDim,iDonor);
      }
      NumChunkDonorPoints += NumPointsInCell;
    }
    ChunkDonorPoints.Resize({{MAX_DIMS,NumChunkDonorPoints}});
    long long iDonorPoint = 0;
    for (long long iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      elem<int,MAX_DIMS> DonorBegin, DonorEnd;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorBegin[iDim] = DonorChunk.Data.Extents(0,iDim,iDonor);
        DonorEnd[iDim] = DonorChunk.Data.Extents(1,iDim,iDonor);
      }
      range DonorRange(NumDims, DonorBegin, DonorEnd);
      bool AwayFromEdge = RangeIncludes(GlobalRange, DonorRange);
      if (AwayFromEdge) {
        for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
          for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
            for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
              elem<int,MAX_DIMS> Point = {i,j,k};
              for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
                ChunkDonorPoints(iDim,iDonorPoint) = Point[iDim];
              }
              ++iDonorPoint;
            }
          }
        }
      } else {
        for (int k = DonorBegin[2]; k < DonorEnd[2]; ++k) {
          for (int j = DonorBegin[1]; j < DonorEnd[1]; ++j) {
            for (int i = DonorBegin[0]; i < DonorEnd[0]; ++i) {
              elem<int,MAX_DIMS> Point = {i,j,k};
              CartPeriodicAdjust(Cart, Point, Point);
              for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
                ChunkDonorPoints(iDim,iDonorPoint) = Point[iDim];
              }
              ++iDonorPoint;
            }
          }
        }
      }
    }
    ChunkDonorPointBinIndices.resize(NumChunkDonorPoints);
    core::MapToPartitionBins(Hash, ChunkDonorPoints, ChunkDonorPointBinIndices);
    for (long long iDonorPoint = 0; iDonorPoint < NumChunkDonorPoints; ++iDonorPoint) {
      int BinIndex = ChunkDonorPointBinIndices[iDonorPoint];
      auto Iter = Bins.lower_bound(BinIndex);
      if (Iter == Bins.end() || Iter->first > BinIndex) {
        Bins.emplace_hint(Iter, BinIndex, core::partition_bin());
      }
    }
  }

  long long NumChunkReceivers = 0;
  std::vector<int> ChunkReceiverBinIndices;
  if (XINTOUTReceivers.HasChunk) {
    const xintout_receiver_chunk &ReceiverChunk = XINTOUTReceivers.Chunk;
    NumChunkReceivers = ReceiverChunk.End - ReceiverChunk.Begin;
    ChunkReceiverBinIndices.resize(NumChunkReceivers);
    core::MapToPartitionBins(Hash, ReceiverChunk.Data.Points, ChunkReceiverBinIndices);
    for (long long iReceiver = 0; iReceiver < NumChunkReceivers; ++iReceiver) {
      int BinIndex = ChunkReceiverBinIndices[iReceiver];
      auto Iter = Bins.lower_bound(BinIndex);
      if (Iter == Bins.end() || Iter->first > BinIndex) {
        Bins.emplace_hint(Iter, BinIndex, core::partition_bin());
      }
    }
  }

  core::EndProfile(Profiler, MapToBinsTime);
  core::StartProfileSync(Profiler, RetrieveBinsTime, Comm);

  core::RetrievePartitionBins(Hash, Bins);

  core::EndProfile(Profiler, RetrieveBinsTime);
  core::StartProfileSync(Profiler, FindRanksTime, Comm);

  std::vector<int> NumChunkDonorRanks;
  std::vector<int *> ChunkDonorRanks;
  std::vector<int> ChunkDonorRanksData;
  if (XINTOUTDonors.HasChunk) {
    const xintout_donor_chunk &DonorChunk = XINTOUTDonors.Chunk;
    NumChunkDonorRanks.resize(NumChunkDonors);
    ChunkDonorRanks.resize(NumChunkDonors);
    ChunkDonorRanksData.resize(NumChunkDonorPoints);
    core::FindPartitions(Hash, Bins, ChunkDonorPoints, ChunkDonorPointBinIndices,
      ChunkDonorRanksData);
    int MaxSize = DonorChunk.Data.MaxSize;
    int MaxPointsInCell = 1;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      MaxPointsInCell *= MaxSize;
    }
    std::vector<int> UniqueRanks(MaxPointsInCell);
    long long iDonorPoint = 0;
    for (long long iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      NumChunkDonorRanks[iDonor] = 0;
      ChunkDonorRanks[iDonor] = ChunkDonorRanksData.data() + iDonorPoint;
      int NumPointsInCell = 1;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        NumPointsInCell *= DonorChunk.Data.Extents(1,iDim,iDonor) -
          DonorChunk.Data.Extents(0,iDim,iDonor);
      }
      for (int iPointInCell = 0; iPointInCell < NumPointsInCell; ++iPointInCell) {
        int Rank = ChunkDonorRanks[iDonor][iPointInCell];
        int iRank = 0;
        while (iRank < NumChunkDonorRanks[iDonor] && UniqueRanks[iRank] != Rank) {
          ++iRank;
        }
        if (iRank == NumChunkDonorRanks[iDonor]) {
          UniqueRanks[iRank] = Rank;
          ++NumChunkDonorRanks[iDonor];
        }
      }
      for (int iRank = 0; iRank < NumChunkDonorRanks[iDonor]; ++iRank) {
        ChunkDonorRanks[iDonor][iRank] = UniqueRanks[iRank];
      }
      iDonorPoint += NumPointsInCell;
    }
    ChunkDonorPoints.Clear();
    ChunkDonorPointBinIndices.clear();
  }

  std::vector<int> ChunkReceiverRanks;
  if (XINTOUTReceivers.HasChunk) {
    const xintout_receiver_chunk &ReceiverChunk = XINTOUTReceivers.Chunk;
    ChunkReceiverRanks.resize(NumChunkReceivers);
    core::FindPartitions(Hash, Bins, ReceiverChunk.Data.Points, ChunkReceiverBinIndices,
      ChunkReceiverRanks);
    ChunkReceiverBinIndices.clear();
  }

  Bins.clear();

  struct donor_send_recv {
    int Count;
    int MaxSize;
    donor_data Data;
    donor_send_recv():
      Count(0),
      MaxSize(0)
    {}
  };

  struct receiver_send_recv {
    int Count;
    receiver_data Data;
    receiver_send_recv():
      Count(0)
    {}
  };

  std::map<int, donor_send_recv> DonorSends;
  std::map<int, receiver_send_recv> ReceiverSends;

  if (XINTOUTDonors.HasChunk) {
    const xintout_donor_chunk &DonorChunk = XINTOUTDonors.Chunk;
    for (long long iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      for (int iRank = 0; iRank < NumChunkDonorRanks[iDonor]; ++iRank) {
        int Rank = ChunkDonorRanks[iDonor][iRank];
        auto Iter = DonorSends.lower_bound(Rank);
        if (Iter == DonorSends.end() || Iter->first > Rank) {
          Iter = DonorSends.emplace_hint(Iter, Rank, donor_send_recv());
        }
        donor_send_recv &Send = Iter->second;
        ++Send.Count;
        Send.MaxSize = std::max(Send.MaxSize, DonorChunk.Data.MaxSize);
      }
    }
  }

  if (XINTOUTReceivers.HasChunk) {
    for (long long iReceiver = 0; iReceiver < NumChunkReceivers; ++iReceiver) {
      int Rank = ChunkReceiverRanks[iReceiver];
      auto Iter = ReceiverSends.lower_bound(Rank);
      if (Iter == ReceiverSends.end() || Iter->first > Rank) {
        Iter = ReceiverSends.emplace_hint(Iter, Rank, receiver_send_recv());
      }
      receiver_send_recv &Send = Iter->second;
      ++Send.Count;
    }
  }

  for (auto &Pair : DonorSends) {
    donor_send_recv &Send = Pair.second;
    CreateDonorData(Send.Data, Send.Count, Send.MaxSize);
    // Reset count for filling in data
    Send.Count = 0;
  }

  for (auto &Pair : ReceiverSends) {
    receiver_send_recv &Send = Pair.second;
    CreateReceiverData(Send.Data, Send.Count);
    // Reset count for filling in data
    Send.Count = 0;
  }

  if (XINTOUTDonors.HasChunk) {
    const xintout_donor_chunk &DonorChunk = XINTOUTDonors.Chunk;
    for (long long iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      for (int iRank = 0; iRank < NumChunkDonorRanks[iDonor]; ++iRank) {
        int Rank = ChunkDonorRanks[iDonor][iRank];
        donor_send_recv &Send = DonorSends[Rank];
        long long iNext = Send.Count;
        int MaxSize = Send.MaxSize;
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send.Data.Extents(0,iDim,iNext) = DonorChunk.Data.Extents(0,iDim,iDonor);
          Send.Data.Extents(1,iDim,iNext) = DonorChunk.Data.Extents(1,iDim,iDonor);
        }
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send.Data.Coords(iDim,iNext) = DonorChunk.Data.Coords(iDim,iDonor);
        }
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          int Size = DonorChunk.Data.Extents(1,iDim,iDonor) -
            DonorChunk.Data.Extents(0,iDim,iDonor);
          for (int iPoint = 0; iPoint < Size; ++iPoint) {
            Send.Data.InterpCoefs(iDim,iPoint,iNext) =
              DonorChunk.Data.InterpCoefs(iDim,iPoint,iDonor);
          }
          // Initialize the rest with zeros since we're technically touching all of the data
          // when sending it
          for (int iPoint = Size; iPoint < MaxSize; ++iPoint) {
            Send.Data.InterpCoefs(iDim,iPoint,iNext) = 0.;
          }
        }
        Send.Data.DestinationGridIDs(iNext) = DonorChunk.Data.DestinationGridIDs(iDonor);
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          Send.Data.DestinationPoints(iDim,iNext) = DonorChunk.Data.DestinationPoints(iDim,iDonor);
        }
        ++Send.Count;
      }
    }
    NumChunkDonorRanks.clear();
    ChunkDonorRanks.clear();
    ChunkDonorRanksData.clear();
  }

  if (XINTOUTReceivers.HasChunk) {
    const xintout_receiver_chunk &ReceiverChunk = XINTOUTReceivers.Chunk;
    for (long long iReceiver = 0; iReceiver < NumChunkReceivers; ++iReceiver) {
      int Rank = ChunkReceiverRanks[iReceiver];
      receiver_send_recv &Send = ReceiverSends[Rank];
      long long iNext = Send.Count;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Send.Data.Points(iDim,iNext) = ReceiverChunk.Data.Points(iDim,iReceiver);
      }
      Send.Data.SourceGridIDs(iNext) = ReceiverChunk.Data.SourceGridIDs(iReceiver);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Send.Data.SourceCells(iDim,iNext) = ReceiverChunk.Data.SourceCells(iDim,iReceiver);
      }
      ++Send.Count;
    }
    ChunkReceiverRanks.clear();
  }

  int NumDonorSends = DonorSends.size();
  int NumReceiverSends = ReceiverSends.size();

  std::vector<int> DonorSendToRanks;
  for (auto &Pair : DonorSends) {
    int Rank = Pair.first;
    DonorSendToRanks.push_back(Rank);
  }

  std::vector<int> ReceiverSendToRanks;
  for (auto &Pair : ReceiverSends) {
    int Rank = Pair.first;
    ReceiverSendToRanks.push_back(Rank);
  }

  core::EndProfile(Profiler, FindRanksTime);
  core::StartProfileSync(Profiler, HandshakeTime, Comm);

  array<int> DonorRecvFromRanks, ReceiverRecvFromRanks;

  core::DynamicHandshake(Comm, DonorSendToRanks, DonorRecvFromRanks);
  core::DynamicHandshake(Comm, ReceiverSendToRanks, ReceiverRecvFromRanks);

  DonorSendToRanks.clear();
  ReceiverSendToRanks.clear();

  core::EndProfile(Profiler, HandshakeTime);
  core::StartProfileSync(Profiler, SendDataTime, Comm);

  std::map<int, donor_send_recv> DonorRecvs;
  for (int Rank : DonorRecvFromRanks) {
    DonorRecvs.emplace(Rank, donor_send_recv());
  }

  std::map<int, receiver_send_recv> ReceiverRecvs;
  for (int Rank : ReceiverRecvFromRanks) {
    ReceiverRecvs.emplace(Rank, receiver_send_recv());
  }

  DonorRecvFromRanks.Clear();
  ReceiverRecvFromRanks.Clear();

  int NumDonorRecvs = DonorRecvs.size();
  int NumReceiverRecvs = ReceiverRecvs.size();

  std::vector<MPI_Request> Requests;

  // 5 requests per donor send/recv (1 for each of Extents, Coords, interp coefs, destination grid
  // IDs, and destination points)
  // 3 requests per receiver send/recv (1 for each of points, source grid IDs, and source cells)
  Requests.reserve(5*(NumDonorSends+NumDonorRecvs)+3*(NumReceiverSends+NumReceiverRecvs));

  auto Isend = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int DestRank, int Tag,
    MPI_Comm Comm) {
    Requests.emplace_back();
    MPI_Isend(Buffer, Count, DataType, DestRank, Tag, Comm, &Requests.back());

  };
  auto Irecv = [&Requests](void *Buffer, int Count, MPI_Datatype DataType, int SourceRank, int Tag,
    MPI_Comm Comm) {
    Requests.emplace_back();
    MPI_Irecv(Buffer, Count, DataType, SourceRank, Tag, Comm, &Requests.back());
  };

  for (auto &Pair : DonorRecvs) {
    int Rank = Pair.first;
    donor_send_recv &Recv = Pair.second;
    Irecv(&Recv.Count, 1, MPI_LONG_LONG, Rank, 0, Comm);
    Irecv(&Recv.MaxSize, 1, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverRecvs) {
    int Rank = Pair.first;
    receiver_send_recv &Recv = Pair.second;
    Irecv(&Recv.Count, 1, MPI_LONG_LONG, Rank, 1, Comm);
  }

  for (auto &Pair : DonorSends) {
    int Rank = Pair.first;
    donor_send_recv &Send = Pair.second;
    Isend(&Send.Count, 1, MPI_LONG_LONG, Rank, 0, Comm);
    Isend(&Send.MaxSize, 1, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverSends) {
    int Rank = Pair.first;
    receiver_send_recv &Send = Pair.second;
    Isend(&Send.Count, 1, MPI_LONG_LONG, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.size(), Requests.data(), MPI_STATUSES_IGNORE);

  Requests.clear();

  for (auto &Pair : DonorRecvs) {
    int Rank = Pair.first;
    donor_send_recv &Recv = Pair.second;
    CreateDonorData(Recv.Data, Recv.Count, Recv.MaxSize);
    Irecv(Recv.Data.Extents.Data(), 2*MAX_DIMS*Recv.Count, MPI_INT, Rank, 0, Comm);
    Irecv(Recv.Data.Coords.Data(), MAX_DIMS*Recv.Count, MPI_DOUBLE, Rank, 0, Comm);
    Irecv(Recv.Data.InterpCoefs.Data(), MAX_DIMS*Recv.MaxSize*Recv.Count, MPI_DOUBLE, Rank, 0,
      Comm);
    Irecv(Recv.Data.DestinationGridIDs.Data(), Recv.Count, MPI_INT, Rank, 0, Comm);
    Irecv(Recv.Data.DestinationPoints.Data(), MAX_DIMS*Recv.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverRecvs) {
    int Rank = Pair.first;
    receiver_send_recv &Recv = Pair.second;
    CreateReceiverData(Recv.Data, Recv.Count);
    Irecv(Recv.Data.Points.Data(), MAX_DIMS*Recv.Count, MPI_INT, Rank, 1, Comm);
    Irecv(Recv.Data.SourceGridIDs.Data(), Recv.Count, MPI_INT, Rank, 1, Comm);
    Irecv(Recv.Data.SourceCells.Data(), MAX_DIMS*Recv.Count, MPI_INT, Rank, 1, Comm);
    
  }

  for (auto &Pair : DonorSends) {
    int Rank = Pair.first;
    donor_send_recv &Send = Pair.second;
    Isend(Send.Data.Extents.Data(), 2*MAX_DIMS*Send.Count, MPI_INT, Rank, 0, Comm);
    Isend(Send.Data.Coords.Data(), MAX_DIMS*Send.Count, MPI_DOUBLE, Rank, 0, Comm);
    Isend(Send.Data.InterpCoefs.Data(), MAX_DIMS*Send.MaxSize*Send.Count, MPI_DOUBLE, Rank, 0,
      Comm);
    Isend(Send.Data.DestinationGridIDs.Data(), Send.Count, MPI_INT, Rank, 0, Comm);
    Isend(Send.Data.DestinationPoints.Data(), MAX_DIMS*Send.Count, MPI_INT, Rank, 0, Comm);
  }

  for (auto &Pair : ReceiverSends) {
    int Rank = Pair.first;
    receiver_send_recv &Send = Pair.second;
    Isend(Send.Data.Points.Data(), MAX_DIMS*Send.Count, MPI_INT, Rank, 1, Comm);
    Isend(Send.Data.SourceGridIDs.Data(), Send.Count, MPI_INT, Rank, 1, Comm);
    Isend(Send.Data.SourceCells.Data(), MAX_DIMS*Send.Count, MPI_INT, Rank, 1, Comm);
  }

  MPI_Waitall(Requests.size(), Requests.data(), MPI_STATUSES_IGNORE);

  Requests.clear();

  DonorSends.clear();
  ReceiverSends.clear();

  long long NumLocalDonors = 0;
  int MaxSize = 0;
  long long NumLocalReceivers = 0;

  for (auto &Pair : DonorRecvs) {
    donor_send_recv &Recv = Pair.second;
    NumLocalDonors += Recv.Count;
    MaxSize = std::max(MaxSize, Recv.MaxSize);
  }

  for (auto &Pair : ReceiverRecvs) {
    receiver_send_recv &Recv = Pair.second;
    NumLocalReceivers += Recv.Count;
  }

  CreateDonorData(DonorData, NumLocalDonors, MaxSize);

  long long iNext;

  iNext = 0;
  for (auto &Pair : DonorRecvs) {
    donor_send_recv &Recv = Pair.second;
    for (long long iDonor = 0; iDonor < Recv.Count; ++iDonor) {
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorData.Extents(0,iDim,iNext) = Recv.Data.Extents(0,iDim,iDonor);
        DonorData.Extents(1,iDim,iNext) = Recv.Data.Extents(1,iDim,iDonor);
      }
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorData.Coords(iDim,iNext) = Recv.Data.Coords(iDim,iDonor);
      }
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        int Size = Recv.Data.Extents(1,iDim,iDonor) - Recv.Data.Extents(0,iDim,iDonor);
        for (int iPoint = 0; iPoint < Size; ++iPoint) {
          DonorData.InterpCoefs(iDim,iPoint,iNext) = Recv.Data.InterpCoefs(iDim,iPoint,iDonor);
        }
      }
      DonorData.DestinationGridIDs(iNext) = Recv.Data.DestinationGridIDs(iDonor);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorData.DestinationPoints(iDim,iNext) = Recv.Data.DestinationPoints(iDim,iDonor);
      }
      ++iNext;
    }
  }

  CreateReceiverData(ReceiverData, NumLocalReceivers);

  iNext = 0;
  for (auto &Pair : ReceiverRecvs) {
    receiver_send_recv &Recv = Pair.second;
    for (long long iReceiver = 0; iReceiver < Recv.Count; ++iReceiver) {
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        ReceiverData.Points(iDim,iNext) = Recv.Data.Points(iDim,iReceiver);
      }
      ReceiverData.SourceGridIDs(iNext) = Recv.Data.SourceGridIDs(iReceiver);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        ReceiverData.SourceCells(iDim,iNext) = Recv.Data.SourceCells(iDim,iReceiver);
      }
      ++iNext;
    }
  }

  core::EndProfile(Profiler, SendDataTime);

}

void ImportConnectivityData(int NumGrids, int NumLocalGrids, const std::vector<int> &LocalGridIDs,
  const std::vector<donor_data> &LocalDonorData, const std::vector<receiver_data>
  &LocalReceiverData, domain &Domain) {

  std::string DomainName;
  GetDomainName(Domain, DomainName);

  const core::comm &Comm = core::GetDomainComm(Domain);

  MPI_Barrier(Comm);

  core::logger &Logger = core::GetDomainLogger(Domain);

  core::LogStatus(Logger, Comm.Rank() == 0, 0, "Importing connectivity data into domain %s...",
    DomainName);

  std::vector<long long *> NumConnections(NumGrids);
  std::vector<long long> NumConnectionsData(NumGrids*NumGrids, 0);
  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    NumConnections[iGrid] = NumConnectionsData.data() + iGrid*NumGrids;
  }

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    int jGrid = LocalGridIDs[iLocalGrid]-1;
    for (long long iReceiver = 0; iReceiver < LocalReceiverData[iLocalGrid].Count; ++iReceiver) {
      int iGrid = LocalReceiverData[iLocalGrid].SourceGridIDs(iReceiver)-1;
      ++NumConnections[iGrid][jGrid];
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, NumConnectionsData.data(), NumGrids*NumGrids, MPI_LONG_LONG, MPI_SUM,
    Comm);

  struct local_connectivity {
    connectivity *Connectivity;
    connectivity_d *Donors;
    connectivity_r *Receivers;
  };

  std::map<int, local_connectivity> LocalConnectivities;

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections[iGrid][jGrid] > 0) {
        int DonorGridID = iGrid+1;
        int ReceiverGridID = jGrid+1;
        bool DonorGridIsLocal = RankHasGrid(Domain, DonorGridID);
        bool ReceiverGridIsLocal = RankHasGrid(Domain, ReceiverGridID);
        if (DonorGridIsLocal || ReceiverGridIsLocal) {
          int iPair = NumGrids*iGrid + jGrid;
          local_connectivity &LocalConnectivity = LocalConnectivities[iPair];
          EditConnectivityLocal(Domain, DonorGridID, ReceiverGridID,
            LocalConnectivity.Connectivity);
          connectivity &Connectivity = *LocalConnectivity.Connectivity;
          LocalConnectivity.Donors = nullptr;
          if (DonorGridIsLocal) {
            EditConnectivityDonorSideLocal(Connectivity, LocalConnectivity.Donors);
          } else {
            EditConnectivityDonorSideRemote(Connectivity);
          }
          LocalConnectivity.Receivers = nullptr;
          if (ReceiverGridIsLocal) {
            EditConnectivityReceiverSideLocal(Connectivity, LocalConnectivity.Receivers);
          } else {
            EditConnectivityReceiverSideRemote(Connectivity);
          }
        } else {
          EditConnectivityRemote(Domain, DonorGridID, ReceiverGridID);
        }
      }
    }
  }

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {

    int GridID = LocalGridIDs[iLocalGrid];
    int iGrid = GridID-1;

    const grid *Grid;
    GetGrid(Domain, GridID, Grid);
    const core::comm &GridComm = core::GetGridComm(*Grid);

    int NumReceiverGrids = 0;
    int NumDonorGrids = 0;
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections[iGrid][jGrid] > 0) {
        ++NumReceiverGrids;
      }
      if (NumConnections[jGrid][iGrid] > 0) {
        ++NumDonorGrids;
      }
    }

    std::vector<connectivity_d *> Donors;
    std::vector<connectivity_r *> Receivers;

    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections[iGrid][jGrid] > 0) {
        int iPair = NumGrids*iGrid + jGrid;
        local_connectivity &LocalConnectivity = LocalConnectivities[iPair];
        Donors.push_back(LocalConnectivity.Donors);
      }
      if (NumConnections[jGrid][iGrid] > 0) {
        int iPair = NumGrids*jGrid + iGrid;
        local_connectivity &LocalConnectivity = LocalConnectivities[iPair];
        Receivers.push_back(LocalConnectivity.Receivers);
      }
    }

    ImportDonors(LocalDonorData[iLocalGrid], GridComm, NumReceiverGrids, Donors);
    ImportReceivers(LocalReceiverData[iLocalGrid], GridComm, NumDonorGrids, Receivers);

  }

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections[iGrid][jGrid] > 0) {
        int DonorGridID = iGrid+1;
        int ReceiverGridID = jGrid+1;
        bool DonorGridIsLocal = RankHasGrid(Domain, DonorGridID);
        bool ReceiverGridIsLocal = RankHasGrid(Domain, ReceiverGridID);
        if (DonorGridIsLocal || ReceiverGridIsLocal) {
          int iPair = NumGrids*iGrid + jGrid;
          local_connectivity &LocalConnectivity = LocalConnectivities[iPair];
          connectivity &Connectivity = *LocalConnectivity.Connectivity;
          if (DonorGridIsLocal) {
            ReleaseConnectivityDonorSideLocal(Connectivity, LocalConnectivity.Donors);
          } else {
            ReleaseConnectivityDonorSideRemote(Connectivity);
          }
          if (ReceiverGridIsLocal) {
            ReleaseConnectivityReceiverSideLocal(Connectivity, LocalConnectivity.Receivers);
          } else {
            ReleaseConnectivityReceiverSideRemote(Connectivity);
          }
          ReleaseConnectivityLocal(Domain, DonorGridID, ReceiverGridID,
            LocalConnectivity.Connectivity);
        } else {
          ReleaseConnectivityRemote(Domain, DonorGridID, ReceiverGridID);
        }
      }
    }
  }

  MPI_Barrier(Comm);

  core::LogStatus(Logger, Comm.Rank() == 0, 0, "Finished importing connectivity data into domain %s.",
    DomainName);

}

void ImportDonors(const donor_data &GridDonors, const core::comm &Comm, int NumDestinationGrids,
  const std::vector<connectivity_d *> &OverkitDonors) {

  if (NumDestinationGrids == 0) return;

  std::vector<int> DestinationGridIDs(NumDestinationGrids);

  for (int iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    connectivity_d &Donors = *OverkitDonors[iDestinationGrid];
    GetConnectivityDonorSideDestinationGridID(Donors, DestinationGridIDs[iDestinationGrid]);
  }

  struct donor_edit {
    long long Count;
    int MaxSize;
    int *Extents[2][MAX_DIMS];
//     static_array<int *,2*MAX_DIMS,2> Extents;
    double *Coords[MAX_DIMS];
//     static_array<double *,MAX_DIMS> Coords;
    array<double *,2> InterpCoefs;
    int *DestinationPoints[MAX_DIMS];
//     static_array<int *,MAX_DIMS> DestinationPoints;
    donor_edit():
      Count(0),
      MaxSize(0)
    {}
  };

  std::map<int, donor_edit> DonorEdits;

  for (int iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    int DestinationGridID = DestinationGridIDs[iDestinationGrid];
    DonorEdits.emplace(DestinationGridID, donor_edit());
  }

  for (long long iDonor = 0; iDonor < GridDonors.Count; ++iDonor) {
    int DestinationGridID = GridDonors.DestinationGridIDs[iDonor];
    donor_edit &DonorEdit = DonorEdits[DestinationGridID];
    ++DonorEdit.Count;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      int Size = GridDonors.Extents(1,iDim,iDonor) - GridDonors.Extents(0,iDim,iDonor);
      DonorEdit.MaxSize = std::max(DonorEdit.MaxSize, Size);
    }
  }

  for (auto &Pair : DonorEdits) {
    donor_edit &DonorEdit = Pair.second;
    MPI_Allreduce(MPI_IN_PLACE, &DonorEdit.MaxSize, 1, MPI_INT, MPI_MAX, Comm);
//     DonorEdit.Extents.Resize({{2,MAX_DIMS}});
//     DonorEdit.Coords.Resize({MAX_DIMS});
    DonorEdit.InterpCoefs.Resize({{MAX_DIMS,DonorEdit.MaxSize}});
//     DonorEdit.DestinationPoints.Resize({MAX_DIMS});
  }

  for (int iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    int DestinationGridID = DestinationGridIDs[iDestinationGrid];
    connectivity_d &Donors = *OverkitDonors[iDestinationGrid];
    donor_edit &DonorEdit = DonorEdits[DestinationGridID];
    int MaxSize = DonorEdit.MaxSize;
    ResizeDonors(Donors, DonorEdit.Count, MaxSize);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      EditDonorExtents(Donors, iDim, DonorEdit.Extents[0][iDim], DonorEdit.Extents[1][iDim]);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      EditDonorCoords(Donors, iDim, DonorEdit.Coords[iDim]);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (int iPoint = 0; iPoint < MaxSize; ++iPoint) {
        EditDonorInterpCoefs(Donors, iDim, iPoint, DonorEdit.InterpCoefs(iDim,iPoint));
      }
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      EditDonorDestinations(Donors, iDim, DonorEdit.DestinationPoints[iDim]);
    }
    // Reset count to 0 for filling in data
    DonorEdit.Count = 0;
  }

  for (long long iDonor = 0; iDonor < GridDonors.Count; ++iDonor) {
    int DestinationGridID = GridDonors.DestinationGridIDs(iDonor);
    donor_edit &DonorEdit = DonorEdits[DestinationGridID];
    long long &iNext = DonorEdit.Count;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      DonorEdit.Extents[0][iDim][iNext] = GridDonors.Extents(0,iDim,iDonor);
      DonorEdit.Extents[1][iDim][iNext] = GridDonors.Extents(1,iDim,iDonor);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      DonorEdit.Coords[iDim][iNext] = GridDonors.Coords(iDim,iDonor);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      int Size = GridDonors.Extents(1,iDim,iDonor)-GridDonors.Extents(0,iDim,iDonor);
      for (int iPoint = 0; iPoint < Size; ++iPoint) {
        DonorEdit.InterpCoefs(iDim,iPoint)[iNext] = GridDonors.InterpCoefs(iDim,iPoint,iDonor);
      }
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      DonorEdit.DestinationPoints[iDim][iNext] = GridDonors.DestinationPoints(iDim,iDonor);
    }
    ++iNext;
  }

  for (int iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    int DestinationGridID = DestinationGridIDs[iDestinationGrid];
    connectivity_d &Donors = *OverkitDonors[iDestinationGrid];
    donor_edit &DonorEdit = DonorEdits[DestinationGridID];
    int MaxSize = DonorEdit.MaxSize;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReleaseDonorExtents(Donors, iDim, DonorEdit.Extents[0][iDim], DonorEdit.Extents[1][iDim]);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReleaseDonorCoords(Donors, iDim, DonorEdit.Coords[iDim]);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (int iPoint = 0; iPoint < MaxSize; ++iPoint) {
        ReleaseDonorInterpCoefs(Donors, iDim, iPoint, DonorEdit.InterpCoefs(iDim,iPoint));
      }
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReleaseDonorDestinations(Donors, iDim, DonorEdit.DestinationPoints[iDim]);
    }
  }

}

void ImportReceivers(const receiver_data &GridReceivers, const core::comm &Comm, int NumSourceGrids,
  const std::vector<connectivity_r *> &OverkitReceivers) {

  if (NumSourceGrids == 0) return;

  std::vector<int> SourceGridIDs(NumSourceGrids);

  for (int iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    connectivity_r &Receivers = *OverkitReceivers[iSourceGrid];
    GetConnectivityReceiverSideSourceGridID(Receivers, SourceGridIDs[iSourceGrid]);
  }

  struct receiver_edit {
    long long Count;
    int *Points[MAX_DIMS];
//     static_array<int *,MAX_DIMS> Points;
    int *SourceCells[MAX_DIMS];
//     static_array<int *,MAX_DIMS> SourceCells;
    receiver_edit():
      Count(0)
    {}
  };

  std::map<int, receiver_edit> ReceiverEdits;

  for (int iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    int SourceGridID = SourceGridIDs[iSourceGrid];
    ReceiverEdits.emplace(SourceGridID, receiver_edit());
  }

  for (long long iReceiver = 0; iReceiver < GridReceivers.Count; ++iReceiver) {
    int SourceGridID = GridReceivers.SourceGridIDs(iReceiver);
    receiver_edit &ReceiverEdit = ReceiverEdits[SourceGridID];
    ++ReceiverEdit.Count;
  }

//   for (auto &Pair : ReceiverEdits) {
//     receiver_edit &ReceiverEdit = Pair.second;
//     ReceiverEdit.Points.Resize({MAX_DIMS});
//     ReceiverEdit.SourceCells.Resize({MAX_DIMS});
//   }

  for (int iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    int SourceGridID = SourceGridIDs[iSourceGrid];
    connectivity_r &Receivers = *OverkitReceivers[iSourceGrid];
    receiver_edit &ReceiverEdit = ReceiverEdits[SourceGridID];
    ResizeReceivers(Receivers, ReceiverEdit.Count);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      EditReceiverPoints(Receivers, iDim, ReceiverEdit.Points[iDim]);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      EditReceiverSources(Receivers, iDim, ReceiverEdit.SourceCells[iDim]);
    }
    // Reset count to 0 for filling in data
    ReceiverEdit.Count = 0;
  }

  for (long long iReceiver = 0; iReceiver < GridReceivers.Count; ++iReceiver) {
    int SourceGridID = GridReceivers.SourceGridIDs(iReceiver);
    receiver_edit &ReceiverEdit = ReceiverEdits[SourceGridID];
    int iNext = ReceiverEdit.Count;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReceiverEdit.Points[iDim][iNext] = GridReceivers.Points(iDim,iReceiver);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReceiverEdit.SourceCells[iDim][iNext] = GridReceivers.SourceCells(iDim,iReceiver);
    }
    ++ReceiverEdit.Count;
  }

  for (int iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    int SourceGridID = SourceGridIDs[iSourceGrid];
    connectivity_r &Receivers = *OverkitReceivers[iSourceGrid];
    receiver_edit &ReceiverEdit = ReceiverEdits[SourceGridID];
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReleaseReceiverPoints(Receivers, iDim, ReceiverEdit.Points[iDim]);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      ReleaseReceiverSources(Receivers, iDim, ReceiverEdit.SourceCells[iDim]);
    }
  }

}

void CreateDonorData(donor_data &Data, long long Count, int MaxSize) {

  Data.Count = Count;
  Data.MaxSize = MaxSize;

  Data.Extents.Resize({{2,MAX_DIMS,Count}});
  Data.Coords.Resize({{MAX_DIMS,Count}});
  Data.InterpCoefs.Resize({{MAX_DIMS,MaxSize,Count}});
  Data.DestinationGridIDs.Resize({Count});
  Data.DestinationPoints.Resize({{MAX_DIMS,Count}});

}

void CreateReceiverData(receiver_data &Data, long long Count) {

  Data.Count = Count;

  Data.Points.Resize({{MAX_DIMS,Count}});
  Data.SourceGridIDs.Resize({Count});
  Data.SourceCells.Resize({{MAX_DIMS,Count}});

}

void CreateConnectionData(connection_data &Data, long long Count) {

  Data.Count = Count;

  Data.ConnectionIDs.Resize({Count});
  Data.DonorGridIDs.Resize({Count});
  Data.DonorCells.Resize({{MAX_DIMS,Count}});
  Data.ReceiverGridIDs.Resize({Count});
  Data.ReceiverPoints.Resize({{MAX_DIMS,Count}});

}

long long BinDivide(long long N, int NumBins) {

  long long NumBinsLongLong = (long long)(NumBins);

  return (N + NumBinsLongLong - 1)/NumBinsLongLong;

}

void Chunkify(long long Count, int MaxChunks, long long TargetChunkSize, int Adjust,
  int &ChunkRankInterval, int &NumChunks, long long &ChunkSize) {

  ChunkRankInterval = 1;
  NumChunks = MaxChunks;

  while (BinDivide(Count, NumChunks) < TargetChunkSize && NumChunks > 1) {
    ChunkRankInterval = int(std::min((long long)(ChunkRankInterval) << 1, (long long)(MaxChunks)));
    NumChunks = BinDivide((long long)(MaxChunks), ChunkRankInterval);
  }

  if (Adjust > 0) {
    int RemainingAdjustAmount = Adjust;
    while (RemainingAdjustAmount > 0 && NumChunks < MaxChunks) {
      ChunkRankInterval = std::max(ChunkRankInterval >> 1, 1);
      NumChunks = BinDivide((long long)(MaxChunks), ChunkRankInterval);
      --RemainingAdjustAmount;
    }
  } else if (Adjust < 0) {
    int RemainingAdjustAmount = -Adjust;
    while (RemainingAdjustAmount > 0 && NumChunks > 1) {
      ChunkRankInterval = int(std::min((long long)(ChunkRankInterval) << 1, (long long)(MaxChunks)));
      NumChunks = BinDivide((long long)(MaxChunks), ChunkRankInterval);
      --RemainingAdjustAmount;
    }
  }

  ChunkSize = BinDivide(Count, NumChunks);

}

int File_read_all_endian(MPI_File File, void *Buffer, int Count, MPI_Datatype DataType,
  endian Endian, MPI_Status *Status, core::profiler &Profiler, const core::comm &Comm) {

  int MPIIOReadTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Read");
  core::StartProfileSync(Profiler, MPIIOReadTime, Comm);
  int MPIError = MPI_File_read_all(File, Buffer, Count, DataType, Status);
  core::EndProfile(Profiler, MPIIOReadTime);

  if (MPIError == MPI_SUCCESS) {
    if (Endian != MachineEndian()) {
      int Size;
      MPI_Type_size(DataType, &Size);
      SwapEndian(Buffer, Size, Count);
    }
  }

  return MPIError;

}

int File_read_at_endian(MPI_File File, MPI_Offset Offset, void *Buffer, int Count,
  MPI_Datatype DataType, endian Endian, MPI_Status *Status, core::profiler &Profiler) {

  int MPIIOReadTime = core::GetProfilerTimerID(Profiler, "XINTOUT::Import::Read::MPI-IO::Read");
  core::StartProfile(Profiler, MPIIOReadTime);
  int MPIError = MPI_File_read_at(File, Offset, Buffer, Count, DataType, Status);
  core::EndProfile(Profiler, MPIIOReadTime);

  if (MPIError == MPI_SUCCESS) {
    if (Endian != MachineEndian()) {
      int Size;
      MPI_Type_size(DataType, &Size);
      SwapEndian(Buffer, Size, Count);
    }
  }

  return MPIError;

}

endian MachineEndian() {

  elem<unsigned char,2> EndianTest = {1,0};

  if(*reinterpret_cast<short *>(EndianTest.Data()) == 1) {
    return endian::LITTLE;
  } else {
    return endian::BIG;
  }

}

void SwapEndian(void *Data, int ElementSize, int NumElements) {

  auto DataBytes = static_cast<unsigned char *>(Data);

  std::vector<unsigned char> Element(ElementSize);

  for (int i = 0; i < NumElements; ++i) {
    for (int j = 0; j < ElementSize; ++j) {
      Element[j] = DataBytes[ElementSize*i+j];
    }
    for (int j = 0; j < ElementSize; ++j) {
      DataBytes[ElementSize*i+j] = Element[ElementSize-j-1];
    }
  }

}

}}
