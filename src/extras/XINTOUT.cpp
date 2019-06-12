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
#include "ovk/core/ConnectivityComponent.hpp"
#include "ovk/core/ConnectivityM.hpp"
#include "ovk/core/ConnectivityN.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Error.hpp"
#include "ovk/core/IDMap.hpp"
#include "ovk/core/IDSet.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Optional.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/ScopeGuard.hpp"
#include "ovk/core/TextProcessing.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <algorithm>
#include <limits>
#include <string>
#include <utility>

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
  context *Context;
  int ID;
  std::string Name;
  int NumDims;
  core::comm_view Comm;
  tuple<int> GlobalSize;
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
  context *Context;
  int NumDims;
  core::comm_view Comm;
  int NumGrids;
  int NumLocalGrids;
  array<xintout_grid> Grids;
  xintout_connections Connections;
};

constexpr int IMPORT_TIME = core::profiler::XINTOUT_IMPORT_TIME;
constexpr int IMPORT_READ_TIME = core::profiler::XINTOUT_IMPORT_READ_TIME;
constexpr int IMPORT_READ_MPI_IO_OPEN_TIME = core::profiler::XINTOUT_IMPORT_READ_MPI_IO_OPEN_TIME;
constexpr int IMPORT_READ_MPI_IO_CLOSE_TIME = core::profiler::XINTOUT_IMPORT_READ_MPI_IO_CLOSE_TIME;
constexpr int IMPORT_READ_MPI_IO_READ_TIME = core::profiler::XINTOUT_IMPORT_READ_MPI_IO_READ_TIME;
constexpr int IMPORT_READ_MPI_IO_OTHER_TIME = core::profiler::XINTOUT_IMPORT_READ_MPI_IO_OTHER_TIME;
constexpr int IMPORT_MATCH_TIME = core::profiler::XINTOUT_IMPORT_MATCH_TIME;
constexpr int IMPORT_MATCH_MAP_TO_BINS_TIME = core::profiler::XINTOUT_IMPORT_MATCH_MAP_TO_BINS_TIME;
constexpr int IMPORT_MATCH_HANDSHAKE_TIME = core::profiler::XINTOUT_IMPORT_MATCH_HANDSHAKE_TIME;
constexpr int IMPORT_MATCH_SEND_TO_BINS_TIME = core::profiler::XINTOUT_IMPORT_MATCH_SEND_TO_BINS_TIME;
constexpr int IMPORT_MATCH_FILL_CONNECTION_DATA_TIME = core::profiler::XINTOUT_IMPORT_MATCH_FILL_CONNECTION_DATA_TIME;
constexpr int IMPORT_MATCH_RECV_FROM_BINS_TIME = core::profiler::XINTOUT_IMPORT_MATCH_RECV_FROM_BINS_TIME;
constexpr int IMPORT_MATCH_UNPACK_TIME = core::profiler::XINTOUT_IMPORT_MATCH_UNPACK_TIME;
constexpr int IMPORT_DISTRIBUTE_TIME = core::profiler::XINTOUT_IMPORT_DISTRIBUTE_TIME;
constexpr int IMPORT_DISTRIBUTE_MAP_TO_BINS_TIME = core::profiler::XINTOUT_IMPORT_DISTRIBUTE_MAP_TO_BINS_TIME;
constexpr int IMPORT_DISTRIBUTE_RETRIEVE_BINS_TIME = core::profiler::XINTOUT_IMPORT_DISTRIBUTE_RETRIEVE_BINS_TIME;
constexpr int IMPORT_DISTRIBUTE_FIND_RANKS_TIME = core::profiler::XINTOUT_IMPORT_DISTRIBUTE_FIND_RANKS_TIME;
constexpr int IMPORT_DISTRIBUTE_HANDSHAKE_TIME = core::profiler::XINTOUT_IMPORT_DISTRIBUTE_HANDSHAKE_TIME;
constexpr int IMPORT_DISTRIBUTE_SEND_DATA_TIME = core::profiler::XINTOUT_IMPORT_DISTRIBUTE_SEND_DATA_TIME;
constexpr int IMPORT_SET_CONNECTIVITIES_TIME = core::profiler::XINTOUT_IMPORT_SET_CONNECTIVITIES_TIME;

xintout CreateXINTOUT(context &Context, int NumDims, core::comm_view Comm, int NumGrids, int
  NumLocalGrids, const array<int> &LocalGridIDs, const array<std::string> &LocalGridNames, const
  array<core::comm_view> &LocalGridComms, const array<tuple<int>> &LocalGridGlobalSizes);

xintout_grid CreateXINTOUTGrid(context &Context, int ID, const std::string &Name, int NumDims,
  core::comm_view Comm, const tuple<int> &GlobalSize);

void ReadXINTOUT(xintout &XINTOUT, const std::string &HOPath, const std::string &XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo);

void ReadGlobalInfo(const xintout &XINTOUT, const std::string &HOPath, endian &Endian,
  xintout_format &Format, bool &WithIBlank);

bool DetectFormat(MPI_File HOFile, endian &Endian, xintout_format &Format, core::profiler
  &Profiler);

void ReadGridInfo(const xintout_grid &XINTOUTGrid, const std::string &HOPath, const std::string
  &XPath, long long &NumDonors, long long &NumReceivers, long long &StartingConnectionID, MPI_Offset
  &HODonorCellsOffset, MPI_Offset &HODonorCoordsOffset, MPI_Offset &HOReceiverPointsOffset,
  MPI_Offset &HOReceiverConnectionIDsOffset, MPI_Offset &XDonorSizesOffset, MPI_Offset
  &XDonorInterpCoefsOffset, endian Endian, xintout_format Format, bool WithIBlank);

void ReadDonors(xintout_grid &XINTOUTGrid, const std::string &HOPath, const std::string &XPath,
  long long NumDonors, long long StartingConnectionID, MPI_Offset HOCellsOffset,
  MPI_Offset HOCoordsOffset, MPI_Offset XSizesOffset, MPI_Offset XInterpCoefsOffset,
  endian Endian, int ReadGranularityAdjust, MPI_Info MPIInfo);

void ReadReceivers(xintout_grid &XINTOUTGrid, const std::string &HOPath, long long NumReceivers,
  MPI_Offset HOPointsOffset, MPI_Offset HOConnectionIDsOffset, endian Endian, xintout_format Format,
  int ReadGranularityAdjust, MPI_Info MPIInfo);

void MatchDonorsAndReceivers(xintout &XINTOUT);

void DistributeConnectivityData(const xintout &XINTOUT, const array<const grid *> &LocalGrids,
  array<donor_data> &LocalDonorData, array<receiver_data> &LocalReceiverData);

void DistributeGridConnectivityData(const xintout_grid &XINTOUTGrid, const grid &Grid,
  donor_data &DonorData, receiver_data &ReceiverData);

void ImportConnectivityData(int NumGrids, int NumLocalGrids, const array<int> &LocalGridIDs,
  const array<donor_data> &LocalDonorData, const array<receiver_data> &LocalReceiverData, domain
  &Domain, int ConnectivityComponentID);

void ImportDonors(const donor_data &GridDonors, core::comm_view Comm, const array<edit_handle<
  connectivity_m>> &ConnectivityMs);
void ImportReceivers(const receiver_data &GridReceivers, core::comm_view Comm, const array<
  edit_handle<connectivity_n>> &ConnectivityNs);

void CreateDonorData(donor_data &Data, long long Count, int MaxSize);
void CreateReceiverData(receiver_data &Data, long long Count);
void CreateConnectionData(connection_data &Data, long long Count);

long long BinDivide(long long N, int NumChunks);
void Chunkify(long long Count, int MaxChunks, long long TargetChunkSize, int Adjust,
  int &ChunkRankInterval, int &NumChunks, long long &ChunkSize);

int File_read_all_endian(MPI_File File, void *Buffer, int Count, MPI_Datatype DataType,
  endian Endian, MPI_Status *Status, core::profiler &Profiler, core::comm_view Comm);
int File_read_at_endian(MPI_File File, MPI_Offset Offset, void *Buffer, int Count,
  MPI_Datatype DataType, endian Endian, MPI_Status *Status, core::profiler &Profiler);

endian MachineEndian();
void SwapEndian(void *Data, int ElementSize, int NumElements);

}

void ImportXINTOUT(domain &Domain, int ConnectivityComponentID, const std::string &HOPath, const
  std::string &XPath, int ReadGranularityAdjust, MPI_Info MPIInfo) {

  context &Context = Domain.Context();
  core::profiler &Profiler = Context.core_Profiler();

  int NumDims = Domain.Dimension();
  int NumGrids = Domain.GridCount();

  core::comm_view Comm = Domain.core_Comm();

  MPI_Barrier(Comm);

  if (NumGrids > 0) {

    int NumLocalGrids = Domain.LocalGridCount();

    array<const grid *> LocalGrids({NumLocalGrids});
    array<int> LocalGridIDs({NumLocalGrids});
    array<std::string> LocalGridNames({NumLocalGrids});
    array<core::comm_view> LocalGridComms({NumLocalGrids});
    array<tuple<int>> LocalGridGlobalSizes({NumLocalGrids});
    int iLocalGrid = 0;
    for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
      int GridID = iGrid+1;
      if (Domain.GridIsLocal(GridID)) {
        const grid &Grid = Domain.Grid(GridID);
        LocalGrids(iLocalGrid) = &Grid;
        LocalGridIDs(iLocalGrid) = GridID;
        LocalGridNames(iLocalGrid) = Grid.Name();
        LocalGridComms(iLocalGrid) = Grid.core_Comm();
        LocalGridGlobalSizes(iLocalGrid) = Grid.Size();
        ++iLocalGrid;
      }
    }

    Profiler.StartSync(IMPORT_TIME, Comm);
    // ReadXINTOUT might throw
    auto StopImportTimer = core::OnScopeExit([&] { Profiler.Stop(IMPORT_TIME); });

    xintout XINTOUT = CreateXINTOUT(Context, NumDims, Comm, NumGrids, NumLocalGrids, LocalGridIDs,
      LocalGridNames, LocalGridComms, LocalGridGlobalSizes);

    Profiler.StartSync(IMPORT_READ_TIME, Comm);
    // ReadXINTOUT might throw
    auto StopReadTimer = core::OnScopeExit([&] { Profiler.Stop(IMPORT_READ_TIME); });

    ReadXINTOUT(XINTOUT, HOPath, XPath, ReadGranularityAdjust, MPIInfo);

    Profiler.Stop(IMPORT_READ_TIME);

    StopImportTimer.Dismiss();
    StopReadTimer.Dismiss();

    Profiler.StartSync(IMPORT_MATCH_TIME, Comm);
    MatchDonorsAndReceivers(XINTOUT);
    Profiler.Stop(IMPORT_MATCH_TIME);

    Profiler.StartSync(IMPORT_DISTRIBUTE_TIME, Comm);
    array<donor_data> LocalDonors({NumLocalGrids});
    array<receiver_data> LocalReceivers({NumLocalGrids});
    DistributeConnectivityData(XINTOUT, LocalGrids, LocalDonors, LocalReceivers);
    Profiler.Stop(IMPORT_DISTRIBUTE_TIME);

    Profiler.StartSync(IMPORT_SET_CONNECTIVITIES_TIME, Comm);
    ImportConnectivityData(NumGrids, NumLocalGrids, LocalGridIDs, LocalDonors, LocalReceivers,
      Domain, ConnectivityComponentID);
    Profiler.Stop(IMPORT_SET_CONNECTIVITIES_TIME);

    Profiler.Stop(IMPORT_TIME);

  }

}

void ImportXINTOUT(domain &Domain, int ConnectivityComponentID, const std::string &HOPath, const
  std::string &XPath, int ReadGranularityAdjust, MPI_Info MPIInfo, error &Error) {

  Error = error::NONE;

  try {
    ImportXINTOUT(Domain, ConnectivityComponentID, HOPath, XPath, ReadGranularityAdjust, MPIInfo);
  } catch (const exception &Exception) {
    Error = Exception.Error();
  }

}

void ExportXINTOUT(const domain &Domain, int ConnectivityComponentID, const std::string &HOPath,
  const std::string &XPath, xintout_format Format, endian Endian, int WriteGranularityAdjust,
  MPI_Info MPIInfo) {

  OVK_DEBUG_ASSERT(false, "ExportXINTOUT is not yet implemented.");

}

void ExportXINTOUT(const domain &Domain, int ConnectivityComponentID, const std::string &HOPath,
  const std::string &XPath, xintout_format Format, endian Endian, int WriteGranularityAdjust,
  MPI_Info MPIInfo, error &Error) {

  Error = error::NONE;

  try {
    ExportXINTOUT(Domain, ConnectivityComponentID, HOPath, XPath, Format, Endian,
      WriteGranularityAdjust, MPIInfo);
  } catch (const exception &Exception) {
    Error = Exception.Error();
  }

}

namespace {

xintout CreateXINTOUT(context &Context, int NumDims, core::comm_view Comm, int NumGrids, int
  NumLocalGrids, const array<int> &LocalGridIDs, const array<std::string> &LocalGridNames, const
  array<core::comm_view> &LocalGridComms, const array<tuple<int>>
  &LocalGridGlobalSizes) {

  xintout XINTOUT;

  XINTOUT.Comm = Comm;

  MPI_Barrier(XINTOUT.Comm);

  XINTOUT.NumDims = NumDims;

  XINTOUT.Context = &Context;

  XINTOUT.NumGrids = NumGrids;
  XINTOUT.NumLocalGrids = NumLocalGrids;

  XINTOUT.Grids.Resize({NumLocalGrids});
  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    XINTOUT.Grids(iLocalGrid) = CreateXINTOUTGrid(Context, LocalGridIDs(iLocalGrid),
      LocalGridNames(iLocalGrid), NumDims, LocalGridComms(iLocalGrid),
      LocalGridGlobalSizes(iLocalGrid));
  }

  XINTOUT.Connections.Count = 0;
  XINTOUT.Connections.BinSize = 0;
  XINTOUT.Connections.HasBin = false;

  MPI_Barrier(XINTOUT.Comm);

  return XINTOUT;

}

xintout_grid CreateXINTOUTGrid(context &Context, int ID, const std::string &Name, int NumDims,
  core::comm_view Comm, const tuple<int> &GlobalSize) {

  xintout_grid XINTOUTGrid;

  XINTOUTGrid.Comm = Comm;

  MPI_Barrier(XINTOUTGrid.Comm);

  XINTOUTGrid.Context = &Context;

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

  return XINTOUTGrid;

}

void ReadXINTOUT(xintout &XINTOUT, const std::string &HOPath, const std::string &XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo) {

  core::comm_view Comm = XINTOUT.Comm;

  MPI_Barrier(Comm);

  core::logger &Logger = XINTOUT.Context->core_Logger();

  int NumLocalGrids = XINTOUT.NumLocalGrids;

  Logger.LogStatus(Comm.Rank() == 0, 0, "Reading XINTOUT files '%s' and '%s'...", HOPath, XPath);

  endian Endian;
  xintout_format Format;
  bool WithIBlank;
  ReadGlobalInfo(XINTOUT, HOPath, Endian, Format, WithIBlank);

  error ReadError = error::NONE;

  try {

    for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {

      xintout_grid &XINTOUTGrid = XINTOUT.Grids(iLocalGrid);

      long long NumDonors, NumReceivers;
      long long DonorStartingConnectionID;
      MPI_Offset HODonorCellsOffset, HODonorCoordsOffset;
      MPI_Offset XDonorSizesOffset, XDonorInterpCoefsOffset;
      MPI_Offset HOReceiverPointsOffset, HOReceiverConnectionIDsOffset;
      ReadGridInfo(XINTOUTGrid, HOPath, XPath, NumDonors, NumReceivers, DonorStartingConnectionID,
        HODonorCellsOffset, HODonorCoordsOffset, HOReceiverPointsOffset,
        HOReceiverConnectionIDsOffset, XDonorSizesOffset, XDonorInterpCoefsOffset, Endian,
        Format, WithIBlank);

      ReadDonors(XINTOUTGrid, HOPath, XPath, NumDonors, DonorStartingConnectionID,
        HODonorCellsOffset, HODonorCoordsOffset, XDonorSizesOffset, XDonorInterpCoefsOffset, Endian,
        ReadGranularityAdjust, MPIInfo);

      ReadReceivers(XINTOUTGrid, HOPath, NumReceivers, HOReceiverPointsOffset,
        HOReceiverConnectionIDsOffset, Endian, Format, ReadGranularityAdjust, MPIInfo);

    }

  } catch (exception &Exception) {

    ReadError = Exception.Error();

  }

  core::SyncError(ReadError, Comm);
  core::CheckError(ReadError);

  MPI_Barrier(Comm);

  Logger.LogStatus(Comm.Rank() == 0, 0, "Done reading XINTOUT files.");

}

void ReadGlobalInfo(const xintout &XINTOUT, const std::string &HOPath, endian &Endian,
  xintout_format &Format, bool &WithIBlank) {

  core::comm_view Comm = XINTOUT.Comm;

  int NumGrids = XINTOUT.NumGrids;

  core::logger &Logger = XINTOUT.Context->core_Logger();
  core::profiler &Profiler = XINTOUT.Context->core_Profiler();

  error ReadGlobalInfoError = error::NONE;

  try {

    if (Comm.Rank() == 0) {

      MPI_File HOFile;
      MPI_Status Status;
      int MPIError;
      int ReadSize;

      Profiler.Start(IMPORT_READ_MPI_IO_OPEN_TIME);
      // MPI_File_open missing const qualifier for path string on some platforms
      MPIError = MPI_File_open(MPI_COMM_SELF, const_cast<char *>(HOPath.c_str()), MPI_MODE_RDONLY,
        MPI_INFO_NULL, &HOFile);
      Profiler.Stop(IMPORT_READ_MPI_IO_OPEN_TIME);
      if (MPIError != MPI_SUCCESS) {
        Logger.LogError(true, "Unable to open file '%s'.", HOPath);
        throw file_open_error();
      }
      auto CloseHO = core::OnScopeExit([&]() {
        Profiler.Start(IMPORT_READ_MPI_IO_CLOSE_TIME);
        MPI_File_close(&HOFile);
        Profiler.Stop(IMPORT_READ_MPI_IO_CLOSE_TIME);
      });

      bool Success = DetectFormat(HOFile, Endian, Format, Profiler);
      if (!Success) {
        Logger.LogError(true, "Unable to detect format of XINTOUT file '%s'.", HOPath);
        throw file_read_error();
      }

      int RecordWrapperSize = Format == xintout_format::STANDARD ? sizeof(int) : 0;

      // Both formats have a record wrapper at the beginning
      MPI_Offset HOHeaderOffset = sizeof(int);
      int NumGridsInFile;
      File_read_at_endian(HOFile, HOHeaderOffset, &NumGridsInFile, 1, MPI_INT, Endian, &Status,
        Profiler);
      MPI_Get_count(&Status, MPI_INT, &ReadSize);
      if (ReadSize < 1) {
        Logger.LogError(true, "Unable to read header of XINTOUT file '%s'.", HOPath);
        throw file_read_error();
      }
      if (NumGridsInFile != NumGrids) {
        Logger.LogError(Comm.Rank() == 0, "XINTOUT file '%s' has incorrect number of grids.",
          HOPath);
        throw file_read_error();
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
      tuple<int> GridSize;
      if (Format == xintout_format::STANDARD) {
        int Data[7];
        File_read_at_endian(HOFile, HOGridOffset, Data, 7, MPI_INT, Endian, &Status, Profiler);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        if (ReadSize < 7) {
          Logger.LogError(true, "Unable to read grid 1 header of XINTOUT file '%s'.", HOPath);
          throw file_read_error();
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
          Logger.LogError(true, "Unable to read grid 1 header of XINTOUT file '%s'.", HOPath);
          throw file_read_error();
        }
        HOGridOffset += 4*sizeof(long long);
        File_read_at_endian(HOFile, HOGridOffset, IntData, 3, MPI_INT, Endian, &Status, Profiler);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        if (ReadSize < 3) {
          Logger.LogError(true, "Unable to read grid 1 header of XINTOUT file '%s'.", HOPath);
          throw file_read_error();
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
            Logger.LogError(true, "Unable to detect whether XINTOUT file '%s' contains IBlank.",
              HOPath);
            throw file_read_error();
          }
        }
      } else {
        WithIBlank = false;
      }

    }

  } catch (const exception &Exception) {

    ReadGlobalInfoError = Exception.Error();

  }

  core::SyncError(ReadGlobalInfoError, Comm);
  core::CheckError(ReadGlobalInfoError);

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

}

bool DetectFormat(MPI_File HOFile, endian &Endian, xintout_format &Format, core::profiler &Profiler)
  {

  unsigned char InitialBytes[sizeof(int)];
  MPI_Status Status;
  Profiler.Start(IMPORT_READ_MPI_IO_READ_TIME);
  int MPIError = MPI_File_read_at(HOFile, 0, InitialBytes, sizeof(int), MPI_BYTE, &Status);
  Profiler.Stop(IMPORT_READ_MPI_IO_READ_TIME);
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

void ReadGridInfo(const xintout_grid &XINTOUTGrid, const std::string &HOPath, const std::string
  &XPath, long long &NumDonors, long long &NumReceivers, long long &StartingConnectionID, MPI_Offset
  &HODonorCellsOffset, MPI_Offset &HODonorCoordsOffset, MPI_Offset &HOReceiverPointsOffset,
  MPI_Offset &HOReceiverConnectionIDsOffset, MPI_Offset &XDonorSizesOffset, MPI_Offset
  &XDonorInterpCoefsOffset, endian Endian, xintout_format Format, bool WithIBlank) {

  int GridID = XINTOUTGrid.ID;
  int iGrid = GridID-1;

  core::comm_view Comm = XINTOUTGrid.Comm;

  core::logger &Logger = XINTOUTGrid.Context->core_Logger();
  core::profiler &Profiler = XINTOUTGrid.Context->core_Profiler();

  error ReadGridInfoError = error::NONE;

  try {

    if (Comm.Rank() == 0) {

      MPI_File HOFile, XFile;
      MPI_Status Status;
      int MPIError;
      int ReadSize;
      tuple<int> GridSize;
      long long NumInterpCoefs = 0;

      Profiler.Start(IMPORT_READ_MPI_IO_OPEN_TIME);
      // MPI_File_open missing const qualifier for path string on some platforms
      MPIError = MPI_File_open(MPI_COMM_SELF, const_cast<char *>(HOPath.c_str()), MPI_MODE_RDONLY,
        MPI_INFO_NULL, &HOFile);
      Profiler.Stop(IMPORT_READ_MPI_IO_OPEN_TIME);
      if (MPIError != MPI_SUCCESS) {
        Logger.LogError(true, "Unable to open file '%s'.", HOPath);
        throw file_open_error();
      }
      auto CloseHO = core::OnScopeExit([&]() {
        Profiler.Start(IMPORT_READ_MPI_IO_CLOSE_TIME);
        MPI_File_close(&HOFile);
        Profiler.Stop(IMPORT_READ_MPI_IO_CLOSE_TIME);
      });

      Profiler.Start(IMPORT_READ_MPI_IO_OPEN_TIME);
      // MPI_File_open missing const qualifier for path string on some platforms
      MPIError = MPI_File_open(MPI_COMM_SELF, const_cast<char *>(XPath.c_str()), MPI_MODE_RDONLY,
        MPI_INFO_NULL, &XFile);
      Profiler.Stop(IMPORT_READ_MPI_IO_OPEN_TIME);
      if (MPIError != MPI_SUCCESS) {
        Logger.LogError(true, "Unable to open file '%s'.", XPath);
        throw file_open_error();
      }
      auto CloseX = core::OnScopeExit([&]() {
        Profiler.Start(IMPORT_READ_MPI_IO_CLOSE_TIME);
        MPI_File_close(&XFile);
        Profiler.Stop(IMPORT_READ_MPI_IO_CLOSE_TIME);
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
            Logger.LogError(true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
              HOPath);
            throw file_read_error();
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
            Logger.LogError(true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
              HOPath);
            throw file_read_error();
          }
          HOGridOffset += 4*sizeof(long long);
          File_read_at_endian(HOFile, HOGridOffset, IntData, 3, MPI_INT, Endian, &Status, Profiler);
          MPI_Get_count(&Status, MPI_INT, &ReadSize);
          if (ReadSize < 3) {
            Logger.LogError(true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
              HOPath);
            throw file_read_error();
          }
          HOGridOffset += 3*sizeof(int);
          File_read_at_endian(HOFile, HOGridOffset, LongLongData+4, 1, MPI_LONG_LONG, Endian, &Status,
            Profiler);
          MPI_Get_count(&Status, MPI_LONG_LONG, &ReadSize);
          if (ReadSize < 1) {
            Logger.LogError(true, "Unable to read grid %i header of XINTOUT file '%s'.", OtherGridID,
              HOPath);
            throw file_read_error();
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
            Logger.LogError(true, "Unable to read grid %i data from XINTOUT file '%s'.", OtherGridID,
              XPath);
            throw file_read_error();
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
          Logger.LogError(true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
            HOPath);
          throw file_read_error();
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
          Logger.LogError(true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
            HOPath);
          throw file_read_error();
        }
        HOGridOffset += 4*sizeof(long long);
        File_read_at_endian(HOFile, HOGridOffset, IntData, 3, MPI_INT, Endian, &Status, Profiler);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        if (ReadSize < 3) {
          Logger.LogError(true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
            HOPath);
          throw file_read_error();
        }
        HOGridOffset += 3*sizeof(int);
        File_read_at_endian(HOFile, HOGridOffset, LongLongData+4, 1, MPI_LONG_LONG, Endian, &Status,
          Profiler);
        MPI_Get_count(&Status, MPI_LONG_LONG, &ReadSize);
        if (ReadSize < 1) {
          Logger.LogError(true, "Unable to read grid %i header of XINTOUT file '%s'.", GridID,
            HOPath);
          throw file_read_error();
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
        Logger.LogError(true, "Grid %i of XINTOUT file '%s' has incorrect size.", GridID, HOPath);
        throw file_read_error();
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

  } catch (const exception &Exception) {

    ReadGridInfoError = Exception.Error();

  }

  core::SyncError(ReadGridInfoError, Comm);
  core::CheckError(ReadGridInfoError);

  MPI_Bcast(&NumDonors, 1, MPI_LONG_LONG, 0, Comm);
  MPI_Bcast(&NumReceivers, 1, MPI_LONG_LONG, 0, Comm);
  MPI_Bcast(&StartingConnectionID, 1, MPI_LONG_LONG, 0, Comm);
  MPI_Bcast(&HODonorCellsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&HODonorCoordsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&HOReceiverPointsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&HOReceiverConnectionIDsOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&XDonorSizesOffset, 1, MPI_OFFSET, 0, Comm);
  MPI_Bcast(&XDonorInterpCoefsOffset, 1, MPI_OFFSET, 0, Comm);

}

void ReadDonors(xintout_grid &XINTOUTGrid, const std::string &HOPath, const std::string &XPath,
  long long NumDonors, long long StartingConnectionID, MPI_Offset HOCellsOffset,
  MPI_Offset HOCoordsOffset, MPI_Offset XSizesOffset, MPI_Offset XInterpCoefsOffset,
  endian Endian, int ReadGranularityAdjust, MPI_Info MPIInfo) {

  int GridID = XINTOUTGrid.ID;

  core::comm_view Comm = XINTOUTGrid.Comm;

  core::logger &Logger = XINTOUTGrid.Context->core_Logger();
  core::profiler &Profiler = XINTOUTGrid.Context->core_Profiler();

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

  if (Logger.LoggingDebug()) {
    std::string NumDonorsString = core::FormatNumber(NumDonors, "donors", "donor");
    std::string NumRanksString = core::FormatNumber(NumChunks, "I/O ranks", "I/O rank");
    Logger.LogDebug(Comm.Rank() == 0, 0, "Grid %s has %s; using %s.", XINTOUTGrid.Name,
      NumDonorsString, NumRanksString);
  }

  XINTOUTDonors.ChunkSize = ChunkSize;
  XINTOUTDonors.HasChunk = HasChunk;

  core::comm ChunkComm = core::CreateSubsetComm(Comm, HasChunk);

  error ReadDonorsError = error::NONE;

  try {

    if (HasChunk) {

      xintout_donor_chunk &Chunk = XINTOUTDonors.Chunk;

      long long LocalBegin = ChunkSize*ChunkComm.Rank();
      long long LocalEnd = std::min(ChunkSize*(ChunkComm.Rank()+1), NumDonors);
      long long NumLocalDonors = LocalEnd - LocalBegin;

      error ChunkSizeError = error::NONE;
      if (NumLocalDonors*sizeof(double) > std::numeric_limits<int>::max()) {
        ChunkSizeError = error::FILE_READ;
      }
      core::SyncError(ChunkSizeError, ChunkComm);
      if (ChunkSizeError != error::NONE) {
        Logger.LogError(ChunkComm.Rank() == 0, "Donor chunk size too big; increase number of "
          "processes or read granularity.");
        core::ThrowError(ChunkSizeError);
      }
        
      Chunk.Begin = LocalBegin;
      Chunk.End = LocalEnd;
      Chunk.StartingConnectionID = StartingConnectionID + LocalBegin;

      MPI_File HOFile, XFile;
      MPI_Status Status;
      int MPIError;
      MPI_Offset DatasetOffset, ReadOffset;
      int ReadSize;

      Profiler.StartSync(IMPORT_READ_MPI_IO_OPEN_TIME, ChunkComm);
      // MPI_File_open missing const qualifier for path string on some platforms
      MPIError = MPI_File_open(ChunkComm, const_cast<char *>(HOPath.c_str()), MPI_MODE_RDONLY,
        MPIInfo, &HOFile);
      Profiler.Stop(IMPORT_READ_MPI_IO_OPEN_TIME);
      if (MPIError != MPI_SUCCESS) {
        Logger.LogError(true, "Unable to open file '%s'.", HOPath);
        throw file_open_error();
      }
      auto CloseHO = core::OnScopeExit([&]() {
        Profiler.StartSync(IMPORT_READ_MPI_IO_CLOSE_TIME, ChunkComm);
        MPI_File_close(&HOFile);
        Profiler.Stop(IMPORT_READ_MPI_IO_CLOSE_TIME);
      });

      Profiler.StartSync(IMPORT_READ_MPI_IO_OPEN_TIME, ChunkComm);
      // MPI_File_open missing const qualifier for path string on some platforms
      MPIError = MPI_File_open(ChunkComm, const_cast<char *>(XPath.c_str()), MPI_MODE_RDONLY,
        MPIInfo, &XFile);
      Profiler.Stop(IMPORT_READ_MPI_IO_OPEN_TIME);
      if (MPIError != MPI_SUCCESS) {
        Logger.LogError(true, "Unable to open file '%s'.", XPath);
        throw file_open_error();
      }
      auto CloseX = core::OnScopeExit([&]() {
        Profiler.StartSync(IMPORT_READ_MPI_IO_CLOSE_TIME, ChunkComm);
        MPI_File_close(&XFile);
        Profiler.Stop(IMPORT_READ_MPI_IO_CLOSE_TIME);
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
  //       Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i donor sizes from "
  //         "XINTOUT file '%s'.", GridID, XPath);
  //       OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
  //     }

      // The above version causes "Conditional jump or move depends on uninitialized value(s)"
      // errors in valgrind with MPICH 3.2... not sure why
      DatasetOffset = XSizesOffset;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        ReadOffset = DatasetOffset + LocalBegin*sizeof(int);
        Profiler.StartSync(IMPORT_READ_MPI_IO_OTHER_TIME, ChunkComm);
        MPI_File_set_view(XFile, ReadOffset, MPI_INT, MPI_INT, "native", MPIInfo);
        Profiler.Stop(IMPORT_READ_MPI_IO_OTHER_TIME);
        File_read_all_endian(XFile, Sizes.Data(iDim,0), int(NumLocalDonors), MPI_INT, Endian,
          &Status, Profiler, ChunkComm);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        error ReadSizesError = error::NONE;
        if (ReadSize < NumLocalDonors) ReadSizesError = error::FILE_READ;
        core::SyncError(ReadSizesError, ChunkComm);
        if (ReadSizesError != error::NONE) {
          Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i donor sizes from XINTOUT "
            "file '%s'.", GridID, XPath);
          core::ThrowError(ReadSizesError);
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
  //       Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i donor cells from "
  //         XINTOUT file '%s'.", GridID, HOPath);
  //       OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
  //     }

      // The above version causes "Conditional jump or move depends on uninitialized value(s)"
      // errors in valgrind with MPICH 3.2... not sure why
      DatasetOffset = HOCellsOffset;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        ReadOffset = DatasetOffset + LocalBegin*sizeof(int);
        Profiler.StartSync(IMPORT_READ_MPI_IO_OTHER_TIME, ChunkComm);
        MPI_File_set_view(HOFile, ReadOffset, MPI_INT, MPI_INT, "native", MPIInfo);
        Profiler.Stop(IMPORT_READ_MPI_IO_OTHER_TIME);
        File_read_all_endian(HOFile, Data.Extents.Data(0,iDim,0), int(NumLocalDonors), MPI_INT,
          Endian, &Status, Profiler, ChunkComm);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        error ReadCellsError = error::NONE;
        if (ReadSize < NumLocalDonors) ReadCellsError = error::FILE_READ;
        core::SyncError(ReadCellsError, ChunkComm);
        if (ReadCellsError != error::NONE) {
          Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i donor cells from XINTOUT "
            "file '%s'.", GridID, HOPath);
          core::ThrowError(ReadCellsError);
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
  //       Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i donor Coords from "
  //         "XINTOUT file '%s'.", GridID, HOPath);
  //       OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
  //     }

      // The above version causes "Conditional jump or move depends on uninitialized value(s)"
      // errors in valgrind with MPICH 3.2... not sure why
      DatasetOffset = HOCoordsOffset;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        ReadOffset = DatasetOffset + LocalBegin*sizeof(double);
        Profiler.StartSync(IMPORT_READ_MPI_IO_OTHER_TIME, ChunkComm);
        MPI_File_set_view(HOFile, ReadOffset, MPI_DOUBLE, MPI_DOUBLE, "native", MPIInfo);
        Profiler.Stop(IMPORT_READ_MPI_IO_OTHER_TIME);
        File_read_all_endian(HOFile, Data.Coords.Data(iDim,0), int(NumLocalDonors), MPI_DOUBLE,
          Endian, &Status, Profiler, ChunkComm);
        MPI_Get_count(&Status, MPI_DOUBLE, &ReadSize);
        error ReadCoordsError = error::NONE;
        if (ReadSize < NumLocalDonors) ReadCoordsError = error::FILE_READ;
        core::SyncError(ReadCoordsError, ChunkComm);
        if (ReadCoordsError != error::NONE) {
          Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i donor Coords from XINTOUT "
            "file '%s'.", GridID, HOPath);
          core::ThrowError(ReadCoordsError);
        }
        DatasetOffset += NumDonors*sizeof(double);
      }

      tuple<long long> NumLocalInterpCoefs = {0,0,0};
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        for (long long iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
          int Size = Data.Extents(1,iDim,iDonor) - Data.Extents(0,iDim,iDonor);
          NumLocalInterpCoefs[iDim] += Size;
        }
      }

      ChunkSizeError = error::NONE;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        if (NumLocalInterpCoefs[iDim]*sizeof(double) > std::numeric_limits<int>::max()) {
          ChunkSizeError = error::FILE_READ;
          break;
        }
      }
      core::SyncError(ChunkSizeError, ChunkComm);
      if (ChunkSizeError != error::NONE) {
        Logger.LogError(ChunkComm.Rank() == 0, "Donor chunk size too big; increase number of "
          "processes or read granularity.");
        core::ThrowError(ChunkSizeError);
      }

      tuple<long long> NumInterpCoefs;
      tuple<long long> NumInterpCoefsBeforeChunk;
      MPI_Allreduce(&NumLocalInterpCoefs, &NumInterpCoefs, MAX_DIMS, MPI_LONG_LONG, MPI_SUM,
        ChunkComm);
      MPI_Scan(&NumLocalInterpCoefs, &NumInterpCoefsBeforeChunk, MAX_DIMS, MPI_LONG_LONG, MPI_SUM,
        ChunkComm);
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        NumInterpCoefsBeforeChunk[iDim] -= NumLocalInterpCoefs[iDim];
      }

      array<double> InterpCoefs({NumLocalInterpCoefs[0]+NumLocalInterpCoefs[1]+
        NumLocalInterpCoefs[2]});
      double *Buffer = InterpCoefs.Data();
      DatasetOffset = XInterpCoefsOffset;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        ReadOffset = DatasetOffset + NumInterpCoefsBeforeChunk[iDim]*sizeof(double);
        Profiler.StartSync(IMPORT_READ_MPI_IO_OTHER_TIME, ChunkComm);
        MPI_File_set_view(XFile, ReadOffset, MPI_DOUBLE, MPI_DOUBLE, "native", MPIInfo);
        Profiler.Stop(IMPORT_READ_MPI_IO_OTHER_TIME);
        File_read_all_endian(XFile, Buffer, int(NumLocalInterpCoefs[iDim]), MPI_DOUBLE, Endian,
          &Status, Profiler, ChunkComm);
        MPI_Get_count(&Status, MPI_DOUBLE, &ReadSize);
        error ReadInterpCoefsError = error::NONE;
        if (ReadSize < NumLocalInterpCoefs[iDim]) ReadInterpCoefsError = error::FILE_READ;
        core::SyncError(ReadInterpCoefsError, ChunkComm);
        if (ReadInterpCoefsError != error::NONE) {
          Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i donor interpolation "
            "coefficients from XINTOUT file '%s'.", GridID, XPath);
          core::ThrowError(ReadInterpCoefsError);
        }
        Buffer += NumLocalInterpCoefs[iDim];
        DatasetOffset += NumInterpCoefs[iDim]*sizeof(double);
      }

      tuple<long long> iNextCoef;
      iNextCoef[0] = 0;
      iNextCoef[1] = iNextCoef[0] + NumLocalInterpCoefs[0];
      iNextCoef[2] = iNextCoef[1] + NumLocalInterpCoefs[1];
      for (long long iDonor = 0; iDonor < NumLocalDonors; ++iDonor) {
        for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
          int Size = Data.Extents(1,iDim,iDonor) - Data.Extents(0,iDim,iDonor);
          for (int iPoint = 0; iPoint < Size; ++iPoint) {
            Data.InterpCoefs(iDim,iPoint,iDonor) = InterpCoefs(iNextCoef[iDim]);
            ++iNextCoef[iDim];
          }
        }
      }

    }

  } catch (const exception &Exception) {

    ReadDonorsError = Exception.Error();

  }

  core::SyncError(ReadDonorsError, Comm);
  core::CheckError(ReadDonorsError);

}

void ReadReceivers(xintout_grid &XINTOUTGrid, const std::string &HOPath, long long NumReceivers,
  MPI_Offset HOPointsOffset, MPI_Offset HOConnectionIDsOffset, endian Endian, xintout_format Format,
  int ReadGranularityAdjust, MPI_Info MPIInfo) {

  int GridID = XINTOUTGrid.ID;

  core::comm_view Comm = XINTOUTGrid.Comm;

  core::logger &Logger = XINTOUTGrid.Context->core_Logger();
  core::profiler &Profiler = XINTOUTGrid.Context->core_Profiler();

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

  if (Logger.LoggingDebug()) {
    std::string NumReceiversString = core::FormatNumber(NumReceivers, "receivers", "receiver");
    std::string NumRanksString = core::FormatNumber(NumChunks, "I/O ranks", "I/O rank");
    Logger.LogDebug(Comm.Rank() == 0, 0, "Grid %s has %s; using %s.", XINTOUTGrid.Name,
      NumReceiversString, NumRanksString);
  }

  XINTOUTReceivers.ChunkSize = ChunkSize;
  XINTOUTReceivers.HasChunk = HasChunk;

  core::comm ChunkComm = core::CreateSubsetComm(Comm, HasChunk);

  error ReadReceiversError = error::NONE;

  try {

    if (HasChunk) {

      xintout_receiver_chunk &Chunk = XINTOUTReceivers.Chunk;

      long long LocalBegin = ChunkSize*ChunkComm.Rank();
      long long LocalEnd = std::min(ChunkSize*(ChunkComm.Rank()+1), NumReceivers);
      long long NumLocalReceivers = LocalEnd - LocalBegin;

      error ChunkSizeError = error::NONE;
      if (NumLocalReceivers*sizeof(long long) > std::numeric_limits<int>::max()) {
        ChunkSizeError = error::FILE_READ;
      }
      core::SyncError(ChunkSizeError, ChunkComm);
      if (ChunkSizeError != error::NONE) {
        Logger.LogError(ChunkComm.Rank() == 0, "Receiver chunk size too big; increase number of "
          "processes or read granularity.");
        core::ThrowError(ChunkSizeError);
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

      Profiler.StartSync(IMPORT_READ_MPI_IO_OPEN_TIME, ChunkComm);
      // MPI_File_open missing const qualifier for path string on some platforms
      MPIError = MPI_File_open(ChunkComm, const_cast<char *>(HOPath.c_str()), MPI_MODE_RDONLY,
        MPIInfo, &HOFile);
      Profiler.Stop(IMPORT_READ_MPI_IO_OPEN_TIME);
      if (MPIError != MPI_SUCCESS) {
        Logger.LogError(true, "Unable to open file '%s'.", HOPath);
        throw file_open_error();
      }
      auto CloseHO = core::OnScopeExit([&]() {
        Profiler.StartSync(IMPORT_READ_MPI_IO_CLOSE_TIME, ChunkComm);
        MPI_File_close(&HOFile);
        Profiler.Stop(IMPORT_READ_MPI_IO_CLOSE_TIME);
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
  //       Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i receiver points from "
  //         "XINTOUT file '%s'.", GridID, HOPath);
  //       OVK_EH_CHECK_SKIP_TO(ErrorHandler, Error, done_reading);
  //     }

      // The above version causes "Conditional jump or move depends on uninitialized value(s)"
      // errors in valgrind with MPICH 3.2... not sure why
      DatasetOffset = HOPointsOffset;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        ReadOffset = DatasetOffset + LocalBegin*sizeof(int);
        Profiler.StartSync(IMPORT_READ_MPI_IO_OTHER_TIME, ChunkComm);
        MPI_File_set_view(HOFile, ReadOffset, MPI_INT, MPI_INT, "native", MPIInfo);
        Profiler.Stop(IMPORT_READ_MPI_IO_OTHER_TIME);
        File_read_all_endian(HOFile, Data.Points.Data(iDim,0), int(NumLocalReceivers), MPI_INT,
          Endian, &Status, Profiler, ChunkComm);
        MPI_Get_count(&Status, MPI_INT, &ReadSize);
        error ReadPointsError = error::NONE;
        if (ReadSize < NumLocalReceivers) ReadPointsError = error::FILE_READ;
        core::SyncError(ReadPointsError, ChunkComm);
        if (ReadPointsError != error::NONE) {
          Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i receiver points from "
            "XINTOUT file '%s'.", GridID, HOPath);
          core::ThrowError(ReadPointsError);
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

      array<unsigned char> ConnectionIDsBytes({NumLocalReceivers*ConnectionIDSize});

      DatasetOffset = HOConnectionIDsOffset;
      ReadOffset = DatasetOffset+LocalBegin*ConnectionIDSize;
      Profiler.StartSync(IMPORT_READ_MPI_IO_OTHER_TIME, ChunkComm);
      MPI_File_set_view(HOFile, ReadOffset, ConnectionIDType, ConnectionIDType, "native", MPIInfo);
      Profiler.Stop(IMPORT_READ_MPI_IO_OTHER_TIME);
      File_read_all_endian(HOFile, ConnectionIDsBytes.Data(), int(NumLocalReceivers),
        ConnectionIDType, Endian, &Status, Profiler, ChunkComm);
      MPI_Get_count(&Status, ConnectionIDType, &ReadSize);
      error ReadConnectionIDsError = error::NONE;
      if (ReadSize < NumLocalReceivers) ReadConnectionIDsError = error::FILE_READ;
      core::SyncError(ReadConnectionIDsError, ChunkComm);
      if (ReadConnectionIDsError != error::NONE) {
        Logger.LogError(ChunkComm.Rank() == 0, "Unable to read grid %i receiver connection IDs "
          "from XINTOUT file '%s'.", GridID, HOPath);
        core::ThrowError(ReadConnectionIDsError);
      }

      if (Format == xintout_format::STANDARD) {
        array_view<const int> ConnectionIDs(reinterpret_cast<int *>(ConnectionIDsBytes.Data()),
          {NumLocalReceivers});
        for (long long iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
          // Convert to zero-based indexing
          Chunk.ConnectionIDs(iReceiver) = ConnectionIDs(iReceiver)-1;
        }
      } else {
        array_view<const long long> ConnectionIDs(reinterpret_cast<long long *>(
          ConnectionIDsBytes.Data()), {NumLocalReceivers});
        for (long long iReceiver = 0; iReceiver < NumLocalReceivers; ++iReceiver) {
          // Convert to zero-based indexing
          Chunk.ConnectionIDs(iReceiver) = ConnectionIDs(iReceiver)-1;
        }
      }

    }

  } catch (const exception &Exception) {

    ReadReceiversError = Exception.Error();

  }

  core::SyncError(ReadReceiversError, Comm);
  core::CheckError(ReadReceiversError);

}

void MatchDonorsAndReceivers(xintout &XINTOUT) {

  core::comm_view Comm = XINTOUT.Comm;

  MPI_Barrier(Comm);

  core::logger &Logger = XINTOUT.Context->core_Logger();
  core::profiler &Profiler = XINTOUT.Context->core_Profiler();

  Logger.LogStatus(Comm.Rank() == 0, 0, "Matching donors and receivers...");

  int NumGrids = XINTOUT.NumGrids;
  int NumLocalGrids = XINTOUT.NumLocalGrids;

  array<long long> NumPointsOnGrid({NumGrids}, 0);
  array<long long> NumReceiversOnGrid({NumGrids}, 0);

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    xintout_grid &XINTOUTGrid = XINTOUT.Grids(iLocalGrid);
    int GridID = XINTOUTGrid.ID;
    int iGrid = GridID-1;
    NumPointsOnGrid(iGrid) = 1;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      NumPointsOnGrid(iGrid) *= XINTOUTGrid.GlobalSize[iDim];
    }
    NumReceiversOnGrid(iGrid) = XINTOUTGrid.Receivers.Count;
  }

  MPI_Allreduce(MPI_IN_PLACE, NumPointsOnGrid.Data(), NumGrids, MPI_LONG_LONG, MPI_MAX, Comm);
  MPI_Allreduce(MPI_IN_PLACE, NumReceiversOnGrid.Data(), NumGrids, MPI_LONG_LONG, MPI_MAX, Comm);

  long long NumPoints = 0;
  long long NumConnections = 0;
  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    NumPoints += NumPointsOnGrid(iGrid);
    NumConnections += NumReceiversOnGrid(iGrid);
  }

  NumPointsOnGrid.Clear();
  NumReceiversOnGrid.Clear();

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
    long long Count;
    connection_data Data;
    send_recv():
      Count(0)
    {}
  };

  Profiler.StartSync(IMPORT_MATCH_MAP_TO_BINS_TIME, Comm);

  std::map<int, send_recv> DonorSends, ReceiverSends;

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    xintout_grid &XINTOUTGrid = XINTOUT.Grids(iLocalGrid);
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
    xintout_grid &XINTOUTGrid = XINTOUT.Grids(iLocalGrid);
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

  array<int> DonorSendToRanks;
  for (auto &Pair : DonorSends) {
    int Rank = Pair.first;
    DonorSendToRanks.Append(Rank);
  }

  array<int> ReceiverSendToRanks;
  for (auto &Pair : ReceiverSends) {
    int Rank = Pair.first;
    ReceiverSendToRanks.Append(Rank);
  }

  Profiler.Stop(IMPORT_MATCH_MAP_TO_BINS_TIME);
  Profiler.StartSync(IMPORT_MATCH_HANDSHAKE_TIME, Comm);

  array<int> DonorRecvFromRanks = core::DynamicHandshake(Comm, DonorSendToRanks);
  array<int> ReceiverRecvFromRanks = core::DynamicHandshake(Comm, ReceiverSendToRanks);

  DonorSendToRanks.Clear();
  ReceiverSendToRanks.Clear();

  Profiler.Stop(IMPORT_MATCH_HANDSHAKE_TIME);
  Profiler.StartSync(IMPORT_MATCH_SEND_TO_BINS_TIME, Comm);

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

  array<MPI_Request> Requests;

  // 3 requests per send/recv (1 for each of grid IDs, connection IDs, and points/cells)
  Requests.Reserve(3*(NumDonorSends+NumDonorRecvs+NumReceiverSends+NumReceiverRecvs));

  auto Isend = [&Requests](void *Buffer, long long Count, MPI_Datatype DataType, int DestRank,
    int Tag, MPI_Comm SendComm) {
    OVK_DEBUG_ASSERT(Count <= std::numeric_limits<int>::max(), "Send count too large.");
    MPI_Request &Request = Requests.Append();
    MPI_Isend(Buffer, int(Count), DataType, DestRank, Tag, SendComm, &Request);
  };

  auto Irecv = [&Requests](void *Buffer, long long Count, MPI_Datatype DataType, int SourceRank,
    int Tag, MPI_Comm RecvComm) {
    OVK_DEBUG_ASSERT(Count <= std::numeric_limits<int>::max(), "Receive count too large.");
    MPI_Request &Request = Requests.Append();
    MPI_Irecv(Buffer, int(Count), DataType, SourceRank, Tag, RecvComm, &Request);
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

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

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

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  Profiler.Stop(IMPORT_MATCH_SEND_TO_BINS_TIME);
  Profiler.StartSync(IMPORT_MATCH_RECV_FROM_BINS_TIME, Comm);

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

  Profiler.Stop(IMPORT_MATCH_RECV_FROM_BINS_TIME);
  Profiler.StartSync(IMPORT_MATCH_FILL_CONNECTION_DATA_TIME, Comm);

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

  Profiler.Stop(IMPORT_MATCH_FILL_CONNECTION_DATA_TIME);
  Profiler.StartSync(IMPORT_MATCH_RECV_FROM_BINS_TIME, Comm);

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

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

  Profiler.Stop(IMPORT_MATCH_RECV_FROM_BINS_TIME);
  Profiler.StartSync(IMPORT_MATCH_UNPACK_TIME, Comm);

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
    xintout_grid &XINTOUTGrid = XINTOUT.Grids(iLocalGrid);
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

  Profiler.Stop(IMPORT_MATCH_UNPACK_TIME);

  MPI_Barrier(Comm);

  Logger.LogStatus(Comm.Rank() == 0, 0, "Done matching donors and receivers.");

}

void DistributeConnectivityData(const xintout &XINTOUT, const array<const grid *> &LocalGrids,
  array<donor_data> &LocalDonorData, array<receiver_data> &LocalReceiverData) {

  core::comm_view Comm = XINTOUT.Comm;

  MPI_Barrier(Comm);

  core::logger &Logger = XINTOUT.Context->core_Logger();

  Logger.LogStatus(Comm.Rank() == 0, 0, "Distributing connectivity data to ranks...");

  for (int iLocalGrid = 0; iLocalGrid < XINTOUT.NumLocalGrids; ++iLocalGrid) {
    const xintout_grid &XINTOUTGrid = XINTOUT.Grids(iLocalGrid);
    const grid &Grid = *LocalGrids(iLocalGrid);
    DistributeGridConnectivityData(XINTOUTGrid, Grid, LocalDonorData(iLocalGrid),
      LocalReceiverData(iLocalGrid));
  }

  MPI_Barrier(Comm);

  Logger.LogStatus(Comm.Rank() == 0, 0, "Done distributing connectivity data to ranks.");

}

void DistributeGridConnectivityData(const xintout_grid &XINTOUTGrid, const grid &Grid,
  donor_data &DonorData, receiver_data &ReceiverData) {

  int NumDims = XINTOUTGrid.NumDims;
  core::comm_view Comm = XINTOUTGrid.Comm;
  core::profiler &Profiler = XINTOUTGrid.Context->core_Profiler();

  const xintout_donors &XINTOUTDonors = XINTOUTGrid.Donors;
  const xintout_receivers &XINTOUTReceivers = XINTOUTGrid.Receivers;

  const cart &Cart = Grid.Cart();
  const core::partition_hash &Hash = Grid.core_PartitionHash();

  Profiler.StartSync(IMPORT_DISTRIBUTE_MAP_TO_BINS_TIME, Comm);

  std::map<int, core::partition_hash::bin> Bins;

  long long NumChunkDonors = 0;
  long long NumChunkDonorPoints = 0;
  array<int,2> ChunkDonorPoints;
  array<int> ChunkDonorPointBinIndices;
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
    long long iNextDonorPoint = 0;
    for (long long iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      range DonorRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        DonorRange.Begin(iDim) = DonorChunk.Data.Extents(0,iDim,iDonor);
        DonorRange.End(iDim) = DonorChunk.Data.Extents(1,iDim,iDonor);
      }
      bool AwayFromEdge = Cart.Range().Includes(DonorRange);
      if (AwayFromEdge) {
        for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
          for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
            for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
              tuple<int> Point = {i,j,k};
              for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
                ChunkDonorPoints(iDim,iNextDonorPoint) = Point[iDim];
              }
              ++iNextDonorPoint;
            }
          }
        }
      } else {
        for (int k = DonorRange.Begin(2); k < DonorRange.End(2); ++k) {
          for (int j = DonorRange.Begin(1); j < DonorRange.End(1); ++j) {
            for (int i = DonorRange.Begin(0); i < DonorRange.End(0); ++i) {
              tuple<int> Point = Cart.PeriodicAdjust({i,j,k});
              for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
                ChunkDonorPoints(iDim,iNextDonorPoint) = Point[iDim];
              }
              ++iNextDonorPoint;
            }
          }
        }
      }
    }
    ChunkDonorPointBinIndices.Resize({NumChunkDonorPoints});
    Hash.MapToBins(ChunkDonorPoints, ChunkDonorPointBinIndices);
    for (long long iDonorPoint = 0; iDonorPoint < NumChunkDonorPoints; ++iDonorPoint) {
      int BinIndex = ChunkDonorPointBinIndices(iDonorPoint);
      auto Iter = Bins.lower_bound(BinIndex);
      if (Iter == Bins.end() || Iter->first > BinIndex) {
        Bins.emplace_hint(Iter, BinIndex, core::partition_hash::bin());
      }
    }
  }

  long long NumChunkReceivers = 0;
  array<int> ChunkReceiverBinIndices;
  if (XINTOUTReceivers.HasChunk) {
    const xintout_receiver_chunk &ReceiverChunk = XINTOUTReceivers.Chunk;
    NumChunkReceivers = ReceiverChunk.End - ReceiverChunk.Begin;
    ChunkReceiverBinIndices.Resize({NumChunkReceivers});
    Hash.MapToBins(ReceiverChunk.Data.Points, ChunkReceiverBinIndices);
    for (long long iReceiver = 0; iReceiver < NumChunkReceivers; ++iReceiver) {
      int BinIndex = ChunkReceiverBinIndices(iReceiver);
      auto Iter = Bins.lower_bound(BinIndex);
      if (Iter == Bins.end() || Iter->first > BinIndex) {
        Bins.emplace_hint(Iter, BinIndex, core::partition_hash::bin());
      }
    }
  }

  Profiler.Stop(IMPORT_DISTRIBUTE_MAP_TO_BINS_TIME);
  Profiler.StartSync(IMPORT_DISTRIBUTE_RETRIEVE_BINS_TIME, Comm);

  Hash.RetrieveBins(Bins);

  Profiler.Stop(IMPORT_DISTRIBUTE_RETRIEVE_BINS_TIME);
  Profiler.StartSync(IMPORT_DISTRIBUTE_FIND_RANKS_TIME, Comm);

  array<int> NumChunkDonorRanks;
  array<int *> ChunkDonorRanks;
  array<int> ChunkDonorRanksData;
  if (XINTOUTDonors.HasChunk) {
    const xintout_donor_chunk &DonorChunk = XINTOUTDonors.Chunk;
    NumChunkDonorRanks.Resize({NumChunkDonors});
    ChunkDonorRanks.Resize({NumChunkDonors});
    ChunkDonorRanksData.Resize({NumChunkDonorPoints});
    Hash.FindPartitions(Bins, ChunkDonorPoints, ChunkDonorPointBinIndices, ChunkDonorRanksData);
    int MaxSize = DonorChunk.Data.MaxSize;
    int MaxPointsInCell = 1;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      MaxPointsInCell *= MaxSize;
    }
    id_set<1> UniqueRanks;
    UniqueRanks.Reserve(MaxPointsInCell);
    long long iDonorPoint = 0;
    for (long long iDonor = 0; iDonor < NumChunkDonors; ++iDonor) {
      ChunkDonorRanks(iDonor) = ChunkDonorRanksData.Data(iDonorPoint);
      int NumPointsInCell = 1;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        NumPointsInCell *= DonorChunk.Data.Extents(1,iDim,iDonor) -
          DonorChunk.Data.Extents(0,iDim,iDonor);
      }
      for (int iPointInCell = 0; iPointInCell < NumPointsInCell; ++iPointInCell) {
        UniqueRanks.Insert(ChunkDonorRanks(iDonor)[iPointInCell]);
      }
      NumChunkDonorRanks(iDonor) = UniqueRanks.Count();
      for (int iRank = 0; iRank < NumChunkDonorRanks(iDonor); ++iRank) {
        ChunkDonorRanks(iDonor)[iRank] = UniqueRanks[iRank];
      }
      UniqueRanks.Clear();
      iDonorPoint += NumPointsInCell;
    }
    ChunkDonorPoints.Clear();
    ChunkDonorPointBinIndices.Clear();
  }

  array<int> ChunkReceiverRanks;
  if (XINTOUTReceivers.HasChunk) {
    const xintout_receiver_chunk &ReceiverChunk = XINTOUTReceivers.Chunk;
    ChunkReceiverRanks.Resize({NumChunkReceivers});
    Hash.FindPartitions(Bins, ReceiverChunk.Data.Points, ChunkReceiverBinIndices,
      ChunkReceiverRanks);
    ChunkReceiverBinIndices.Clear();
  }

  Bins.clear();

  struct donor_send_recv {
    long long Count;
    int MaxSize;
    donor_data Data;
    donor_send_recv():
      Count(0),
      MaxSize(0)
    {}
  };

  struct receiver_send_recv {
    long long Count;
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
      for (int iRank = 0; iRank < NumChunkDonorRanks(iDonor); ++iRank) {
        int Rank = ChunkDonorRanks(iDonor)[iRank];
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
      int Rank = ChunkReceiverRanks(iReceiver);
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
      for (int iRank = 0; iRank < NumChunkDonorRanks(iDonor); ++iRank) {
        int Rank = ChunkDonorRanks(iDonor)[iRank];
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
    NumChunkDonorRanks.Clear();
    ChunkDonorRanks.Clear();
    ChunkDonorRanksData.Clear();
  }

  if (XINTOUTReceivers.HasChunk) {
    const xintout_receiver_chunk &ReceiverChunk = XINTOUTReceivers.Chunk;
    for (long long iReceiver = 0; iReceiver < NumChunkReceivers; ++iReceiver) {
      int Rank = ChunkReceiverRanks(iReceiver);
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
    ChunkReceiverRanks.Clear();
  }

  int NumDonorSends = DonorSends.size();
  int NumReceiverSends = ReceiverSends.size();

  array<int> DonorSendToRanks;
  for (auto &Pair : DonorSends) {
    int Rank = Pair.first;
    DonorSendToRanks.Append(Rank);
  }

  array<int> ReceiverSendToRanks;
  for (auto &Pair : ReceiverSends) {
    int Rank = Pair.first;
    ReceiverSendToRanks.Append(Rank);
  }

  Profiler.Stop(IMPORT_DISTRIBUTE_FIND_RANKS_TIME);
  Profiler.StartSync(IMPORT_DISTRIBUTE_HANDSHAKE_TIME, Comm);

  array<int> DonorRecvFromRanks = core::DynamicHandshake(Comm, DonorSendToRanks);
  array<int> ReceiverRecvFromRanks = core::DynamicHandshake(Comm, ReceiverSendToRanks);

  DonorSendToRanks.Clear();
  ReceiverSendToRanks.Clear();

  Profiler.Stop(IMPORT_DISTRIBUTE_HANDSHAKE_TIME);
  Profiler.StartSync(IMPORT_DISTRIBUTE_SEND_DATA_TIME, Comm);

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

  array<MPI_Request> Requests;

  // 5 requests per donor send/recv (1 for each of Extents, Coords, interp coefs, destination grid
  // IDs, and destination points)
  // 3 requests per receiver send/recv (1 for each of points, source grid IDs, and source cells)
  Requests.Reserve(5*(NumDonorSends+NumDonorRecvs)+3*(NumReceiverSends+NumReceiverRecvs));

  auto Isend = [&Requests](void *Buffer, long long Count, MPI_Datatype DataType, int DestRank,
    int Tag, MPI_Comm SendComm) {
    OVK_DEBUG_ASSERT(Count <= std::numeric_limits<int>::max(), "Send count too large.");
    MPI_Request &Request = Requests.Append();
    MPI_Isend(Buffer, int(Count), DataType, DestRank, Tag, SendComm, &Request);
  };

  auto Irecv = [&Requests](void *Buffer, long long Count, MPI_Datatype DataType, int SourceRank,
    int Tag, MPI_Comm RecvComm) {
    OVK_DEBUG_ASSERT(Count <= std::numeric_limits<int>::max(), "Receive count too large.");
    MPI_Request &Request = Requests.Append();
    MPI_Irecv(Buffer, Count, DataType, SourceRank, Tag, RecvComm, &Request);
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

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

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

  MPI_Waitall(Requests.Count(), Requests.Data(), MPI_STATUSES_IGNORE);

  Requests.Clear();

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

  Profiler.Stop(IMPORT_DISTRIBUTE_SEND_DATA_TIME);

}

void ImportConnectivityData(int NumGrids, int NumLocalGrids, const array<int> &LocalGridIDs,
  const array<donor_data> &LocalDonorData, const array<receiver_data> &LocalReceiverData, domain
  &Domain, int ConnectivityComponentID) {

  core::comm_view Comm = Domain.core_Comm();

  MPI_Barrier(Comm);

  core::logger &Logger = Domain.Context().core_Logger();

  Logger.LogStatus(Comm.Rank() == 0, 0, "Importing connectivity data into domain %s...",
    Domain.Name());

  OVK_DEBUG_ASSERT(ConnectivityComponentID >= 0, "Invalid connectivity component ID.");
  OVK_DEBUG_ASSERT(Domain.ComponentExists(ConnectivityComponentID), "Component %i does not "
    "exist.", ConnectivityComponentID);

  array<long long,2> NumConnections({{NumGrids,NumGrids}}, 0);

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {
    int jGrid = LocalGridIDs(iLocalGrid)-1;
    for (long long iReceiver = 0; iReceiver < LocalReceiverData(iLocalGrid).Count; ++iReceiver) {
      int iGrid = LocalReceiverData(iLocalGrid).SourceGridIDs(iReceiver)-1;
      ++NumConnections(iGrid,jGrid);
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, NumConnections.Data(), NumConnections.Count(), MPI_LONG_LONG, MPI_SUM,
    Comm);

  auto ConnectivityComponent = Domain.EditComponent<connectivity_component>(
    ConnectivityComponentID);

  for (int iGrid = 0; iGrid < NumGrids; ++iGrid) {
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections(iGrid,jGrid) > 0) {
        int MGridID = iGrid+1;
        int NGridID = jGrid+1;
        ConnectivityComponent->CreateConnectivity(MGridID, NGridID);
      }
    }
  }

  for (int iLocalGrid = 0; iLocalGrid < NumLocalGrids; ++iLocalGrid) {

    int LocalGridID = LocalGridIDs(iLocalGrid);
    int iGrid = LocalGridID-1;

    const grid &Grid = Domain.Grid(LocalGridID);
    core::comm_view GridComm = Grid.core_Comm();

    int NumConnectivityMs = 0;
    int NumConnectivityNs = 0;
    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      if (NumConnections(iGrid,jGrid) > 0) {
        ++NumConnectivityMs;
      }
      if (NumConnections(jGrid,iGrid) > 0) {
        ++NumConnectivityNs;
      }
    }

    array<edit_handle<connectivity_m>> ConnectivityMs;
    array<edit_handle<connectivity_n>> ConnectivityNs;

    ConnectivityMs.Reserve(NumConnectivityMs);
    ConnectivityNs.Reserve(NumConnectivityNs);

    for (int jGrid = 0; jGrid < NumGrids; ++jGrid) {
      int OtherGridID = jGrid+1;
      if (NumConnections(iGrid,jGrid) > 0) {
        ConnectivityMs.Append(ConnectivityComponent->EditConnectivityM(LocalGridID,
          OtherGridID));
      }
      if (NumConnections(jGrid,iGrid) > 0) {
        ConnectivityNs.Append(ConnectivityComponent->EditConnectivityN(OtherGridID,
          LocalGridID));
      }
    }

    ImportDonors(LocalDonorData(iLocalGrid), GridComm, ConnectivityMs);
    ImportReceivers(LocalReceiverData(iLocalGrid), GridComm, ConnectivityNs);

  }

  ConnectivityComponent.Restore();

  MPI_Barrier(Comm);

  Logger.LogStatus(Comm.Rank() == 0, 0, "Done importing connectivity data into domain %s.",
    Domain.Name());

}

void ImportDonors(const donor_data &GridDonors, core::comm_view Comm, const array<edit_handle<
    connectivity_m>> &ConnectivityMs) {

  int NumDestinationGrids = ConnectivityMs.Count();

  if (NumDestinationGrids == 0) return;

  id_map<1,int> DestinationGridIDToIndex;

  for (int iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    connectivity_m &ConnectivityM = *ConnectivityMs(iDestinationGrid);
    DestinationGridIDToIndex.Insert(ConnectivityM.DestinationGridID(), iDestinationGrid);
  }

  array<long long> Counts({NumDestinationGrids}, 0);
  array<int> MaxSizes({NumDestinationGrids}, 1);

  for (long long iDonor = 0; iDonor < GridDonors.Count; ++iDonor) {
    int DestinationGridID = GridDonors.DestinationGridIDs(iDonor);
    int iDestinationGrid = DestinationGridIDToIndex(DestinationGridID);
    ++Counts(iDestinationGrid);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      int Size = GridDonors.Extents(1,iDim,iDonor) - GridDonors.Extents(0,iDim,iDonor);
      MaxSizes(iDestinationGrid) = std::max(MaxSizes(iDestinationGrid), Size);
    }
  }

  for (int iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    connectivity_m &ConnectivityM = *ConnectivityMs(iDestinationGrid);
    ConnectivityM.Resize(Counts(iDestinationGrid), MaxSizes(iDestinationGrid));
  }

  array<edit_handle<array<int,3>>> AllExtents({NumDestinationGrids});
  array<edit_handle<array<double,2>>> AllCoords({NumDestinationGrids});
  array<edit_handle<array<double,3>>> AllInterpCoefs({NumDestinationGrids});
  array<edit_handle<array<int,2>>> AllDestinations({NumDestinationGrids});

  for (int iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    connectivity_m &ConnectivityM = *ConnectivityMs(iDestinationGrid);
    AllExtents(iDestinationGrid) = ConnectivityM.EditExtents();
    AllCoords(iDestinationGrid) = ConnectivityM.EditCoords();
    AllInterpCoefs(iDestinationGrid) = ConnectivityM.EditInterpCoefs();
    AllDestinations(iDestinationGrid) = ConnectivityM.EditDestinations();
  }

  // Reuse counts for filling in data
  for (int iDestinationGrid = 0; iDestinationGrid < NumDestinationGrids; ++iDestinationGrid) {
    Counts(iDestinationGrid) = 0;
  }

  for (long long iDonor = 0; iDonor < GridDonors.Count; ++iDonor) {
    int DestinationGridID = GridDonors.DestinationGridIDs(iDonor);
    int iDestinationGrid = DestinationGridIDToIndex(DestinationGridID);
    array<int,3> &Extents = *AllExtents(iDestinationGrid);
    array<double,2> &Coords = *AllCoords(iDestinationGrid);
    array<double,3> &InterpCoefs = *AllInterpCoefs(iDestinationGrid);
    array<int,2> &Destinations = *AllDestinations(iDestinationGrid);
    long long &iNext = Counts(iDestinationGrid);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Extents(0,iDim,iNext) = GridDonors.Extents(0,iDim,iDonor);
      Extents(1,iDim,iNext) = GridDonors.Extents(1,iDim,iDonor);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Coords(iDim,iNext) = GridDonors.Coords(iDim,iDonor);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      int Size = GridDonors.Extents(1,iDim,iDonor)-GridDonors.Extents(0,iDim,iDonor);
      for (int iPoint = 0; iPoint < Size; ++iPoint) {
        InterpCoefs(iDim,iPoint,iNext) = GridDonors.InterpCoefs(iDim,iPoint,iDonor);
      }
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Destinations(iDim,iNext) = GridDonors.DestinationPoints(iDim,iDonor);
    }
    ++iNext;
  }

}

void ImportReceivers(const receiver_data &GridReceivers, core::comm_view Comm, const array<
  edit_handle<connectivity_n>> &ConnectivityNs) {

  int NumSourceGrids = ConnectivityNs.Count();

  if (NumSourceGrids == 0) return;

  id_map<1,int> SourceGridIDToIndex;

  for (int iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    connectivity_n &ConnectivityN = *ConnectivityNs(iSourceGrid);
    SourceGridIDToIndex.Insert(ConnectivityN.SourceGridID(), iSourceGrid);
  }

  array<long long> Counts({NumSourceGrids}, 0);

  for (long long iReceiver = 0; iReceiver < GridReceivers.Count; ++iReceiver) {
    int SourceGridID = GridReceivers.SourceGridIDs(iReceiver);
    int iSourceGrid = SourceGridIDToIndex(SourceGridID);
    ++Counts(iSourceGrid);
  }

  for (int iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    connectivity_n &ConnectivityN = *ConnectivityNs(iSourceGrid);
    ConnectivityN.Resize(Counts(iSourceGrid));
  }

  array<edit_handle<array<int,2>>> AllPoints({NumSourceGrids});
  array<edit_handle<array<int,2>>> AllSources({NumSourceGrids});

  for (int iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    connectivity_n &ConnectivityN = *ConnectivityNs(iSourceGrid);
    AllPoints(iSourceGrid) = ConnectivityN.EditPoints();
    AllSources(iSourceGrid) = ConnectivityN.EditSources();
  }

  // Reuse counts for filling in data
  for (int iSourceGrid = 0; iSourceGrid < NumSourceGrids; ++iSourceGrid) {
    Counts(iSourceGrid) = 0;
  }

  for (long long iReceiver = 0; iReceiver < GridReceivers.Count; ++iReceiver) {
    int SourceGridID = GridReceivers.SourceGridIDs(iReceiver);
    int iSourceGrid = SourceGridIDToIndex(SourceGridID);
    array<int,2> &Points = *AllPoints(iSourceGrid);
    array<int,2> &Sources = *AllSources(iSourceGrid);
    long long &iNext = Counts(iSourceGrid);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Points(iDim,iNext) = GridReceivers.Points(iDim,iReceiver);
    }
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Sources(iDim,iNext) = GridReceivers.SourceCells(iDim,iReceiver);
    }
    ++iNext;
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
  endian Endian, MPI_Status *Status, core::profiler &Profiler, core::comm_view Comm) {

  Profiler.StartSync(IMPORT_READ_MPI_IO_READ_TIME, Comm);
  int MPIError = MPI_File_read_all(File, Buffer, Count, DataType, Status);
  Profiler.Stop(IMPORT_READ_MPI_IO_READ_TIME);

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

  Profiler.Start(IMPORT_READ_MPI_IO_READ_TIME);
  int MPIError = MPI_File_read_at(File, Offset, Buffer, Count, DataType, Status);
  Profiler.Stop(IMPORT_READ_MPI_IO_READ_TIME);

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

  unsigned char EndianTest[] = {1,0};

  if(*reinterpret_cast<short *>(&EndianTest[0]) == 1) {
    return endian::LITTLE;
  } else {
    return endian::BIG;
  }

}

void SwapEndian(void *Data, int ElementSize, int NumElements) {

  array_view<unsigned char> DataBytes(static_cast<unsigned char *>(Data),
    {NumElements*ElementSize});

  array<unsigned char> Element({ElementSize});

  for (int i = 0; i < NumElements; ++i) {
    for (int j = 0; j < ElementSize; ++j) {
      Element(j) = DataBytes(ElementSize*i+j);
    }
    for (int j = 0; j < ElementSize; ++j) {
      DataBytes(ElementSize*i+j) = Element(ElementSize-j-1);
    }
  }

}

}}
