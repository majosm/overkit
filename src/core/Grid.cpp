// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Grid.hpp"

#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/ErrorHandler.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/PartitionHash.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace ovk {

namespace {
void CreateNeighbors(grid &Grid);
void DestroyNeighbors(grid &Grid);
void PrintGridSummary(const grid &Grid);
void PrintGridDecomposition(const grid &Grid);
}

namespace core {

void CreateGrid(grid &Grid, int ID, const grid_params &Params, core::logger &Logger,
  core::error_handler &ErrorHandler) {

  int NumDims = Params.NumDims_;

  Grid.Comm_ = core::comm(Params.Comm_);

  MPI_Barrier(Grid.Comm_);

  Grid.Logger_ = &Logger;
  Grid.ErrorHandler_ = &ErrorHandler;

  Grid.ID_ = ID;

  if (Params.Name_.length() > 0) {
    Grid.Name_ = Params.Name_;
  } else {
    Grid.Name_ = "Grid" + std::to_string(ID);
  }

  Grid.NumDims_ = NumDims;

  DefaultCart(Grid.Cart_, NumDims);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    if (Params.Periodic_[iDim] && Params.PeriodicStorage_ == periodic_storage::DUPLICATED) {
      Grid.Cart_.Size[iDim] = Params.Size_[iDim]-1;
    } else {
      Grid.Cart_.Size[iDim] = Params.Size_[iDim];
    }
    Grid.Cart_.Periodic[iDim] = Params.Periodic_[iDim];
  }

  Grid.PeriodicStorage_ = Params.PeriodicStorage_;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Grid.PeriodicLength_[iDim] = Params.PeriodicLength_[iDim];
  }

  Grid.GeometryType_ = Params.GeometryType_;

  DefaultRange(Grid.GlobalRange_, NumDims);
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Grid.GlobalRange_.Begin[iDim] = 0;
    Grid.GlobalRange_.End[iDim] = Params.Size_[iDim];
  }

  Grid.LocalRange_ = Params.LocalRange_;

  core::CreatePartitionHash(Grid.PartitionHash_, Grid.NumDims_, Grid.Comm_, Grid.GlobalRange_,
    Grid.LocalRange_);

  CreateNeighbors(Grid);

  if (Grid.Comm_.Rank() == 0) {
    PrintGridSummary(Grid);
  }

  if (OVK_DEBUG) {
    PrintGridDecomposition(Grid);
  }

  MPI_Barrier(Grid.Comm_);

}

void DestroyGrid(grid &Grid) {

  MPI_Barrier(Grid.Comm_);

  core::DestroyPartitionHash(Grid.PartitionHash_);

  DestroyNeighbors(Grid);

  MPI_Barrier(Grid.Comm_);

  LogStatus(*Grid.Logger_, Grid.Comm_.Rank() == 0, 0, "Destroyed grid %s.", Grid.Name_);

  Grid.Comm_.Reset();

}

}

void GetGridID(const grid &Grid, int &ID) {

  ID = Grid.ID_;

}

void GetGridName(const grid &Grid, std::string &Name) {

  Name = Grid.Name_;

}

void GetGridDimension(const grid &Grid, int &NumDims) {

  NumDims = Grid.NumDims_;

}

void GetGridComm(const grid &Grid, MPI_Comm &Comm) {

  Comm = Grid.Comm_.Get();

}

void GetGridCommSize(const grid &Grid, int &CommSize) {

  CommSize = Grid.Comm_.Size();

}

void GetGridCommRank(const grid &Grid, int &CommRank) {

  CommRank = Grid.Comm_.Rank();

}

void GetGridCart(const grid &Grid, cart &Cart) {

  Cart = Grid.Cart_;

}

void GetGridSize(const grid &Grid, int *Size) {

  RangeSize(Grid.GlobalRange_, Size);

}

void GetGridPeriodic(const grid &Grid, bool *Periodic) {

  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Periodic[iDim] = Grid.Cart_.Periodic[iDim];
  }

}

void GetGridPeriodicStorage(const grid &Grid, periodic_storage &PeriodicStorage) {

  PeriodicStorage = Grid.PeriodicStorage_;

}

void GetGridPeriodicLength(const grid &Grid, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = Grid.PeriodicLength_[iDim];
  }

}

void GetGridGeometryType(const grid &Grid, geometry_type &GeometryType) {

  GeometryType = Grid.GeometryType_;

}

void GetGridGlobalRange(const grid &Grid, range &GlobalRange) {

  GlobalRange = Grid.GlobalRange_;

}

void GetGridLocalRange(const grid &Grid, range &LocalRange) {

  LocalRange = Grid.LocalRange_;

}

void GetGridGlobalCount(const grid &Grid, long long &NumGlobal) {

  RangeCount(Grid.GlobalRange_, NumGlobal);

}

void GetGridLocalCount(const grid &Grid, long long &NumLocal) {

  RangeCount(Grid.LocalRange_, NumLocal);

}

namespace {

void CreateNeighbors(grid &Grid) {

  int NumDims = Grid.NumDims_;

  range GlobalRange = Grid.GlobalRange_;
  range LocalRange = Grid.LocalRange_;
  int Periodic[MAX_DIMS] = {
    Grid.Cart_.Periodic[0],
    Grid.Cart_.Periodic[1],
    Grid.Cart_.Periodic[2]
  };

  bool HasNeighborsBefore[MAX_DIMS];
  bool HasNeighborsAfter[MAX_DIMS];
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    bool HasPartitionBefore = LocalRange.Begin[iDim] > GlobalRange.Begin[iDim];
    bool HasPartitionAfter = LocalRange.End[iDim] < GlobalRange.End[iDim];
    if (Periodic[iDim]) {
      HasNeighborsBefore[iDim] = HasPartitionBefore || HasPartitionAfter;
      HasNeighborsAfter[iDim] = HasPartitionBefore || HasPartitionAfter;
    } else {
      HasNeighborsBefore[iDim] = HasPartitionBefore;
      HasNeighborsAfter[iDim] = HasPartitionAfter;
    }
  }

  int NumNeighborFaces = 0;
  range NeighborFaces[2*MAX_DIMS];
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    if (HasNeighborsBefore[iDim]) {
      range &Face = NeighborFaces[NumNeighborFaces];
      Face = LocalRange;
      Face.Begin[iDim] -= 1;
      Face.End[iDim] = Face.Begin[iDim]+1;
      for (int jDim = iDim+1; jDim < NumDims; ++jDim) {
        Face.Begin[jDim] -= int(HasNeighborsBefore[jDim]);
        Face.End[jDim] += int(HasNeighborsAfter[jDim]);
      }
      ++NumNeighborFaces;
    }
    if (HasNeighborsAfter[iDim]) {
      range &Face = NeighborFaces[NumNeighborFaces];
      Face = LocalRange;
      Face.End[iDim] += 1;
      Face.Begin[iDim] = Face.End[iDim]-1;
      for (int jDim = iDim+1; jDim < NumDims; ++jDim) {
        Face.Begin[jDim] -= int(HasNeighborsBefore[jDim]);
        Face.End[jDim] += int(HasNeighborsAfter[jDim]);
      }
      ++NumNeighborFaces;
    }
  }

  long long NumNeighborPoints = 0;
  for (int iFace = 0; iFace < NumNeighborFaces; ++iFace) {
    long long NumPoints;
    RangeCount(NeighborFaces[iFace], NumPoints);
    NumNeighborPoints += NumPoints;
  }

  std::vector<int> NeighborPointsFlat(MAX_DIMS*NumNeighborPoints);
  int *NeighborPoints[MAX_DIMS] = {
    NeighborPointsFlat.data(),
    NeighborPointsFlat.data() + NumNeighborPoints,
    NeighborPointsFlat.data() + 2*NumNeighborPoints
  };

  long long iNextPoint = 0;
  for (int iFace = 0; iFace < NumNeighborFaces; ++iFace) {
    range &Face = NeighborFaces[iFace];
    for (int k = Face.Begin[2]; k < Face.End[2]; ++k) {
      for (int j = Face.Begin[1]; j < Face.End[1]; ++j) {
        for (int i = Face.Begin[0]; i < Face.End[0]; ++i) {
          int Point[MAX_DIMS] = {i, j, k};
          CartPeriodicAdjust(Grid.Cart_, Point, Point);
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            NeighborPoints[iDim][iNextPoint] = Point[iDim];
          }
          ++iNextPoint;
        }
      }
    }
  }

  std::vector<int> NeighborPointBinIndices(NumNeighborPoints);

  core::MapToPartitionBins(Grid.PartitionHash_, NumNeighborPoints, NeighborPoints,
    NeighborPointBinIndices.data());

  std::map<int, core::partition_bin> Bins;

  for (long long iPoint = 0; iPoint < NumNeighborPoints; ++iPoint) {
    Bins.emplace(NeighborPointBinIndices[iPoint], core::partition_bin());
  }

  core::RetrievePartitionBins(Grid.PartitionHash_, Bins);

  std::vector<int> NeighborPointRanks(NumNeighborPoints);

  core::FindPartitions(Grid.PartitionHash_, Bins, NumNeighborPoints, NeighborPoints,
    NeighborPointBinIndices.data(), NeighborPointRanks.data());

  Bins.clear();

  std::set<int> NeighborRanks;

  for (long long iPoint = 0; iPoint < NumNeighborPoints; ++iPoint) {
    NeighborRanks.insert(NeighborPointRanks[iPoint]);
  }

  NeighborPointsFlat.clear();
  NeighborPointBinIndices.clear();
  NeighborPointRanks.clear();

  int NumNeighbors = NeighborRanks.size();

  Grid.Neighbors_.resize(NumNeighbors);

  auto RankIter = NeighborRanks.begin();
  for (auto &Neighbor : Grid.Neighbors_) {
    Neighbor.Rank = *RankIter;
    ++RankIter;
  }

  std::vector<int> NeighborRangesFlat(2*MAX_DIMS*NumNeighbors);

  int LocalRangeFlat[2*MAX_DIMS];
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    LocalRangeFlat[iDim] = LocalRange.Begin[iDim];
    LocalRangeFlat[MAX_DIMS+iDim] = LocalRange.End[iDim];
  }

  int *NeighborRangeFlat;

  std::vector<MPI_Request> Requests(2*NumNeighbors);
  NeighborRangeFlat = NeighborRangesFlat.data();
  MPI_Request *Request = Requests.data();
  for (auto &Neighbor : Grid.Neighbors_) {
    MPI_Irecv(NeighborRangeFlat, 2*MAX_DIMS, MPI_INT, Neighbor.Rank, 0, Grid.Comm_, Request);
    MPI_Isend(LocalRangeFlat, 2*MAX_DIMS, MPI_INT, Neighbor.Rank, 0, Grid.Comm_, Request+1);
    NeighborRangeFlat += 2*MAX_DIMS;
    Request += 2;
  }
  MPI_Waitall(2*NumNeighbors, Requests.data(), MPI_STATUSES_IGNORE);

  NeighborRangeFlat = NeighborRangesFlat.data();
  for (auto &Neighbor : Grid.Neighbors_) {
    DefaultRange(Neighbor.LocalRange, NumDims);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Neighbor.LocalRange.Begin[iDim] = NeighborRangeFlat[iDim];
      Neighbor.LocalRange.End[iDim] = NeighborRangeFlat[MAX_DIMS+iDim];
    }
    NeighborRangeFlat += 2*MAX_DIMS;
  }

}

void DestroyNeighbors(grid &Grid) {

  Grid.Neighbors_.clear();

}

void PrintGridSummary(const grid &Grid) {

  std::string GlobalSizeString;
  for (int iDim = 0; iDim < Grid.NumDims_; ++iDim) {
    GlobalSizeString += core::FormatNumber(Grid.GlobalRange_.End[iDim]);
    if (iDim != Grid.NumDims_-1) GlobalSizeString += " x ";
  }

  long long TotalPoints;
  RangeCount(Grid.GlobalRange_, TotalPoints);

  std::string TotalPointsString = core::FormatNumber(TotalPoints, "points", "point");
  std::string ProcessesString = core::FormatNumber(Grid.Comm_.Size(), "processes", "process");

  LogStatus(*Grid.Logger_, true, 0, "Created grid %s (ID=%i): %s (%s) on %s.", Grid.Name_, Grid.ID_,
    GlobalSizeString, TotalPointsString, ProcessesString);

}

void PrintGridDecomposition(const grid &Grid) {

  const char *DimNames[3] = {"i", "j", "k"};
  std::string LocalRangeString;
  for (int iDim = 0; iDim < Grid.NumDims_; ++iDim) {
    std::string LocalBeginString = core::FormatNumber(Grid.LocalRange_.Begin[iDim]);
    std::string LocalEndString = core::FormatNumber(Grid.LocalRange_.End[iDim]);
    LocalRangeString += std::string(DimNames[iDim]) + '=' + LocalBeginString + ':' + LocalEndString;
    if (iDim != Grid.NumDims_-1) LocalRangeString += ", ";
  }

  long long TotalLocalPoints;
  RangeCount(Grid.LocalRange_, TotalLocalPoints);
  std::string TotalLocalPointsString = core::FormatNumber(TotalLocalPoints, "points", "point");

  std::string NeighborRanksString;
  for (int iNeighbor = 0; iNeighbor < int(Grid.Neighbors_.size()); ++iNeighbor) {
    // List separated by commas, so don't add thousands separators
    NeighborRanksString += std::to_string(Grid.Neighbors_[iNeighbor].Rank);
    if (iNeighbor != int(Grid.Neighbors_.size())-1) NeighborRanksString += ", ";
  }

  LogStatus(*Grid.Logger_, Grid.Comm_.Rank() == 0, 0, "Grid %s decomposition info:", Grid.Name_);

  MPI_Barrier(Grid.Comm_);

  for (int OtherRank = 0; OtherRank < Grid.Comm_.Size(); ++OtherRank) {
    if (OtherRank == Grid.Comm_.Rank()) {
      std::string RankString = core::FormatNumber(Grid.Comm_.Rank());
      LogStatus(*Grid.Logger_, true, 1, "Rank %s (global rank @rank@) contains %s (%s).",
        RankString, LocalRangeString, TotalLocalPointsString);
      if (Grid.Neighbors_.size() > 0) {
        LogStatus(*Grid.Logger_, true, 1, "Rank %s has neighbors: %s", RankString,
          NeighborRanksString);
      }
    }
    MPI_Barrier(Grid.Comm_);
  }

}

}

namespace core {

const comm &GetGridComm(const grid &Grid) {

  return Grid.Comm_;

}

const std::vector<grid_neighbor> &GetGridNeighbors(const grid &Grid) {

  return Grid.Neighbors_;

}

const partition_hash &GetGridPartitionHash(const grid &Grid) {

  return Grid.PartitionHash_;

}

}

void CreateGridParams(grid_params &Params, int NumDims) {

  OVK_DEBUG_ASSERT(NumDims == 2 || NumDims == 3, "Invalid dimension.");

  Params.NumDims_ = NumDims;
  Params.Comm_ = MPI_COMM_NULL;

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    Params.Size_[iDim] = 0;
  }
  for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
    Params.Size_[iDim] = 1;
  }

  Params.Periodic_[0] = false;
  Params.Periodic_[1] = false;
  Params.Periodic_[2] = false;
  Params.PeriodicStorage_ = periodic_storage::UNIQUE;
  Params.PeriodicLength_[0] = 0.;
  Params.PeriodicLength_[1] = 0.;
  Params.PeriodicLength_[2] = 0.;
  Params.GeometryType_ = geometry_type::CURVILINEAR;

  DefaultRange(Params.LocalRange_, NumDims);

}

void DestroyGridParams(grid_params &) {}

void GetGridParamName(const grid_params &Params, std::string &Name) {

  Name = Params.Name_;

}

void SetGridParamName(grid_params &Params, std::string Name) {

  Params.Name_ = std::move(Name);

}

void GetGridParamDimension(const grid_params &Params, int &NumDims) {

  NumDims = Params.NumDims_;

}

void GetGridParamComm(const grid_params &Params, MPI_Comm &Comm) {

  Comm = Params.Comm_;

}

void SetGridParamComm(grid_params &Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Params.Comm_ = Comm;

}

void GetGridParamSize(const grid_params &Params, int *Size) {

  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Size[iDim] = Params.Size_[iDim];
  }

}

void SetGridParamSize(grid_params &Params, const int *Size) {

  OVK_DEBUG_ASSERT(Size, "Invalid size pointer.");

  for (int iDim = 0; iDim < Params.NumDims_; ++iDim) {
    OVK_DEBUG_ASSERT(Size[iDim] > 0, "Size must be greater than 0 in each dimension.");
    Params.Size_[iDim] = Size[iDim];
  }

}

void GetGridParamPeriodic(const grid_params &Params, bool *Periodic) {

  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Periodic[iDim] = Params.Periodic_[iDim];
  }

}

void SetGridParamPeriodic(grid_params &Params, const bool *Periodic) {

  OVK_DEBUG_ASSERT(Periodic, "Invalid periodic pointer.");

  for (int iDim = 0; iDim < Params.NumDims_; ++iDim) {
    Params.Periodic_[iDim] = Periodic[iDim];
  }

}

void GetGridParamPeriodicStorage(const grid_params &Params, periodic_storage &PeriodicStorage) {

  PeriodicStorage = Params.PeriodicStorage_;

}

void SetGridParamPeriodicStorage(grid_params &Params, periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(ValidPeriodicStorage(PeriodicStorage), "Invalid periodic storage.");

  OVK_DEBUG_ASSERT(PeriodicStorage == periodic_storage::UNIQUE, "Duplicated periodic storage "
    "is not currently supported.");

  Params.PeriodicStorage_ = PeriodicStorage;

}

void GetGridParamPeriodicLength(const grid_params &Params, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = Params.PeriodicLength_[iDim];
  }

}

void SetGridParamPeriodicLength(grid_params &Params, const double *PeriodicLength) {

  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  for (int iDim = 0; iDim < Params.NumDims_; ++iDim) {
    OVK_DEBUG_ASSERT(PeriodicLength[iDim] >= 0., "Periodic length must be nonnegative.");
    Params.PeriodicLength_[iDim] = PeriodicLength[iDim];
  }

}

void GetGridParamGeometryType(const grid_params &Params, geometry_type &GeometryType) {

  GeometryType = Params.GeometryType_;

}

void SetGridParamGeometryType(grid_params &Params, geometry_type GeometryType) {

  OVK_DEBUG_ASSERT(ValidGeometryType(GeometryType), "Invalid geometry type.");

  Params.GeometryType_ = GeometryType;

}

void GetGridParamLocalRange(const grid_params &Params, range &LocalRange) {

  LocalRange = Params.LocalRange_;

}

void SetGridParamLocalRange(grid_params &Params, const range &LocalRange) {

  Params.LocalRange_ = LocalRange;

}

namespace core {

void CreateGridInfo(grid_info &Info, const grid *Grid, const comm &Comm) {

  bool IsLocal = Grid != nullptr;
  bool IsRoot = false;
  if (IsLocal) {
    IsRoot = Grid->Comm_.Rank() == 0;
  }

  int RootRank;
  if (IsRoot) RootRank = Comm.Rank();
  core::BroadcastAnySource(&RootRank, 1, MPI_INT, IsRoot, Comm);

  if (IsRoot) {
    Info.ID_ = Grid->ID_;
    Info.NumDims_ = Grid->NumDims_;
  }
  MPI_Bcast(&Info.ID_, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info.NumDims_, 1, MPI_INT, RootRank, Comm);

  int NameLength;
  if (IsRoot) NameLength = Grid->Name_.length();
  MPI_Bcast(&NameLength, 1, MPI_INT, RootRank, Comm);
  std::vector<char> NameChars(NameLength);
  if (IsRoot) NameChars.assign(Grid->Name_.begin(), Grid->Name_.end());
  MPI_Bcast(NameChars.data(), NameLength, MPI_CHAR, RootRank, Comm);
  Info.Name_.assign(NameChars.begin(), NameChars.end());

  Info.RootRank_ = RootRank;

  int PeriodicInt[MAX_DIMS];
  if (IsRoot) {
    Info.Cart_ = Grid->Cart_;
    PeriodicInt[0] = int(Info.Cart_.Periodic[0]);
    PeriodicInt[1] = int(Info.Cart_.Periodic[1]);
    PeriodicInt[2] = int(Info.Cart_.Periodic[2]);
  }
  MPI_Bcast(&Info.Cart_.NumDims, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info.Cart_.Size, MAX_DIMS, MPI_INT, RootRank, Comm);
  MPI_Bcast(&PeriodicInt, MAX_DIMS, MPI_INT, RootRank, Comm);
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Info.Cart_.Periodic[iDim] = bool(PeriodicInt[iDim]);
  }

  if (IsRoot) {
    Info.PeriodicLength_[0] = Grid->PeriodicLength_[0];
    Info.PeriodicLength_[1] = Grid->PeriodicLength_[1];
    Info.PeriodicLength_[2] = Grid->PeriodicLength_[2];
  }
  MPI_Bcast(&Info.PeriodicLength_, MAX_DIMS, MPI_DOUBLE, RootRank, Comm);

  int GeometryTypeInt;
  if (IsRoot) {
    GeometryTypeInt = int(Grid->GeometryType_);
  }
  MPI_Bcast(&GeometryTypeInt, 1, MPI_INT, RootRank, Comm);
  Info.GeometryType_ = ovk::geometry_type(GeometryTypeInt);

  if (IsRoot) {
    Info.GlobalRange_ = Grid->GlobalRange_;
  }
  MPI_Bcast(&Info.GlobalRange_.NumDims, 1, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info.GlobalRange_.Begin, MAX_DIMS, MPI_INT, RootRank, Comm);
  MPI_Bcast(&Info.GlobalRange_.End, MAX_DIMS, MPI_INT, RootRank, Comm);

  Info.IsLocal_ = IsLocal;

}

void DestroyGridInfo(grid_info &Info) {

  Info.Name_.clear();

}

}

void GetGridInfoID(const grid_info &Info, int &ID) {

  ID = Info.ID_;

}

void GetGridInfoName(const grid_info &Info, std::string &Name) {

  Name = Info.Name_;

}

void GetGridInfoDimension(const grid_info &Info, int &NumDims) {

  NumDims = Info.NumDims_;

}

void GetGridInfoRootRank(const grid_info &Info, int &RootRank) {

  RootRank = Info.RootRank_;

}

void GetGridInfoCart(const grid_info &Info, cart &Cart) {

  Cart = Info.Cart_;

}


void GetGridInfoPeriodicLength(const grid_info &Info, double *PeriodicLength) {

  OVK_DEBUG_ASSERT(PeriodicLength, "Invalid periodic length pointer.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    PeriodicLength[iDim] = Info.PeriodicLength_[iDim];
  }

}

void GetGridInfoGeometryType(const grid_info &Info, geometry_type &GeometryType) {

  GeometryType = Info.GeometryType_;

}

void GetGridInfoGlobalRange(const grid_info &Info, range &GlobalRange) {

  GlobalRange = Info.GlobalRange_;

}

void GetGridInfoIsLocal(const grid_info &Info, bool &IsLocal) {

  IsLocal = Info.IsLocal_;

}

}
