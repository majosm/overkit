// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Domain.hpp"

#include "ovk/core/ArrayView.hpp"
#include "ovk/core/AssemblyOptions.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/ErrorHandler.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <map>
#include <string>
#include <utility>

namespace ovk {

namespace {

void CreateGridGlobal(domain &Domain, int GridID, const grid_params *Params, bool IsLocal);

void EditGridGlobal(domain &Domain, int GridID, grid *&Grid, bool IsLocal);
void ReleaseGridGlobal(domain &Domain, int GridID, grid *&Grid, bool IsLocal);

void EnableConnectivityComponent(domain &Domain);
void DisableConnectivityComponent(domain &Domain);

void CreateConnectivitiesForGrid(domain &Domain, int GridID);
void DestroyConnectivitiesForGrid(domain &Domain, int GridID);

void CreateConnectivityGlobal(domain &Domain, int DonorGridID, int ReceiverGridID);
void DestroyConnectivityGlobal(domain &Domain, int DonorGridID, int ReceiverGridID);

void EditConnectivityGlobal(domain &Domain, int DonorGridID, int ReceiverGridID,
  connectivity *&Connectivity, bool IsLocal);
void ReleaseConnectivityGlobal(domain &Domain, int DonorGridID, int ReceiverGridID,
  connectivity *&Connectivity, bool IsLocal);

void EnableExchangeComponent(domain &Domain);
void DisableExchangeComponent(domain &Domain);

void CreateExchangesForGrid(domain &Domain, int GridID);
void DestroyExchangesForGrid(domain &Domain, int GridID);

void CreateExchangeGlobal(domain &Domain, int DonorGridID, int ReceiverGridID);
void DestroyExchangeGlobal(domain &Domain, int DonorGridID, int ReceiverGridID);

bool EditingGrid(const domain &Domain, int GridID);
bool EditingConnectivity(const domain &Domain, int DonorGridID, int ReceiverGridID);

void AssembleExchange(domain &Domain);

void ResetAllConnectivityEdits(domain &Domain);

void CreateDomainGridInfo(domain::grid_info &GridInfo, grid *Grid, const core::comm &Comm);
void DestroyDomainGridInfo(domain::grid_info &GridInfo);

void CreateDomainConnectivityInfo(domain::connectivity_info &ConnectivityInfo, connectivity
  *Connectivity, const core::comm &Comm);
void DestroyDomainConnectivityInfo(domain::connectivity_info &ConnectivityInfo);

}

namespace core {

void CreateDomain(domain &Domain, const domain_params &Params, logger &Logger, error_handler
  &ErrorHandler) {

  Domain.Comm_ = core::comm(Params.Comm_);

  MPI_Barrier(Domain.Comm_);

  Domain.Logger_ = &Logger;
  Domain.ErrorHandler_ = &ErrorHandler;

  core::CreateProfiler(Domain.Profiler_, Domain.Comm_);
  if (OVK_TIMERS) core::EnableProfiler(Domain.Profiler_);

  if (Params.Name_.length() > 0) {
    Domain.Name_ = Params.Name_;
  } else {
    Domain.Name_ = "Domain";
  }

  Domain.NumDims_ = Params.NumDims_;

  Domain.Config_ = domain_config::NONE;

  Domain.NumGrids_ = 0;
  Domain.AllGridsEditRefCount_ = 0;
  Domain.AllConnectivitiesEditRefCount_ = 0;

  if (Domain.Comm_.Rank() == 0) {
    std::string ProcessesString = FormatNumber(Domain.Comm_.Size(), "processes", "process");
    core::LogStatus(*Domain.Logger_, true, 0, "Created %1iD domain %s on %s.", Domain.NumDims_,
      Domain.Name_, ProcessesString);
  }

  MPI_Barrier(Domain.Comm_);

}

void DestroyDomain(domain &Domain) {

  MPI_Barrier(Domain.Comm_);

  if ((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE) {
    DisableExchangeComponent(Domain);
  }

  if ((Domain.Config_ & domain_config::CONNECTIVITY) != domain_config::NONE) {
    DisableConnectivityComponent(Domain);
  }

  for (auto &Pair : Domain.GridInfo_) {
    domain::grid_info &Info = Pair.second;
    DestroyDomainGridInfo(Info);
  }
  Domain.GridInfo_.clear();

  for (auto &Pair : Domain.LocalGrids_) {
    grid &Grid = Pair.second;
    DestroyGrid(Grid);
  }
  Domain.LocalGrids_.clear();

  std::string ProfileTimesString = core::WriteProfileTimes(Domain.Profiler_);
  if (Domain.Comm_.Rank() == 0) {
    printf("%s", ProfileTimesString.c_str());
  }
  core::DestroyProfiler(Domain.Profiler_);

  MPI_Barrier(Domain.Comm_);

  core::LogStatus(*Domain.Logger_, Domain.Comm_.Rank() == 0, 0, "Destroyed domain %s.",
    Domain.Name_);

  Domain.Comm_.Reset();

}

}

void GetDomainName(const domain &Domain, std::string &Name) {

  Name = Domain.Name_;

}

void GetDomainDimension(const domain &Domain, int &NumDims) {

  NumDims = Domain.NumDims_;

}

void GetDomainComm(const domain &Domain, MPI_Comm &Comm) {

  Comm = Domain.Comm_.Get();

}

void GetDomainCommSize(const domain &Domain, int &CommSize) {

  CommSize = Domain.Comm_.Size();

}

void GetDomainCommRank(const domain &Domain, int &CommRank) {

  CommRank = Domain.Comm_.Rank();

}

void ConfigureDomain(domain &Domain, domain_config Config) {

  OVK_DEBUG_ASSERT(ValidDomainConfig(Config), "Invalid domain config.");
  OVK_DEBUG_ASSERT(!EditingGrid(Domain, ALL_GRIDS), "Cannot configure domain while editing grids.");
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot configure "
    "domain while editing connectivities.");

  MPI_Barrier(Domain.Comm_);

  domain_config OldConfig = Domain.Config_;

  bool HasGeometry = (Config & domain_config::GEOMETRY) != domain_config::NONE;
  bool HasOverlap = (Config & domain_config::OVERLAP) != domain_config::NONE;
  bool HasConnectivity = (Config & domain_config::CONNECTIVITY) != domain_config::NONE;
  bool HasExchange = (Config & domain_config::EXCHANGE) != domain_config::NONE;

  OVK_DEBUG_ASSERT(!HasOverlap || HasGeometry, "Domain overlap component requires geometry "
    "component.");
  OVK_DEBUG_ASSERT(!HasExchange || HasConnectivity, "Domain exchange component requires "
    "connectivity component.");

  Domain.Config_ = Config;

  domain_config ConfigAdded = ~OldConfig & Config;
  domain_config ConfigRemoved = OldConfig & ~Config;

  if ((ConfigAdded & domain_config::CONNECTIVITY) != domain_config::NONE) {
    EnableConnectivityComponent(Domain);
  } else if ((ConfigRemoved & domain_config::CONNECTIVITY) != domain_config::NONE) {
    DisableConnectivityComponent(Domain);
  }

  if ((ConfigAdded & domain_config::EXCHANGE) != domain_config::NONE) {
    EnableExchangeComponent(Domain);
  } else if ((ConfigRemoved & domain_config::EXCHANGE) != domain_config::NONE) {
    DisableExchangeComponent(Domain);
  }

  MPI_Barrier(Domain.Comm_);

}

void GetDomainConfiguration(const domain &Domain, domain_config &Config) {

  Config = Domain.Config_;

}

void GetDomainGridCount(const domain &Domain, int &NumGrids) {

  NumGrids = Domain.NumGrids_;

}

void GetDomainGridIDs(const domain &Domain, int *GridIDs) {

  int iGrid = 0;
  for (auto &Pair : Domain.GridInfo_) {
    GridIDs[iGrid] = Pair.first;
    ++iGrid;
  }

}

void GetNextAvailableGridID(const domain &Domain, int &GridID) {

  auto Iter = Domain.GridInfo_.cbegin();

  if (Iter == Domain.GridInfo_.cend() || Iter->first > 0) {
    GridID = 0;
    return;
  }

  auto PrevIter = Iter;
  ++Iter;
  while (Iter != Domain.GridInfo_.end()) {
    if (Iter->first - PrevIter->first > 1) {
      GridID = PrevIter->first + 1;
      return;
    }
    PrevIter = Iter;
    ++Iter;
  }

  GridID = PrevIter->first + 1;

}

namespace core {

const comm &GetDomainComm(const domain &Domain) {

  return Domain.Comm_;

}

logger &GetDomainLogger(const domain &Domain) {

  return *Domain.Logger_;

}

error_handler &GetDomainErrorHandler(const domain &Domain) {

  return *Domain.ErrorHandler_;

}

profiler &GetDomainProfiler(const domain &Domain) {

  return Domain.Profiler_;

}

}

void CreateGridLocal(domain &Domain, int GridID, const grid_params &Params) {

  CreateGridGlobal(Domain, GridID, &Params, true);

}

void CreateGridRemote(domain &Domain, int GridID) {

  CreateGridGlobal(Domain, GridID, nullptr, false);

}

namespace {

void CreateGridGlobal(domain &Domain, int GridID, const grid_params *Params, bool IsLocal) {

  MPI_Barrier(Domain.Comm_);

  OVK_DEBUG_ASSERT(!EditingGrid(Domain, ALL_GRIDS), "Cannot create grid while editing other grids.");
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot create grid "
    "while editing connectivities.");

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(!GridExists(Domain, GridID), "Grid %i already exists.", GridID);

  if (OVK_DEBUG) {
    int IsLocalInt = IsLocal ? 1 : 0;
    int AtLeastOneLocal;
    MPI_Allreduce(&IsLocalInt, &AtLeastOneLocal, 1, MPI_INT, MPI_LOR, Domain.Comm_);
    OVK_DEBUG_ASSERT(AtLeastOneLocal, "Grid must exist on at least one process.");
  }

  domain::grid_info GridInfo;

  if (IsLocal) {
    grid Grid;
    core::CreateGrid(Grid, GridID, *Params, *Domain.Logger_, *Domain.ErrorHandler_);
    CreateDomainGridInfo(GridInfo, &Grid, Domain.Comm_);
    Domain.LocalGrids_.emplace(GridID, std::move(Grid));
  } else {
    CreateDomainGridInfo(GridInfo, nullptr, Domain.Comm_);
  }

  Domain.GridInfo_.emplace(GridID, std::move(GridInfo));

  ++Domain.NumGrids_;

  if ((Domain.Config_ & domain_config::CONNECTIVITY) != domain_config::NONE) {
    CreateConnectivitiesForGrid(Domain, GridID);
  }

  if ((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE) {
    CreateExchangesForGrid(Domain, GridID);
  }

  MPI_Barrier(Domain.Comm_);

}

}

void DestroyGrid(domain &Domain, int GridID) {

  MPI_Barrier(Domain.Comm_);

  OVK_DEBUG_ASSERT(!EditingGrid(Domain, ALL_GRIDS), "Cannot destroy grid while editing other grids.");
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot destroy grid "
    "while editing connectivities.");

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, GridID), "Grid %i does not exist.", GridID);

  bool IsLocal = RankHasGrid(Domain, GridID);

  if ((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE) {
    DestroyExchangesForGrid(Domain, GridID);
  }

  if ((Domain.Config_ & domain_config::CONNECTIVITY) != domain_config::NONE) {
    DestroyConnectivitiesForGrid(Domain, GridID);
  }

  DestroyDomainGridInfo(Domain.GridInfo_[GridID]);

  Domain.GridInfo_.erase(GridID);

  if (IsLocal) {
    core::DestroyGrid(Domain.LocalGrids_[GridID]);
    Domain.LocalGrids_.erase(GridID);
  }

  --Domain.NumGrids_;

  MPI_Barrier(Domain.Comm_);

}

bool GridExists(const domain &Domain, int GridID) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");

  return Domain.GridInfo_.find(GridID) != Domain.GridInfo_.end();

}

void GetGridInfo(const domain &Domain, int GridID, const grid_info *&GridInfo) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, GridID), "Grid %i does not exist.", GridID);

  GridInfo = &Domain.GridInfo_.at(GridID);

}

bool RankHasGrid(const domain &Domain, int GridID) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, GridID), "Grid %i does not exist.", GridID);

  return Domain.LocalGrids_.find(GridID) != Domain.LocalGrids_.end();

}

void GetGrid(const domain &Domain, int GridID, const grid *&Grid) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, GridID), "Grid %i does not exist.", GridID);
  if (OVK_DEBUG) {
    const domain::grid_info &GridInfo = Domain.GridInfo_.at(GridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, GridID), "Grid %s does not have local data on rank @rank@.",
      GridName);
  }

  Grid = &Domain.LocalGrids_.at(GridID);

}

void EditGridLocal(domain &Domain, int GridID, grid *&Grid) {

  EditGridGlobal(Domain, GridID, Grid, true);

}

void EditGridRemote(domain &Domain, int GridID) {

  grid *IgnoredGrid = nullptr;
  EditGridGlobal(Domain, GridID, IgnoredGrid, false);

}

void ReleaseGridLocal(domain &Domain, int GridID, grid *&Grid) {

  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  ReleaseGridGlobal(Domain, GridID, Grid, true);

}

void ReleaseGridRemote(domain &Domain, int GridID) {

  grid *IgnoredGrid = nullptr;
  ReleaseGridGlobal(Domain, GridID, IgnoredGrid, false);

}

namespace {

void EditGridGlobal(domain &Domain, int GridID, grid *&Grid, bool IsLocal) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot edit "
    "grid while editing connectivities.");

  domain::grid_info &GridInfo = Domain.GridInfo_[GridID];

  if (OVK_DEBUG) {
    if (IsLocal) {
      std::string GridName;
      GetGridInfoName(GridInfo, GridName);
      OVK_DEBUG_ASSERT(RankHasGrid(Domain, GridID), "Grid %s does not have local data on rank "
        "@rank@.", GridName);
    }
  }

//   bool StartEditAll = Domain.AllGridsEditRefCount_ == 0;
  ++Domain.AllGridsEditRefCount_;
  bool StartEdit = GridInfo.EditRefCount_ == 0;
  ++GridInfo.EditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Domain.Comm_);
  }

  if (IsLocal) {
    Grid = &Domain.LocalGrids_[GridID];
  }

}

void ReleaseGridGlobal(domain &Domain, int GridID, grid *&Grid, bool IsLocal) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, GridID), "Grid %i does not exist.", GridID);

  domain::grid_info &GridInfo = Domain.GridInfo_[GridID];

  if (OVK_DEBUG) {
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(GridInfo.EditRefCount_ > 0, "Unable to release grid %s; not currently being "
      "edited.", GridName);
    if (IsLocal) {
      OVK_DEBUG_ASSERT(RankHasGrid(Domain, GridID), "Grid %s does not have local data on rank "
        "@rank@.", GridName);
      OVK_DEBUG_ASSERT(Grid == &Domain.LocalGrids_[GridID], "Invalid grid pointer.");
    }
  }

  --GridInfo.EditRefCount_;
  bool EndEdit = GridInfo.EditRefCount_ == 0;
  --Domain.AllGridsEditRefCount_;
//   bool EndEditAll = Domain.AllGridsEditRefCount_ == 0;

  if (IsLocal) {
    Grid = nullptr;
  }

  if (EndEdit) {
    MPI_Barrier(Domain.Comm_);
  }

}

void EnableConnectivityComponent(domain &Domain) {

  for (auto &GridNPair : Domain.GridInfo_) {
    int ReceiverGridID = GridNPair.first;
    for (auto &GridMPair : Domain.GridInfo_) {
      int DonorGridID = GridMPair.first;
      if (DonorGridID != ReceiverGridID) {
        CreateConnectivityGlobal(Domain, DonorGridID, ReceiverGridID);
      }
    }
  }

}

void DisableConnectivityComponent(domain &Domain) {

  for (auto &GridNPair : Domain.GridInfo_) {
    int ReceiverGridID = GridNPair.first;
    for (auto &GridMPair : Domain.GridInfo_) {
      int DonorGridID = GridMPair.first;
      if (DonorGridID != ReceiverGridID) {
        DestroyConnectivityGlobal(Domain, DonorGridID, ReceiverGridID);
      }
    }
  }

}

void CreateConnectivitiesForGrid(domain &Domain, int GridID) {

  for (auto &Pair : Domain.GridInfo_) {
    int OtherGridID = Pair.first;
    if (OtherGridID != GridID) {
      CreateConnectivityGlobal(Domain, GridID, OtherGridID);
      CreateConnectivityGlobal(Domain, OtherGridID, GridID);
    }
  }

}

void DestroyConnectivitiesForGrid(domain &Domain, int GridID) {

  for (auto &Pair : Domain.GridInfo_) {
    int OtherGridID = Pair.first;
    if (OtherGridID != GridID) {
      DestroyConnectivityGlobal(Domain, GridID, OtherGridID);
      DestroyConnectivityGlobal(Domain, OtherGridID, GridID);
    }
  }

}

void CreateConnectivityGlobal(domain &Domain, int DonorGridID, int ReceiverGridID) {

  MPI_Barrier(Domain.Comm_);

  bool DonorGridIsLocal = RankHasGrid(Domain, DonorGridID);
  bool ReceiverGridIsLocal = RankHasGrid(Domain, ReceiverGridID);

  bool IsLocal = DonorGridIsLocal || ReceiverGridIsLocal;

  core::comm ConnectivityComm = core::CreateSubsetComm(Domain.Comm_, IsLocal);

  domain::connectivity_info ConnectivityInfo;

  if (IsLocal) {
    const grid *DonorGrid = nullptr;
    if (DonorGridIsLocal) DonorGrid = &Domain.LocalGrids_[DonorGridID];
    const grid *ReceiverGrid = nullptr;
    if (ReceiverGridIsLocal) ReceiverGrid = &Domain.LocalGrids_[ReceiverGridID];
    connectivity Connectivity;
    core::CreateConnectivity(Connectivity, Domain.NumDims_, std::move(ConnectivityComm), DonorGrid,
      ReceiverGrid, *Domain.Logger_, *Domain.ErrorHandler_);
    CreateDomainConnectivityInfo(ConnectivityInfo, &Connectivity, Domain.Comm_);
    Domain.LocalConnectivities_[DonorGridID].emplace(ReceiverGridID, std::move(Connectivity));
  } else {
    CreateDomainConnectivityInfo(ConnectivityInfo, nullptr, Domain.Comm_);
  }

  Domain.ConnectivityInfo_[DonorGridID].emplace(ReceiverGridID, std::move(ConnectivityInfo));

  MPI_Barrier(Domain.Comm_);

}

void DestroyConnectivityGlobal(domain &Domain, int DonorGridID, int ReceiverGridID) {

  MPI_Barrier(Domain.Comm_);

  bool IsLocal = RankHasConnectivity(Domain, DonorGridID, ReceiverGridID);

  DestroyDomainConnectivityInfo(Domain.ConnectivityInfo_[DonorGridID][ReceiverGridID]);

  Domain.ConnectivityInfo_[DonorGridID].erase(ReceiverGridID);
  if (Domain.ConnectivityInfo_[DonorGridID].size() == 0) {
    Domain.ConnectivityInfo_.erase(DonorGridID);
  }

  if (IsLocal) {
    core::DestroyConnectivity(Domain.LocalConnectivities_[DonorGridID][ReceiverGridID]);
    Domain.LocalConnectivities_[DonorGridID].erase(ReceiverGridID);
    if (Domain.LocalConnectivities_[DonorGridID].size() == 0) {
      Domain.LocalConnectivities_.erase(DonorGridID);
    }
  }

  MPI_Barrier(Domain.Comm_);

}

}

bool ConnectivityExists(const domain &Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(GridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);

  auto RowIter = Domain.ConnectivityInfo_.find(DonorGridID);
  if (RowIter != Domain.ConnectivityInfo_.end()) {
    const std::map<int, domain::connectivity_info> &Row = RowIter->second;
    return Row.find(ReceiverGridID) != Row.end();
  } else {
    return false;
  }

}

void GetConnectivityInfo(const domain &Domain, int DonorGridID, int ReceiverGridID, const
  connectivity_info *&ConnectivityInfo) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  ConnectivityInfo = &Domain.ConnectivityInfo_.at(DonorGridID).at(ReceiverGridID);

}

bool RankHasConnectivity(const domain &Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(GridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  auto RowIter = Domain.LocalConnectivities_.find(DonorGridID);
  if (RowIter != Domain.LocalConnectivities_.end()) {
    const std::map<int, connectivity> &Row = RowIter->second;
    return Row.find(ReceiverGridID) != Row.end();
  } else {
    return false;
  }

}

void GetConnectivity(const domain &Domain, int DonorGridID, int ReceiverGridID, const connectivity
  *&Connectivity) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);
  if (OVK_DEBUG) {
    const domain::connectivity_info &ConnectivityInfo = Domain.ConnectivityInfo_.at(DonorGridID)
      .at(ReceiverGridID);
    std::string Name;
    GetConnectivityInfoName(ConnectivityInfo, Name);
    OVK_DEBUG_ASSERT(RankHasConnectivity(Domain, DonorGridID, ReceiverGridID), "Connectivity %s "
      "does not have local data on rank @rank@.", Name);
  }

  Connectivity = &Domain.LocalConnectivities_.at(DonorGridID).at(ReceiverGridID);

}

void EditConnectivityLocal(domain &Domain, int DonorGridID, int ReceiverGridID,
  connectivity *&Connectivity) {

  EditConnectivityGlobal(Domain, DonorGridID, ReceiverGridID, Connectivity, true);

}

void EditConnectivityRemote(domain &Domain, int DonorGridID, int ReceiverGridID) {

  connectivity *IgnoredConnectivity = nullptr;
  EditConnectivityGlobal(Domain, DonorGridID, ReceiverGridID, IgnoredConnectivity, false);

}

void ReleaseConnectivityLocal(domain &Domain, int DonorGridID, int ReceiverGridID,
  connectivity *&Connectivity) {

  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  ReleaseConnectivityGlobal(Domain, DonorGridID, ReceiverGridID, Connectivity, true);

}

void ReleaseConnectivityRemote(domain &Domain, int DonorGridID, int ReceiverGridID) {

  connectivity *IgnoredConnectivity = nullptr;
  ReleaseConnectivityGlobal(Domain, DonorGridID, ReceiverGridID, IgnoredConnectivity, false);

}

namespace {

void EditConnectivityGlobal(domain &Domain, int DonorGridID, int ReceiverGridID,
  connectivity *&Connectivity, bool IsLocal) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity "
    "(%i,%i) does not exist.", DonorGridID, ReceiverGridID);

  OVK_DEBUG_ASSERT(!EditingGrid(Domain, OVK_ALL_GRIDS), "Cannot edit connectivity while editing "
    "grids.");

  domain::connectivity_info &ConnectivityInfo = Domain.ConnectivityInfo_[DonorGridID][ReceiverGridID];

  if (OVK_DEBUG) {
    if (IsLocal) {
      std::string Name;
      GetConnectivityInfoName(ConnectivityInfo, Name);
      OVK_DEBUG_ASSERT(RankHasConnectivity(Domain, DonorGridID, ReceiverGridID), "Connectivity %s "
        "does not have local data on rank @rank@.", Name);
    }
  }

//   bool StartEditAll = Domain.AllConnectivitiesEditRefCount_ == 0;
  ++Domain.AllConnectivitiesEditRefCount_;
  bool StartEdit = ConnectivityInfo.EditRefCount_ == 0;
  ++ConnectivityInfo.EditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Domain.Comm_);
  }

  if (IsLocal) {
    Connectivity = &Domain.LocalConnectivities_[DonorGridID][ReceiverGridID];
  }

}

void ReleaseConnectivityGlobal(domain &Domain, int DonorGridID, int ReceiverGridID, connectivity
  *&Connectivity, bool IsLocal) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(GridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  domain::connectivity_info &ConnectivityInfo = Domain.ConnectivityInfo_[DonorGridID][ReceiverGridID];

  if (OVK_DEBUG) {
    std::string Name;
    GetConnectivityInfoName(ConnectivityInfo, Name);
    OVK_DEBUG_ASSERT(ConnectivityInfo.EditRefCount_ > 0, "Unable to release connectivity %s; not "
      "currently being edited.", Name);
    if (IsLocal) {
      OVK_DEBUG_ASSERT(RankHasConnectivity(Domain, DonorGridID, ReceiverGridID), "Connectivity %s "
        "does not have local data on rank @rank@.", Name);
      OVK_DEBUG_ASSERT(Connectivity == &Domain.LocalConnectivities_[DonorGridID][ReceiverGridID],
        "Invalid connectivity pointer.");
    }
  }

  --ConnectivityInfo.EditRefCount_;
  bool EndEdit = ConnectivityInfo.EditRefCount_ == 0;
  --Domain.AllConnectivitiesEditRefCount_;
//   bool EndEditAll = Domain.AllConnectivitiesEditRefCount_ == 0;

  if (IsLocal) {
    Connectivity = nullptr;
  }

  if (EndEdit) {
    MPI_Barrier(Domain.Comm_);
  }

}

void EnableExchangeComponent(domain &Domain) {

  for (auto &GridNPair : Domain.GridInfo_) {
    int ReceiverGridID = GridNPair.first;
    for (auto &GridMPair : Domain.GridInfo_) {
      int DonorGridID = GridMPair.first;
      if (DonorGridID != ReceiverGridID) {
        CreateExchangeGlobal(Domain, DonorGridID, ReceiverGridID);
      }
    }
  }

}

void DisableExchangeComponent(domain &Domain) {

  for (auto &GridNPair : Domain.GridInfo_) {
    int ReceiverGridID = GridNPair.first;
    for (auto &GridMPair : Domain.GridInfo_) {
      int DonorGridID = GridMPair.first;
      if (DonorGridID != ReceiverGridID) {
        DestroyExchangeGlobal(Domain, DonorGridID, ReceiverGridID);
      }
    }
  }

}

void CreateExchangesForGrid(domain &Domain, int GridID) {

  for (auto &Pair : Domain.GridInfo_) {
    int OtherGridID = Pair.first;
    if (OtherGridID != GridID) {
      CreateExchangeGlobal(Domain, GridID, OtherGridID);
      CreateExchangeGlobal(Domain, OtherGridID, GridID);
    }
  }

}

void DestroyExchangesForGrid(domain &Domain, int GridID) {

  for (auto &Pair : Domain.GridInfo_) {
    int OtherGridID = Pair.first;
    if (OtherGridID != GridID) {
      DestroyExchangeGlobal(Domain, GridID, OtherGridID);
      DestroyExchangeGlobal(Domain, OtherGridID, GridID);
    }
  }

}

void CreateExchangeGlobal(domain &Domain, int DonorGridID, int ReceiverGridID) {

  MPI_Barrier(Domain.Comm_);

  bool IsLocal = RankHasConnectivity(Domain, DonorGridID, ReceiverGridID);

  exchange_info ExchangeInfo;

  if (IsLocal) {
    const connectivity &Connectivity = Domain.LocalConnectivities_[DonorGridID][ReceiverGridID];
    exchange Exchange;
    core::CreateExchange(Exchange, Connectivity, *Domain.Logger_, *Domain.ErrorHandler_,
      Domain.Profiler_);
    core::CreateExchangeInfo(ExchangeInfo, &Exchange, Domain.Comm_);
    Domain.LocalExchanges_[DonorGridID].emplace(ReceiverGridID, std::move(Exchange));
  } else {
    core::CreateExchangeInfo(ExchangeInfo, nullptr, Domain.Comm_);
  }

  Domain.ExchangeInfo_[DonorGridID].emplace(ReceiverGridID, std::move(ExchangeInfo));

  MPI_Barrier(Domain.Comm_);

}

void DestroyExchangeGlobal(domain &Domain, int DonorGridID, int ReceiverGridID) {

  MPI_Barrier(Domain.Comm_);

  bool IsLocal = RankHasExchange(Domain, DonorGridID, ReceiverGridID);

  core::DestroyExchangeInfo(Domain.ExchangeInfo_[DonorGridID][ReceiverGridID]);

  Domain.ExchangeInfo_[DonorGridID].erase(ReceiverGridID);
  if (Domain.ExchangeInfo_[DonorGridID].size() == 0) {
    Domain.ExchangeInfo_.erase(DonorGridID);
  }

  if (IsLocal) {
    core::DestroyExchange(Domain.LocalExchanges_[DonorGridID][ReceiverGridID]);
    Domain.LocalExchanges_[DonorGridID].erase(ReceiverGridID);
    if (Domain.LocalExchanges_[DonorGridID].size() == 0) {
      Domain.LocalExchanges_.erase(DonorGridID);
    }
  }

  MPI_Barrier(Domain.Comm_);

}

}

bool ExchangeExists(const domain &Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(GridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);

  auto RowIter = Domain.ExchangeInfo_.find(DonorGridID);
  if (RowIter != Domain.ExchangeInfo_.end()) {
    const std::map<int, exchange_info> &Row = RowIter->second;
    return Row.find(ReceiverGridID) != Row.end();
  } else {
    return false;
  }

}

void GetExchangeInfo(const domain &Domain, int DonorGridID, int ReceiverGridID, const exchange_info
  *&ExchangeInfo) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
    "exist.", DonorGridID, ReceiverGridID);

  ExchangeInfo = &Domain.ExchangeInfo_.at(DonorGridID).at(ReceiverGridID);

}

bool RankHasExchange(const domain &Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(GridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(GridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
    "exist.", DonorGridID, ReceiverGridID);

  auto RowIter = Domain.LocalExchanges_.find(DonorGridID);
  if (RowIter != Domain.LocalExchanges_.end()) {
    const std::map<int, exchange> &Row = RowIter->second;
    return Row.find(ReceiverGridID) != Row.end();
  } else {
    return false;
  }

}

void GetExchange(const domain &Domain, int DonorGridID, int ReceiverGridID, const exchange
  *&Exchange) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
    "exist.", DonorGridID, ReceiverGridID);
  if (OVK_DEBUG) {
    const exchange_info &ExchangeInfo = Domain.ExchangeInfo_.at(DonorGridID).at(ReceiverGridID);
    std::string Name;
    GetExchangeInfoName(ExchangeInfo, Name);
    OVK_DEBUG_ASSERT(RankHasExchange(Domain, DonorGridID, ReceiverGridID), "Exchange %s does not "
      "have local data on rank @rank@.", Name);
  }

  Exchange = &Domain.LocalExchanges_.at(DonorGridID).at(ReceiverGridID);

}

namespace {

bool EditingGrid(const domain &Domain, int GridID) {

  if (GridID == ALL_GRIDS) {
    return Domain.AllGridsEditRefCount_ > 0;
  } else {
    return Domain.GridInfo_.at(GridID).EditRefCount_ > 0;
  }

}

bool EditingConnectivity(const domain &Domain, int DonorGridID, int ReceiverGridID) {

  if (DonorGridID == OVK_ALL_GRIDS && ReceiverGridID == OVK_ALL_GRIDS) {
    return Domain.AllConnectivitiesEditRefCount_ > 0;
  } else {
    return Domain.ConnectivityInfo_.at(DonorGridID).at(ReceiverGridID).EditRefCount_ > 0;
  }

}

}

void GetLocalDonorCount(const domain &Domain, int DonorGridID, int ReceiverGridID, long long
  &NumDonors) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (RankHasConnectivity(Domain, DonorGridID, ReceiverGridID)) {
    const connectivity *Connectivity;
    GetConnectivity(Domain, DonorGridID, ReceiverGridID, Connectivity);
    if (RankHasConnectivityDonorSide(*Connectivity)) {
      const connectivity_d *Donors;
      GetConnectivityDonorSide(*Connectivity, Donors);
      GetConnectivityDonorSideCount(*Donors, NumDonors);
    } else {
      NumDonors = 0;
    }
  } else {
    NumDonors = 0;
  }

}

void GetLocalReceiverCount(const domain &Domain, int DonorGridID, int ReceiverGridID,
  long long &NumReceivers) {

  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (RankHasConnectivity(Domain, DonorGridID, ReceiverGridID)) {
    const connectivity *Connectivity;
    GetConnectivity(Domain, DonorGridID, ReceiverGridID, Connectivity);
    if (RankHasConnectivityReceiverSide(*Connectivity)) {
      const connectivity_r *Receivers;
      GetConnectivityReceiverSide(*Connectivity, Receivers);
      GetConnectivityReceiverSideCount(*Receivers, NumReceivers);
    } else {
      NumReceivers = 0;
    }
  } else {
    NumReceivers = 0;
  }

}

void Assemble(domain &Domain, const assembly_options &Options) {

  bool IsDomainRoot = Domain.Comm_.Rank() == 0;

  core::LogStatus(*Domain.Logger_, IsDomainRoot, 0, "Beginning overset assembly on domain %s.",
    Domain.Name_);

//   bool HasOverlap = (Domain.Config_ & domain_config::OVERLAP) != domain_config::NONE;
  bool HasConnectivity = (Domain.Config_ & domain_config::CONNECTIVITY) != domain_config::NONE;
  bool HasExchange = (Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE;

//   if (HasOverlap) {
// //     AssembleOverlap(Domain);
//   }

//   if (HasOverlap && HasConnectivity) {
// //     AssembleConnectivity(Domain);
//   }

  if (HasExchange) {
    AssembleExchange(Domain);
  }

  if (HasConnectivity) {
    ResetAllConnectivityEdits(Domain);
  }

  core::LogStatus(*Domain.Logger_, IsDomainRoot, 0, "Finished overset assembly on domain %s.",
    Domain.Name_);

}

namespace {

void AssembleExchange(domain &Domain) {

  for (auto &MPair : Domain.LocalExchanges_) {
    for (auto &NPair : MPair.second) {
      exchange &Exchange = NPair.second;
      core::UpdateExchange(Exchange);
    }
  }

}

void ResetAllConnectivityEdits(domain &Domain) {

  for (auto &MPair : Domain.LocalConnectivities_) {
    for (auto &NPair : MPair.second) {
      connectivity &Connectivity = NPair.second;
      core::ResetConnectivityEdits(Connectivity);
    }
  }

}

}

void Collect(const domain &Domain, int DonorGridID, int ReceiverGridID, data_type ValueType,
  int Count, collect_op CollectOp, const range &GridValuesRange, array_layout GridValuesLayout,
  const void * const *GridValues, void **DonorValues) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(ValidCollectOp(CollectOp), "Invalid collect operation.");
  OVK_DEBUG_ASSERT(ValidArrayLayout(GridValuesLayout), "Invalid grid values layout.");
  OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
    "exist.", DonorGridID, ReceiverGridID);
  if (OVK_DEBUG) {
    const exchange_info &ExchangeInfo = Domain.ExchangeInfo_.at(DonorGridID).at(ReceiverGridID);
    std::string Name;
    GetExchangeInfoName(ExchangeInfo, Name);
    OVK_DEBUG_ASSERT(RankHasExchange(Domain, DonorGridID, ReceiverGridID), "Exchange %s does not "
      "have local data on rank @rank@.", Name);
  }

  int CollectTime = core::GetProfilerTimerID(Domain.Profiler_, "Collect");
  const core::comm &Comm = core::GetGridComm(Domain.LocalGrids_.at(DonorGridID));
  core::StartProfileSync(Domain.Profiler_, CollectTime, Comm);

  const exchange &Exchange = Domain.LocalExchanges_.at(DonorGridID).at(ReceiverGridID);

  core::Collect(Exchange, ValueType, Count, CollectOp, GridValuesRange, GridValuesLayout,
    GridValues, DonorValues);

  core::EndProfile(Domain.Profiler_, CollectTime);

}

void Send(const domain &Domain, int DonorGridID, int ReceiverGridID, data_type ValueType, int Count,
  const void * const *DonorValues, int Tag, request &Request) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(Tag >= 0, "Invalid tag.");
  OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
    "exist.", DonorGridID, ReceiverGridID);
  if (OVK_DEBUG) {
    const exchange_info &ExchangeInfo = Domain.ExchangeInfo_.at(DonorGridID).at(ReceiverGridID);
    std::string Name;
    GetExchangeInfoName(ExchangeInfo, Name);
    OVK_DEBUG_ASSERT(RankHasExchange(Domain, DonorGridID, ReceiverGridID), "Exchange %s does not "
      "have local data on rank @rank@.", Name);
  }

  const exchange &Exchange = Domain.LocalExchanges_.at(DonorGridID).at(ReceiverGridID);

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  core::Send(Exchange, ValueType, Count, DonorValues, Tag, Request);

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

void Receive(const domain &Domain, int DonorGridID, int ReceiverGridID, data_type ValueType,
  int Count, void **ReceiverValues, int Tag, request &Request) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(Tag >= 0, "Invalid tag.");
  OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
    "exist.", DonorGridID, ReceiverGridID);
  if (OVK_DEBUG) {
    const exchange_info &ExchangeInfo = Domain.ExchangeInfo_.at(DonorGridID).at(ReceiverGridID);
    std::string Name;
    GetExchangeInfoName(ExchangeInfo, Name);
    OVK_DEBUG_ASSERT(RankHasExchange(Domain, DonorGridID, ReceiverGridID), "Exchange %s does not "
      "have local data on rank @rank@.", Name);
  }

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  const exchange &Exchange = Domain.LocalExchanges_.at(DonorGridID).at(ReceiverGridID);

  core::Receive(Exchange, ValueType, Count, ReceiverValues, Tag, Request);

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

void Wait(const domain &Domain, request &Request) {

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  if (Request) {
    Request.Wait();
  }

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

void WaitAll(const domain &Domain, array_view<request> Requests) {

  OVK_DEBUG_ASSERT(Requests || Requests.Count() == 0, "Invalid requests array.");

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  RequestWaitAll(Requests);

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

void WaitAny(const domain &Domain, array_view<request> Requests, int &Index) {

  OVK_DEBUG_ASSERT(Requests || Requests.Count() == 0, "Invalid requests array.");

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  RequestWaitAny(Requests, Index);

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

void WaitAll(const domain &Domain, array_view<request *> Requests) {

  OVK_DEBUG_ASSERT(Requests || Requests.Count() == 0, "Invalid requests array.");

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  RequestWaitAll(Requests);

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

void WaitAny(const domain &Domain, array_view<request *> Requests, int &Index) {

  OVK_DEBUG_ASSERT(Requests || Requests.Count() == 0, "Invalid requests array.");

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  RequestWaitAny(Requests, Index);

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

void Disperse(const domain &Domain, int DonorGridID, int ReceiverGridID, data_type ValueType,
  int Count, disperse_op DisperseOp, const void * const *ReceiverValues, const range
  &GridValuesRange, array_layout GridValuesLayout, void **GridValues) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(ValidDisperseOp(DisperseOp), "Invalid disperse operation.");
  OVK_DEBUG_ASSERT(ValidArrayLayout(GridValuesLayout), "Invalid grid values layout.");
  OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
    "exist.", DonorGridID, ReceiverGridID);
  if (OVK_DEBUG) {
    const exchange_info &ExchangeInfo = Domain.ExchangeInfo_.at(DonorGridID).at(ReceiverGridID);
    std::string Name;
    GetExchangeInfoName(ExchangeInfo, Name);
    OVK_DEBUG_ASSERT(RankHasExchange(Domain, DonorGridID, ReceiverGridID), "Exchange %s does not "
      "have local data on rank @rank@.", Name);
  }

  int DisperseTime = core::GetProfilerTimerID(Domain.Profiler_, "Disperse");
  core::StartProfile(Domain.Profiler_, DisperseTime);

  const exchange &Exchange = Domain.LocalExchanges_.at(DonorGridID).at(ReceiverGridID);

  core::Disperse(Exchange, ValueType, Count, DisperseOp, ReceiverValues, GridValuesRange,
    GridValuesLayout, GridValues);

  core::EndProfile(Domain.Profiler_, DisperseTime);

}

namespace {

void CreateDomainGridInfo(domain::grid_info &GridInfo, grid *Grid, const core::comm &Comm) {

  core::CreateGridInfo(GridInfo, Grid, Comm);
  GridInfo.EditRefCount_ = 0;

}

void DestroyDomainGridInfo(domain::grid_info &GridInfo) {

  core::DestroyGridInfo(GridInfo);

}

void CreateDomainConnectivityInfo(domain::connectivity_info &ConnectivityInfo, connectivity
  *Connectivity, const core::comm &Comm) {

  core::CreateConnectivityInfo(ConnectivityInfo, Connectivity, Comm);
  ConnectivityInfo.EditRefCount_ = 0;

}

void DestroyDomainConnectivityInfo(domain::connectivity_info &ConnectivityInfo) {

  core::DestroyConnectivityInfo(ConnectivityInfo);

}

}

void CreateDomainParams(domain_params &Params, int NumDims) {

  OVK_DEBUG_ASSERT(NumDims == 2 || NumDims == 3, "Invalid dimension.");

  Params.NumDims_ = NumDims;
  Params.Comm_ = MPI_COMM_NULL;

}

void DestroyDomainParams(domain_params &) {}

void GetDomainParamName(const domain_params &Params, std::string &Name) {

  Name = Params.Name_;

}

void SetDomainParamName(domain_params &Params, std::string Name) {

  Params.Name_ = std::move(Name);

}

void GetDomainParamDimension(const domain_params &Params, int &NumDims) {

  NumDims = Params.NumDims_;

}

void GetDomainParamComm(const domain_params &Params, MPI_Comm &Comm) {

  Comm = Params.Comm_;

}

void SetDomainParamComm(domain_params &Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Params.Comm_ = Comm;

}

}
