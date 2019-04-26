// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Domain.hpp"

#include "ovk/core/ArrayView.hpp"
#include "ovk/core/AssemblyOptions.hpp"
#include "ovk/core/Collect.hpp"
#include "ovk/core/CollectMap.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/ConnectivityD.hpp"
#include "ovk/core/ConnectivityR.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Disperse.hpp"
#include "ovk/core/ErrorHandler.hpp"
// #include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Recv.hpp"
#include "ovk/core/RecvMap.hpp"
#include "ovk/core/Request.hpp"
#include "ovk/core/Send.hpp"
#include "ovk/core/SendMap.hpp"
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

// void CreateExchangeGlobal(domain &Domain, int DonorGridID, int ReceiverGridID);
// void DestroyExchangeGlobal(domain &Domain, int DonorGridID, int ReceiverGridID);

bool EditingGrid(const domain &Domain, int GridID);
bool EditingConnectivity(const domain &Domain, int DonorGridID, int ReceiverGridID);

void AssembleConnectivity(domain &Domain);
void UpdateSourceDestRanks(connectivity &Connectivity);

void AssembleExchange(domain &Domain);

void ResetAllConnectivityEdits(domain &Domain);

template <typename T> int GetNextAvailableID(const std::map<int, T> &Map);

array<long long> GetSendRecvOrder(const array<int,2> &ReceiverPoints, const range
  &ReceiverGridGlobalRange);

void CreateDomainGridInfo(domain::grid_info &GridInfo, grid *Grid, core::comm_view Comm);
void DestroyDomainGridInfo(domain::grid_info &GridInfo);

void CreateDomainConnectivityInfo(domain::connectivity_info &ConnectivityInfo, connectivity
  *Connectivity, core::comm_view Comm);
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
    Domain.Logger_->LogStatus(true, 0, "Created %1iD domain %s on %s.", Domain.NumDims_,
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

  Domain.LocalGrids_.clear();

  std::string ProfileTimesString = core::WriteProfileTimes(Domain.Profiler_);
  if (Domain.Comm_.Rank() == 0) {
    printf("%s", ProfileTimesString.c_str());
  }
  core::DestroyProfiler(Domain.Profiler_);

  MPI_Barrier(Domain.Comm_);

  Domain.Logger_->LogStatus(Domain.Comm_.Rank() == 0, 0, "Destroyed domain %s.", Domain.Name_);

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

  GridID = GetNextAvailableID(Domain.GridInfo_);

}

namespace core {

comm_view GetDomainComm(const domain &Domain) {

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
    grid Grid(GridID, *Params, *Domain.Logger_, *Domain.ErrorHandler_, Domain.Profiler_);
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
    Grid = &Domain.LocalGrids_.at(GridID);
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
      OVK_DEBUG_ASSERT(Grid == &Domain.LocalGrids_.at(GridID), "Invalid grid pointer.");
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
    if (DonorGridIsLocal) DonorGrid = &Domain.LocalGrids_.at(DonorGridID);
    const grid *ReceiverGrid = nullptr;
    if (ReceiverGridIsLocal) ReceiverGrid = &Domain.LocalGrids_.at(ReceiverGridID);
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
  if (Domain.ConnectivityInfo_[DonorGridID].empty()) {
    Domain.ConnectivityInfo_.erase(DonorGridID);
  }

  if (IsLocal) {
    core::DestroyConnectivity(Domain.LocalConnectivities_[DonorGridID][ReceiverGridID]);
    Domain.LocalConnectivities_[DonorGridID].erase(ReceiverGridID);
    if (Domain.LocalConnectivities_[DonorGridID].empty()) {
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

  AddProfilerTimer(Domain.Profiler_, "Collect");
  AddProfilerTimer(Domain.Profiler_, "Collect::MemAlloc");
  AddProfilerTimer(Domain.Profiler_, "Collect::MPI");
  AddProfilerTimer(Domain.Profiler_, "Collect::Pack");
  AddProfilerTimer(Domain.Profiler_, "Collect::Reduce");
  AddProfilerTimer(Domain.Profiler_, "SendRecv");
  AddProfilerTimer(Domain.Profiler_, "SendRecv::MemAlloc");
  AddProfilerTimer(Domain.Profiler_, "SendRecv::Pack");
  AddProfilerTimer(Domain.Profiler_, "SendRecv::MPI");
  AddProfilerTimer(Domain.Profiler_, "SendRecv::Unpack");
  AddProfilerTimer(Domain.Profiler_, "Disperse");

//   for (auto &GridNPair : Domain.GridInfo_) {
//     int ReceiverGridID = GridNPair.first;
//     for (auto &GridMPair : Domain.GridInfo_) {
//       int DonorGridID = GridMPair.first;
//       if (DonorGridID != ReceiverGridID) {
//         CreateExchangeGlobal(Domain, DonorGridID, ReceiverGridID);
//       }
//     }
//   }

}

void DisableExchangeComponent(domain &Domain) {

//   for (auto &GridNPair : Domain.GridInfo_) {
//     int ReceiverGridID = GridNPair.first;
//     for (auto &GridMPair : Domain.GridInfo_) {
//       int DonorGridID = GridMPair.first;
//       if (DonorGridID != ReceiverGridID) {
//         DestroyExchangeGlobal(Domain, DonorGridID, ReceiverGridID);
//       }
//     }
//   }

  int CollectTime = core::GetProfilerTimerID(Domain.Profiler_, "Collect");
  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  int DisperseTime = core::GetProfilerTimerID(Domain.Profiler_, "Disperse");

  core::StartProfile(Domain.Profiler_, CollectTime);
  Domain.CollectData_.clear();
  core::EndProfile(Domain.Profiler_, CollectTime);
  core::StartProfile(Domain.Profiler_, SendRecvTime);
  Domain.SendData_.clear();
  Domain.RecvData_.clear();
  core::EndProfile(Domain.Profiler_, SendRecvTime);
  core::StartProfile(Domain.Profiler_, DisperseTime);
  Domain.DisperseData_.clear();
  core::EndProfile(Domain.Profiler_, DisperseTime);

}

void CreateExchangesForGrid(domain &Domain, int GridID) {

//   for (auto &Pair : Domain.GridInfo_) {
//     int OtherGridID = Pair.first;
//     if (OtherGridID != GridID) {
//       CreateExchangeGlobal(Domain, GridID, OtherGridID);
//       CreateExchangeGlobal(Domain, OtherGridID, GridID);
//     }
//   }

}

void DestroyExchangesForGrid(domain &Domain, int GridID) {

//   for (auto &Pair : Domain.GridInfo_) {
//     int OtherGridID = Pair.first;
//     if (OtherGridID != GridID) {
//       DestroyExchangeGlobal(Domain, GridID, OtherGridID);
//       DestroyExchangeGlobal(Domain, OtherGridID, GridID);
//     }
//   }

  int CollectTime = core::GetProfilerTimerID(Domain.Profiler_, "Collect");
  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  int DisperseTime = core::GetProfilerTimerID(Domain.Profiler_, "Disperse");

  core::StartProfile(Domain.Profiler_, CollectTime);
  Domain.CollectData_.erase(GridID);
  auto CollectDataRowIter = Domain.CollectData_.begin();
  while (CollectDataRowIter != Domain.CollectData_.end()) {
    std::map<int, domain::collect_data> &CollectDataRow = CollectDataRowIter->second;
    CollectDataRow.erase(GridID);
    if (CollectDataRow.empty()) {
      CollectDataRowIter = Domain.CollectData_.erase(CollectDataRowIter);
    } else {
      ++CollectDataRowIter;
    }
  }
  core::EndProfile(Domain.Profiler_, CollectTime);

  core::StartProfile(Domain.Profiler_, SendRecvTime);
  Domain.SendData_.erase(GridID);
  auto SendDataRowIter = Domain.SendData_.begin();
  while (SendDataRowIter != Domain.SendData_.end()) {
    std::map<int, domain::send_data> &SendDataRow = SendDataRowIter->second;
    SendDataRow.erase(GridID);
    if (SendDataRow.empty()) {
      SendDataRowIter = Domain.SendData_.erase(SendDataRowIter);
    } else {
      ++SendDataRowIter;
    }
  }
  Domain.RecvData_.erase(GridID);
  auto RecvDataRowIter = Domain.RecvData_.begin();
  while (RecvDataRowIter != Domain.RecvData_.end()) {
    std::map<int, domain::recv_data> &RecvDataRow = RecvDataRowIter->second;
    RecvDataRow.erase(GridID);
    if (RecvDataRow.empty()) {
      RecvDataRowIter = Domain.RecvData_.erase(RecvDataRowIter);
    } else {
      ++RecvDataRowIter;
    }
  }
  core::EndProfile(Domain.Profiler_, SendRecvTime);

  core::StartProfile(Domain.Profiler_, DisperseTime);
  Domain.DisperseData_.erase(GridID);
  auto DisperseDataRowIter = Domain.DisperseData_.begin();
  while (DisperseDataRowIter != Domain.DisperseData_.end()) {
    std::map<int, domain::disperse_data> &DisperseDataRow = DisperseDataRowIter->second;
    DisperseDataRow.erase(GridID);
    if (DisperseDataRow.empty()) {
      DisperseDataRowIter = Domain.DisperseData_.erase(DisperseDataRowIter);
    } else {
      ++DisperseDataRowIter;
    }
  }
  core::EndProfile(Domain.Profiler_, DisperseTime);

}

// void CreateExchangeGlobal(domain &Domain, int DonorGridID, int ReceiverGridID) {

//   MPI_Barrier(Domain.Comm_);

//   bool IsLocal = RankHasConnectivity(Domain, DonorGridID, ReceiverGridID);

//   exchange_info ExchangeInfo;

//   if (IsLocal) {
//     const connectivity &Connectivity = Domain.LocalConnectivities_[DonorGridID][ReceiverGridID];
//     exchange Exchange;
//     core::CreateExchange(Exchange, Connectivity, *Domain.Logger_, *Domain.ErrorHandler_,
//       Domain.Profiler_);
//     core::CreateExchangeInfo(ExchangeInfo, &Exchange, Domain.Comm_);
//     Domain.LocalExchanges_[DonorGridID].emplace(ReceiverGridID, std::move(Exchange));
//     if (RankHasConnectivityDonorSide(Connectivity)) {
//       const connectivity_d *Donors;
//       GetConnectivityDonorSide(Connectivity, Donors);
//       const grid *DonorGrid;
//       GetConnectivityDonorSideGrid(*Donors, DonorGrid);
//     }
//   } else {
//     core::CreateExchangeInfo(ExchangeInfo, nullptr, Domain.Comm_);
//   }

//   Domain.ExchangeInfo_[DonorGridID].emplace(ReceiverGridID, std::move(ExchangeInfo));

//   MPI_Barrier(Domain.Comm_);

// }

// void DestroyExchangeGlobal(domain &Domain, int DonorGridID, int ReceiverGridID) {

//   MPI_Barrier(Domain.Comm_);

//   bool IsLocal = RankHasExchange(Domain, DonorGridID, ReceiverGridID);

//   core::DestroyExchangeInfo(Domain.ExchangeInfo_[DonorGridID][ReceiverGridID]);

//   Domain.ExchangeInfo_[DonorGridID].erase(ReceiverGridID);
//   if (Domain.ExchangeInfo_[DonorGridID].empty()) {
//     Domain.ExchangeInfo_.erase(DonorGridID);
//   }

//   if (IsLocal) {
//     core::DestroyExchange(Domain.LocalExchanges_[DonorGridID][ReceiverGridID]);
//     Domain.LocalExchanges_[DonorGridID].erase(ReceiverGridID);
//     if (Domain.LocalExchanges_[DonorGridID].empty()) {
//       Domain.LocalExchanges_.erase(DonorGridID);
//     }
//   }

//   MPI_Barrier(Domain.Comm_);

// }

}

// bool ExchangeExists(const domain &Domain, int DonorGridID, int ReceiverGridID) {

//   OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
//   OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
//   OVK_DEBUG_ASSERT(GridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
//   OVK_DEBUG_ASSERT(GridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);

//   auto RowIter = Domain.ExchangeInfo_.find(DonorGridID);
//   if (RowIter != Domain.ExchangeInfo_.end()) {
//     const std::map<int, exchange_info> &Row = RowIter->second;
//     return Row.find(ReceiverGridID) != Row.end();
//   } else {
//     return false;
//   }

// }

// void GetExchangeInfo(const domain &Domain, int DonorGridID, int ReceiverGridID, const exchange_info
//   *&ExchangeInfo) {

//   OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
//   OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
//   OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
//     "exist.", DonorGridID, ReceiverGridID);

//   ExchangeInfo = &Domain.ExchangeInfo_.at(DonorGridID).at(ReceiverGridID);

// }

// bool RankHasExchange(const domain &Domain, int DonorGridID, int ReceiverGridID) {

//   OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
//   OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
//   OVK_DEBUG_ASSERT(GridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
//   OVK_DEBUG_ASSERT(GridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
//   OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
//     "exist.", DonorGridID, ReceiverGridID);

//   auto RowIter = Domain.LocalExchanges_.find(DonorGridID);
//   if (RowIter != Domain.LocalExchanges_.end()) {
//     const std::map<int, exchange> &Row = RowIter->second;
//     return Row.find(ReceiverGridID) != Row.end();
//   } else {
//     return false;
//   }

// }

// void GetExchange(const domain &Domain, int DonorGridID, int ReceiverGridID, const exchange
//   *&Exchange) {

//   OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
//   OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
//   OVK_DEBUG_ASSERT(ExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) does not "
//     "exist.", DonorGridID, ReceiverGridID);
//   if (OVK_DEBUG) {
//     const exchange_info &ExchangeInfo = Domain.ExchangeInfo_.at(DonorGridID).at(ReceiverGridID);
//     std::string Name;
//     GetExchangeInfoName(ExchangeInfo, Name);
//     OVK_DEBUG_ASSERT(RankHasExchange(Domain, DonorGridID, ReceiverGridID), "Exchange %s does not "
//       "have local data on rank @rank@.", Name);
//   }

//   Exchange = &Domain.LocalExchanges_.at(DonorGridID).at(ReceiverGridID);

// }

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

  Domain.Logger_->LogStatus(IsDomainRoot, 0, "Beginning overset assembly on domain %s.",
    Domain.Name_);

//   bool HasOverlap = (Domain.Config_ & domain_config::OVERLAP) != domain_config::NONE;
  bool HasConnectivity = (Domain.Config_ & domain_config::CONNECTIVITY) != domain_config::NONE;
  bool HasExchange = (Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE;

//   if (HasOverlap) {
// //     AssembleOverlap(Domain);
//   }

  if (HasConnectivity) {
    AssembleConnectivity(Domain);
  }

  if (HasExchange) {
    AssembleExchange(Domain);
  }

  if (HasConnectivity) {
    ResetAllConnectivityEdits(Domain);
  }

  Domain.Logger_->LogStatus(IsDomainRoot, 0, "Finished overset assembly on domain %s.",
    Domain.Name_);

}

namespace {

void AssembleConnectivity(domain &Domain) {

  MPI_Barrier(Domain.Comm_);

  // TODO: Try to find a way to do this that doesn't block for each grid pair (serializes work
  // that could be done in parallel)
  for (auto &MPair : Domain.LocalConnectivities_) {
    for (auto &NPair : MPair.second) {
      connectivity &Connectivity = NPair.second;
      UpdateSourceDestRanks(Connectivity);
    }
  }

  MPI_Barrier(Domain.Comm_);

}

void UpdateSourceDestRanks(connectivity &Connectivity) {

  int NumDims;
  GetConnectivityDimension(Connectivity, NumDims);

  core::comm_view Comm = core::GetConnectivityComm(Connectivity);

  bool DonorGridIsLocal = RankHasConnectivityDonorSide(Connectivity);
  bool ReceiverGridIsLocal = RankHasConnectivityReceiverSide(Connectivity);

  const connectivity_d *Donors;
  const grid *DonorGrid;
  long long NumDonors;
  if (DonorGridIsLocal) {
    GetConnectivityDonorSide(Connectivity, Donors);
    GetConnectivityDonorSideGrid(*Donors, DonorGrid);
    GetConnectivityDonorSideCount(*Donors, NumDonors);
  }

  const connectivity_r *Receivers;
  const grid *ReceiverGrid;
  long long NumReceivers;
  if (ReceiverGridIsLocal) {
    GetConnectivityReceiverSide(Connectivity, Receivers);
    GetConnectivityReceiverSideGrid(*Receivers, ReceiverGrid);
    GetConnectivityReceiverSideCount(*Receivers, NumReceivers);
  }

  const connectivity::edits *Edits;
  core::GetConnectivityEdits(Connectivity, Edits);

  if (Edits->DonorDestinations_) {

    array<bool> DonorCommunicates;
    if (DonorGridIsLocal) {
      const cart &Cart = DonorGrid->Cart();
      const range &LocalRange = DonorGrid->LocalRange();
      DonorCommunicates.Resize({NumDonors});
      for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
        tuple<int> DonorLower = Cart.PeriodicAdjust({
          Donors->Extents_(0,0,iDonor),
          Donors->Extents_(0,1,iDonor),
          Donors->Extents_(0,2,iDonor)
        });
        DonorCommunicates(iDonor) = LocalRange.Contains(DonorLower);
      }
    }

    long long NumUnmapped = 0;
    if (DonorGridIsLocal) {
      for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
        if (DonorCommunicates(iDonor) && Donors->DestinationRanks_(iDonor) < 0) {
          ++NumUnmapped;
        }
      }
    }

    int GenerateDestinations = NumUnmapped > 0;
    MPI_Allreduce(MPI_IN_PLACE, &GenerateDestinations, 1, MPI_INT, MPI_MAX, Comm);

    if (GenerateDestinations) {

      const grid_info *ReceiverGridInfo;
      GetConnectivityReceiverGridInfo(Connectivity, ReceiverGridInfo);

      range ReceiverGridGlobalRange;
      GetGridInfoGlobalRange(*ReceiverGridInfo, ReceiverGridGlobalRange);

      range ReceiverGridLocalRange = MakeEmptyRange(NumDims);
      if (ReceiverGridIsLocal) {
        ReceiverGridLocalRange = ReceiverGrid->LocalRange();
      }

      core::partition_hash DestinationHash(NumDims, Comm, ReceiverGridGlobalRange,
        ReceiverGridLocalRange);

      array<int,2> Destinations;
      array<int> DestinationBinIndices;
      std::map<int, core::partition_hash::bin> Bins;

      if (DonorGridIsLocal) {
        Destinations.Resize({{MAX_DIMS,NumUnmapped}});
        DestinationBinIndices.Resize({NumUnmapped});
        long long iUnmapped = 0;
        for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
          if (DonorCommunicates(iDonor) && Donors->DestinationRanks_(iDonor) < 0) {
            for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
              Destinations(iDim,iUnmapped) = Donors->Destinations_(iDim,iDonor);
            }
            ++iUnmapped;
          }
        }
        DestinationHash.MapToBins(Destinations, DestinationBinIndices);
        for (long long iUnmapped = 0; iUnmapped < NumUnmapped; ++iUnmapped) {
          int BinIndex = DestinationBinIndices(iUnmapped);
          auto Iter = Bins.lower_bound(BinIndex);
          if (Iter == Bins.end() || Iter->first > BinIndex) {
            Bins.emplace_hint(Iter, BinIndex, core::partition_hash::bin());
          }
        }
      }

      DestinationHash.RetrieveBins(Bins);

      if (DonorGridIsLocal) {
        array<int> DestinationRanks({NumUnmapped});
        DestinationHash.FindPartitions(Bins, Destinations, DestinationBinIndices,
          DestinationRanks);
        connectivity_d *DonorsEdit;
        EditConnectivityDonorSideLocal(Connectivity, DonorsEdit);
        int *DestinationRanksEdit;
        EditDonorDestinationRanks(*DonorsEdit, DestinationRanksEdit);
        long long iUnmapped = 0;
        for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
          if (DonorCommunicates(iDonor) && Donors->DestinationRanks_(iDonor) < 0) {
            DestinationRanksEdit[iDonor] = DestinationRanks(iUnmapped);
            ++iUnmapped;
          }
        }
        ReleaseDonorDestinationRanks(*DonorsEdit, DestinationRanksEdit);
        ReleaseConnectivityDonorSideLocal(Connectivity, DonorsEdit);
      } else {
        EditConnectivityDonorSideRemote(Connectivity);
        ReleaseConnectivityDonorSideRemote(Connectivity);
      }

    }

  }

  if (Edits->ReceiverSources_) {

    long long NumUnmapped = 0;
    if (ReceiverGridIsLocal) {
      for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
        if (Receivers->SourceRanks_(iReceiver) < 0) {
          ++NumUnmapped;
        }
      }
    }

    int GenerateSources = NumUnmapped > 0;
    MPI_Allreduce(MPI_IN_PLACE, &GenerateSources, 1, MPI_INT, MPI_MAX, Comm);

    if (GenerateSources) {

      const grid_info *DonorGridInfo;
      GetConnectivityDonorGridInfo(Connectivity, DonorGridInfo);

      range DonorGridGlobalRange;
      GetGridInfoGlobalRange(*DonorGridInfo, DonorGridGlobalRange);

      range DonorGridLocalRange = MakeEmptyRange(NumDims);
      if (DonorGridIsLocal) {
        DonorGridLocalRange = DonorGrid->LocalRange();
      }

      core::partition_hash SourceHash(NumDims, Comm, DonorGridGlobalRange, DonorGridLocalRange);

      array<int,2> Sources;
      array<int> SourceBinIndices;
      std::map<int, core::partition_hash::bin> Bins;

      if (ReceiverGridIsLocal) {
        Sources.Resize({{MAX_DIMS,NumUnmapped}});
        SourceBinIndices.Resize({NumUnmapped});
        long long iUnmapped = 0;
        for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
          if (Receivers->SourceRanks_(iReceiver) < 0) {
            for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
              Sources(iDim,iUnmapped) = Receivers->Sources_(iDim,iReceiver);
            }
            ++iUnmapped;
          }
        }
        SourceHash.MapToBins(Sources, SourceBinIndices);
        for (long long iUnmapped = 0; iUnmapped < NumUnmapped; ++iUnmapped) {
          int BinIndex = SourceBinIndices(iUnmapped);
          auto Iter = Bins.lower_bound(BinIndex);
          if (Iter == Bins.end() || Iter->first > BinIndex) {
            Bins.emplace_hint(Iter, BinIndex, core::partition_hash::bin());
          }
        }
      }

      SourceHash.RetrieveBins(Bins);

      if (ReceiverGridIsLocal) {
        array<int> SourceRanks({NumUnmapped});
        SourceHash.FindPartitions(Bins, Sources, SourceBinIndices,
          SourceRanks);
        connectivity_r *ReceiversEdit;
        EditConnectivityReceiverSideLocal(Connectivity, ReceiversEdit);
        int *SourceRanksEdit;
        EditReceiverSourceRanks(*ReceiversEdit, SourceRanksEdit);
        long long iUnmapped = 0;
        for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
          if (Receivers->SourceRanks_(iReceiver) < 0) {
            SourceRanksEdit[iReceiver] = SourceRanks(iUnmapped);
            ++iUnmapped;
          }
        }
        ReleaseReceiverSourceRanks(*ReceiversEdit, SourceRanksEdit);
        ReleaseConnectivityReceiverSideLocal(Connectivity, ReceiversEdit);
      } else {
        EditConnectivityReceiverSideRemote(Connectivity);
        ReleaseConnectivityReceiverSideRemote(Connectivity);
      }

    }

  }

}

void AssembleExchange(domain &Domain) {

  MPI_Barrier(Domain.Comm_);

//   for (auto &MPair : Domain.LocalExchanges_) {
//     for (auto &NPair : MPair.second) {
//       exchange &Exchange = NPair.second;
//       core::UpdateExchange(Exchange);
//     }
//   }

  int CollectTime = core::GetProfilerTimerID(Domain.Profiler_, "Collect");
  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  int DisperseTime = core::GetProfilerTimerID(Domain.Profiler_, "Disperse");

  for (auto &MPair : Domain.LocalConnectivities_) {
    int DonorGridID = MPair.first;
    for (auto &NPair : MPair.second) {
      int ReceiverGridID = NPair.first;
      const connectivity &Connectivity = NPair.second;
      if (RankHasConnectivityDonorSide(Connectivity)) {
        const connectivity::edits *Edits;
        core::GetConnectivityEdits(Connectivity, Edits);
        core::StartProfile(Domain.Profiler_, CollectTime);
        auto CollectDataRowIter = Domain.CollectData_.find(DonorGridID);
        if (CollectDataRowIter != Domain.CollectData_.end()) {
          std::map<int, domain::collect_data> &CollectDataRow = CollectDataRowIter->second;
          CollectDataRow.erase(ReceiverGridID);
          if (CollectDataRow.empty()) {
            Domain.CollectData_.erase(CollectDataRowIter);
          }
        }
        core::EndProfile(Domain.Profiler_, CollectTime);
        core::StartProfile(Domain.Profiler_, SendRecvTime);
        auto SendDataRowIter = Domain.SendData_.find(DonorGridID);
        if (SendDataRowIter != Domain.SendData_.end()) {
          std::map<int, domain::send_data> &SendDataRow = SendDataRowIter->second;
          SendDataRow.erase(ReceiverGridID);
          if (SendDataRow.empty()) {
            Domain.SendData_.erase(SendDataRowIter);
          }
        }
        auto RecvDataRowIter = Domain.RecvData_.find(DonorGridID);
        if (RecvDataRowIter != Domain.RecvData_.end()) {
          std::map<int, domain::recv_data> &RecvDataRow = RecvDataRowIter->second;
          RecvDataRow.erase(ReceiverGridID);
          if (RecvDataRow.empty()) {
            Domain.RecvData_.erase(RecvDataRowIter);
          }
        }
        core::EndProfile(Domain.Profiler_, SendRecvTime);
        core::StartProfile(Domain.Profiler_, DisperseTime);
        auto DisperseDataRowIter = Domain.DisperseData_.find(DonorGridID);
        if (DisperseDataRowIter != Domain.DisperseData_.end()) {
          std::map<int, domain::disperse_data> &DisperseDataRow = DisperseDataRowIter->second;
          DisperseDataRow.erase(ReceiverGridID);
          if (DisperseDataRow.empty()) {
            Domain.DisperseData_.erase(DisperseDataRowIter);
          }
        }
        core::EndProfile(Domain.Profiler_, DisperseTime);
      }
    }
  }

  MPI_Barrier(Domain.Comm_);

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

void CreateCollect(domain &Domain, int DonorGridID, int ReceiverGridID, int CollectID,
  collect_op CollectOp, data_type ValueType, int Count, const range &GridValuesRange, array_layout
  GridValuesLayout) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");
  OVK_DEBUG_ASSERT(ValidCollectOp(CollectOp), "Invalid collect operation.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(ValidArrayLayout(GridValuesLayout), "Invalid grid values layout.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  const grid &DonorGrid = Domain.LocalGrids_.at(DonorGridID);

  OVK_DEBUG_ASSERT(GridValuesRange.Includes(DonorGrid.LocalRange()), "Invalid grid values range.");

  const core::comm &GridComm = DonorGrid.core_Comm();

  MPI_Barrier(GridComm);

  int CollectTime = core::GetProfilerTimerID(Domain.Profiler_, "Collect");
  core::StartProfile(Domain.Profiler_, CollectTime);

  const connectivity &Connectivity = Domain.LocalConnectivities_.at(DonorGridID).at(ReceiverGridID);

  const connectivity_d *Donors;
  GetConnectivityDonorSide(Connectivity, Donors);

  std::map<int, domain::collect_data> &CollectDataRow = Domain.CollectData_[DonorGridID];
  auto CollectDataIter = CollectDataRow.lower_bound(ReceiverGridID);
  if (CollectDataIter == CollectDataRow.end() || CollectDataIter->first > ReceiverGridID) {
    CollectDataIter = CollectDataRow.emplace_hint(CollectDataIter, ReceiverGridID,
      domain::collect_data());
    CollectDataIter->second.Map = core::collect_map(DonorGrid.Cart(), DonorGrid.core_Partition(),
      Donors->Extents_);
  }
  domain::collect_data &CollectData = CollectDataIter->second;
  const core::collect_map &CollectMap = CollectData.Map;
  std::map<int, core::collect> &Collects = CollectData.Collects;

  auto Iter = Collects.lower_bound(CollectID);

  OVK_DEBUG_ASSERT(Iter == Collects.end() || Iter->first > CollectID, "Collect %i already "
    "exists.", CollectID);

  core::collect Collect;

  const cart &Cart = DonorGrid.Cart();
  const range &LocalRange = DonorGrid.LocalRange();

  switch (CollectOp) {
  case collect_op::NONE:
    Collect = core::MakeCollectNone(GridComm, Cart, LocalRange, CollectMap, ValueType, Count,
      GridValuesRange, GridValuesLayout, Domain.Profiler_);
    break;
  case collect_op::ANY:
    Collect = core::MakeCollectAny(GridComm, Cart, LocalRange, CollectMap, ValueType, Count,
      GridValuesRange, GridValuesLayout, Domain.Profiler_);
    break;
  case collect_op::NOT_ALL:
    Collect = core::MakeCollectNotAll(GridComm, Cart, LocalRange, CollectMap, ValueType, Count,
      GridValuesRange, GridValuesLayout, Domain.Profiler_);
    break;
  case collect_op::ALL:
    Collect = core::MakeCollectAll(GridComm, Cart, LocalRange, CollectMap, ValueType, Count,
      GridValuesRange, GridValuesLayout, Domain.Profiler_);
    break;
  case collect_op::INTERPOLATE:
    Collect = core::MakeCollectInterp(GridComm, Cart, LocalRange, CollectMap, ValueType, Count,
      GridValuesRange, GridValuesLayout, Domain.Profiler_, Donors->InterpCoefs_);
    break;
  }

  Collects.emplace_hint(Iter, CollectID, std::move(Collect));

  core::EndProfile(Domain.Profiler_, CollectTime);

  MPI_Barrier(GridComm);

}

void DestroyCollect(domain &Domain, int DonorGridID, int ReceiverGridID, int CollectID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  const grid &DonorGrid = Domain.LocalGrids_.at(DonorGridID);

  const core::comm &Comm = DonorGrid.core_Comm();

  MPI_Barrier(Comm);

  int CollectTime = core::GetProfilerTimerID(Domain.Profiler_, "Collect");
  core::StartProfile(Domain.Profiler_, CollectTime);

  auto CollectDataRowIter = Domain.CollectData_.find(DonorGridID);

  OVK_DEBUG_ASSERT(CollectDataRowIter != Domain.CollectData_.end(), "Collect %i does not exist.",
    CollectID);

  std::map<int, domain::collect_data> &CollectDataRow = CollectDataRowIter->second;

  auto CollectDataIter = CollectDataRow.find(ReceiverGridID);

  OVK_DEBUG_ASSERT(CollectDataIter != CollectDataRow.end(), "Collect %i does not exist.",
    CollectID);

  domain::collect_data &CollectData = CollectDataIter->second;
  std::map<int, core::collect> &Collects = CollectData.Collects;

  auto Iter = Collects.find(CollectID);

  OVK_DEBUG_ASSERT(Iter != Collects.end(), "Collect %i does not exist.", CollectID);

  Collects.erase(Iter);

  if (Collects.empty()) {
    CollectDataRow.erase(CollectDataIter);
    if (CollectDataRow.empty()) {
      Domain.CollectData_.erase(CollectDataRowIter);
    }
  }

  core::EndProfile(Domain.Profiler_, CollectTime);

  MPI_Barrier(Comm);

}

namespace {

const std::map<int, core::collect> *FindCollects(const domain &Domain, int DonorGridID, int
  ReceiverGridID) {

  const std::map<int, core::collect> *Collects = nullptr;

  auto CollectDataRowIter = Domain.CollectData_.find(DonorGridID);
  if (CollectDataRowIter != Domain.CollectData_.end()) {
    const std::map<int, domain::collect_data> &CollectDataRow = CollectDataRowIter->second;
    auto CollectDataIter = CollectDataRow.find(ReceiverGridID);
    if (CollectDataIter != CollectDataRow.end()) {
      const domain::collect_data &CollectData = CollectDataIter->second;
      Collects = &CollectData.Collects;
    }
  }

  return Collects;

}

std::map<int, core::collect> *FindCollects(domain &Domain, int DonorGridID, int ReceiverGridID) {

  std::map<int, core::collect> *Collects = nullptr;

  auto CollectDataRowIter = Domain.CollectData_.find(DonorGridID);
  if (CollectDataRowIter != Domain.CollectData_.end()) {
    std::map<int, domain::collect_data> &CollectDataRow = CollectDataRowIter->second;
    auto CollectDataIter = CollectDataRow.find(ReceiverGridID);
    if (CollectDataIter != CollectDataRow.end()) {
      domain::collect_data &CollectData = CollectDataIter->second;
      Collects = &CollectData.Collects;
    }
  }

  return Collects;

}

}

void Collect(domain &Domain, int DonorGridID, int ReceiverGridID, int CollectID, const void * const
  *GridValues, void **DonorValues) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  int CollectTime = core::GetProfilerTimerID(Domain.Profiler_, "Collect");
  core::StartProfile(Domain.Profiler_, CollectTime);

  std::map<int, core::collect> *CollectsPtr = FindCollects(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(CollectsPtr, "Collect %i does not exist.", CollectID);
  auto Iter = CollectsPtr->find(CollectID);
  OVK_DEBUG_ASSERT(Iter != CollectsPtr->end(), "Collect %i does not exist.", CollectID);
  core::collect &Collect = Iter->second;

  Collect.Collect(GridValues, DonorValues);

  core::EndProfile(Domain.Profiler_, CollectTime);

}

bool CollectExists(const domain &Domain, int DonorGridID, int ReceiverGridID, int CollectID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(CollectID >= 0, "Invalid collect ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  bool Exists = false;

  const std::map<int, core::collect> *CollectsPtr = FindCollects(Domain, DonorGridID,
    ReceiverGridID);

  if (CollectsPtr) {
    auto Iter = CollectsPtr->find(CollectID);
    if (Iter != CollectsPtr->end()) {
      Exists = true;
    }
  }

  return Exists;

}

void GetNextAvailableCollectID(const domain &Domain, int DonorGridID, int ReceiverGridID, int
  &CollectID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  const std::map<int, core::collect> *CollectsPtr = FindCollects(Domain, DonorGridID,
    ReceiverGridID);

  if (CollectsPtr) {
    CollectID = GetNextAvailableID(*CollectsPtr);
  } else {
    CollectID = 0;
  }

}

void CreateSend(domain &Domain, int DonorGridID, int ReceiverGridID, int SendID, data_type
  ValueType, int Count, int Tag) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(Tag >= 0, "Invalid tag.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  const connectivity &Connectivity = Domain.LocalConnectivities_.at(DonorGridID).at(ReceiverGridID);

  const core::comm &ConnectivityComm = core::GetConnectivityComm(Connectivity);

  const connectivity_d *Donors;
  GetConnectivityDonorSide(Connectivity, Donors);

  std::map<int, domain::send_data> &SendDataRow = Domain.SendData_[DonorGridID];
  auto SendDataIter = SendDataRow.lower_bound(ReceiverGridID);
  if (SendDataIter == SendDataRow.end() || SendDataIter->first > ReceiverGridID) {
    const grid_info *ReceiverGridInfo;
    GetConnectivityReceiverGridInfo(Connectivity, ReceiverGridInfo);
    range ReceiverGridGlobalRange;
    GetGridInfoGlobalRange(*ReceiverGridInfo, ReceiverGridGlobalRange);
    array<long long> Order = GetSendRecvOrder(Donors->Destinations_, ReceiverGridGlobalRange);
    SendDataIter = SendDataRow.emplace_hint(SendDataIter, ReceiverGridID, domain::send_data());
    SendDataIter->second.Map = core::send_map(Donors->Count_, Order, Donors->DestinationRanks_);
  }
  domain::send_data &SendData = SendDataIter->second;
  const core::send_map &SendMap = SendData.Map;
  std::map<int, core::send> &Sends = SendData.Sends;

  auto Iter = Sends.lower_bound(SendID);

  OVK_DEBUG_ASSERT(Iter == Sends.end() || Iter->first > SendID, "Send %i already exists.", SendID);

  core::send Send = core::MakeSend(ConnectivityComm, SendMap, ValueType, Count, Tag,
    Domain.Profiler_);

  Sends.emplace_hint(Iter, SendID, std::move(Send));

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

void DestroySend(domain &Domain, int DonorGridID, int ReceiverGridID, int SendID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  auto SendDataRowIter = Domain.SendData_.find(DonorGridID);

  OVK_DEBUG_ASSERT(SendDataRowIter != Domain.SendData_.end(), "Send %i does not exist.", SendID);

  std::map<int, domain::send_data> &SendDataRow = SendDataRowIter->second;

  auto SendDataIter = SendDataRow.find(ReceiverGridID);

  OVK_DEBUG_ASSERT(SendDataIter != SendDataRow.end(), "Send %i does not exist.", SendID);

  domain::send_data &SendData = SendDataIter->second;
  std::map<int, core::send> &Sends = SendData.Sends;

  auto Iter = Sends.find(SendID);

  OVK_DEBUG_ASSERT(Iter != Sends.end(), "Send %i does not exist.", SendID);

  Sends.erase(Iter);

  if (Sends.empty()) {
    SendDataRow.erase(SendDataIter);
    if (SendDataRow.empty()) {
      Domain.SendData_.erase(SendDataRowIter);
    }
  }

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

namespace {

const std::map<int, core::send> *FindSends(const domain &Domain, int DonorGridID, int
  ReceiverGridID) {

  const std::map<int, core::send> *Sends = nullptr;

  auto SendDataRowIter = Domain.SendData_.find(DonorGridID);
  if (SendDataRowIter != Domain.SendData_.end()) {
    const std::map<int, domain::send_data> &SendDataRow = SendDataRowIter->second;
    auto SendDataIter = SendDataRow.find(ReceiverGridID);
    if (SendDataIter != SendDataRow.end()) {
      const domain::send_data &SendData = SendDataIter->second;
      Sends = &SendData.Sends;
    }
  }

  return Sends;

}

std::map<int, core::send> *FindSends(domain &Domain, int DonorGridID, int ReceiverGridID) {

  std::map<int, core::send> *Sends = nullptr;

  auto SendDataRowIter = Domain.SendData_.find(DonorGridID);
  if (SendDataRowIter != Domain.SendData_.end()) {
    std::map<int, domain::send_data> &SendDataRow = SendDataRowIter->second;
    auto SendDataIter = SendDataRow.find(ReceiverGridID);
    if (SendDataIter != SendDataRow.end()) {
      domain::send_data &SendData = SendDataIter->second;
      Sends = &SendData.Sends;
    }
  }

  return Sends;

}

}

request Send(domain &Domain, int DonorGridID, int ReceiverGridID, int SendID, const void * const
  *DonorValues) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  std::map<int, core::send> *SendsPtr = FindSends(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(SendsPtr, "Send %i does not exist.", SendID);
  auto Iter = SendsPtr->find(SendID);
  OVK_DEBUG_ASSERT(Iter != SendsPtr->end(), "Send %i does not exist.", SendID);
  core::send &Send = Iter->second;

  request Request = Send.Send(DonorValues);

  core::EndProfile(Domain.Profiler_, SendRecvTime);

  return Request;

}

bool SendExists(const domain &Domain, int DonorGridID, int ReceiverGridID, int SendID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(SendID >= 0, "Invalid send ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  bool Exists = false;

  const std::map<int, core::send> *SendsPtr = FindSends(Domain, DonorGridID, ReceiverGridID);

  if (SendsPtr) {
    auto Iter = SendsPtr->find(SendID);
    if (Iter != SendsPtr->end()) {
      Exists = true;
    }
  }

  return Exists;

}

void GetNextAvailableSendID(const domain &Domain, int DonorGridID, int ReceiverGridID, int
  &SendID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(DonorGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, DonorGridID), "Grid %s does not have local data on rank "
      "@rank@.", GridName);
  }

  const std::map<int, core::send> *SendsPtr = FindSends(Domain, DonorGridID, ReceiverGridID);

  if (SendsPtr) {
    SendID = GetNextAvailableID(*SendsPtr);
  } else {
    SendID = 0;
  }

}

void CreateReceive(domain &Domain, int DonorGridID, int ReceiverGridID, int RecvID, data_type
  ValueType, int Count, int Tag) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(Tag >= 0, "Invalid tag.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  const connectivity &Connectivity = Domain.LocalConnectivities_.at(DonorGridID).at(ReceiverGridID);

  const core::comm &ConnectivityComm = core::GetConnectivityComm(Connectivity);

  const connectivity_r *Receivers;
  GetConnectivityReceiverSide(Connectivity, Receivers);

  std::map<int, domain::recv_data> &RecvDataRow = Domain.RecvData_[DonorGridID];
  auto RecvDataIter = RecvDataRow.lower_bound(ReceiverGridID);
  if (RecvDataIter == RecvDataRow.end() || RecvDataIter->first > ReceiverGridID) {
    const grid *ReceiverGrid;
    GetConnectivityReceiverSideGrid(*Receivers, ReceiverGrid);
    array<long long> Order = GetSendRecvOrder(Receivers->Points_, ReceiverGrid->GlobalRange());
    RecvDataIter = RecvDataRow.emplace_hint(RecvDataIter, ReceiverGridID, domain::recv_data());
    RecvDataIter->second.Map = core::recv_map(Receivers->Count_, Order, Receivers->SourceRanks_);
  }
  domain::recv_data &RecvData = RecvDataIter->second;
  const core::recv_map &RecvMap = RecvData.Map;
  std::map<int, core::recv> &Recvs = RecvData.Recvs;

  auto Iter = Recvs.lower_bound(RecvID);

  OVK_DEBUG_ASSERT(Iter == Recvs.end() || Iter->first > RecvID, "Receive %i already exists.",
    RecvID);

  core::recv Recv = core::MakeRecv(ConnectivityComm, RecvMap, ValueType, Count, Tag,
    Domain.Profiler_);

  Recvs.emplace_hint(Iter, RecvID, std::move(Recv));

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

void DestroyReceive(domain &Domain, int DonorGridID, int ReceiverGridID, int RecvID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  auto RecvDataRowIter = Domain.RecvData_.find(DonorGridID);

  OVK_DEBUG_ASSERT(RecvDataRowIter != Domain.RecvData_.end(), "Receive %i does not exist.", RecvID);

  std::map<int, domain::recv_data> &RecvDataRow = RecvDataRowIter->second;

  auto RecvDataIter = RecvDataRow.find(ReceiverGridID);

  OVK_DEBUG_ASSERT(RecvDataIter != RecvDataRow.end(), "Receive %i does not exist.", RecvID);

  domain::recv_data &RecvData = RecvDataIter->second;
  std::map<int, core::recv> &Recvs = RecvData.Recvs;

  auto Iter = Recvs.find(RecvID);

  OVK_DEBUG_ASSERT(Iter != Recvs.end(), "Receive %i does not exist.", RecvID);

  Recvs.erase(Iter);

  if (Recvs.empty()) {
    RecvDataRow.erase(RecvDataIter);
    if (RecvDataRow.empty()) {
      Domain.RecvData_.erase(RecvDataRowIter);
    }
  }

  core::EndProfile(Domain.Profiler_, SendRecvTime);

}

namespace {

const std::map<int, core::recv> *FindReceives(const domain &Domain, int DonorGridID, int
  ReceiverGridID) {

  const std::map<int, core::recv> *Recvs = nullptr;

  auto RecvDataRowIter = Domain.RecvData_.find(DonorGridID);
  if (RecvDataRowIter != Domain.RecvData_.end()) {
    const std::map<int, domain::recv_data> &RecvDataRow = RecvDataRowIter->second;
    auto RecvDataIter = RecvDataRow.find(ReceiverGridID);
    if (RecvDataIter != RecvDataRow.end()) {
      const domain::recv_data &RecvData = RecvDataIter->second;
      Recvs = &RecvData.Recvs;
    }
  }

  return Recvs;

}

std::map<int, core::recv> *FindReceives(domain &Domain, int DonorGridID, int ReceiverGridID) {

  std::map<int, core::recv> *Recvs = nullptr;

  auto RecvDataRowIter = Domain.RecvData_.find(DonorGridID);
  if (RecvDataRowIter != Domain.RecvData_.end()) {
    std::map<int, domain::recv_data> &RecvDataRow = RecvDataRowIter->second;
    auto RecvDataIter = RecvDataRow.find(ReceiverGridID);
    if (RecvDataIter != RecvDataRow.end()) {
      domain::recv_data &RecvData = RecvDataIter->second;
      Recvs = &RecvData.Recvs;
    }
  }

  return Recvs;

}

}

request Receive(domain &Domain, int DonorGridID, int ReceiverGridID, int RecvID, void
  **ReceiverValues) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  core::StartProfile(Domain.Profiler_, SendRecvTime);

  std::map<int, core::recv> *RecvsPtr = FindReceives(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(RecvsPtr, "Receive %i does not exist.", RecvID);
  auto Iter = RecvsPtr->find(RecvID);
  OVK_DEBUG_ASSERT(Iter != RecvsPtr->end(), "Receive %i does not exist.", RecvID);
  core::recv &Recv = Iter->second;

  request Request = Recv.Recv(ReceiverValues);

  core::EndProfile(Domain.Profiler_, SendRecvTime);

  return Request;

}

bool ReceiveExists(const domain &Domain, int DonorGridID, int ReceiverGridID, int RecvID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(RecvID >= 0, "Invalid receive ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  bool Exists = false;

  const std::map<int, core::recv> *RecvsPtr = FindReceives(Domain, DonorGridID, ReceiverGridID);

  if (RecvsPtr) {
    auto Iter = RecvsPtr->find(RecvID);
    if (Iter != RecvsPtr->end()) {
      Exists = true;
    }
  }

  return Exists;

}

void GetNextAvailableReceiveID(const domain &Domain, int DonorGridID, int ReceiverGridID, int
  &RecvID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  const std::map<int, core::recv> *RecvsPtr = FindReceives(Domain, DonorGridID, ReceiverGridID);

  if (RecvsPtr) {
    RecvID = GetNextAvailableID(*RecvsPtr);
  } else {
    RecvID = 0;
  }

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

void CreateDisperse(domain &Domain, int DonorGridID, int ReceiverGridID, int DisperseID,
  disperse_op DisperseOp, data_type ValueType, int Count, const range &GridValuesRange, array_layout
  GridValuesLayout) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");
  OVK_DEBUG_ASSERT(ValidDisperseOp(DisperseOp), "Invalid disperse operation.");
  OVK_DEBUG_ASSERT(ValidDataType(ValueType), "Invalid value type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(ValidArrayLayout(GridValuesLayout), "Invalid grid values layout.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  const grid &ReceiverGrid = Domain.LocalGrids_.at(ReceiverGridID);

  OVK_DEBUG_ASSERT(GridValuesRange.Includes(ReceiverGrid.LocalRange()), "Invalid grid values "
    "range.");

  int DisperseTime = core::GetProfilerTimerID(Domain.Profiler_, "Disperse");
  core::StartProfile(Domain.Profiler_, DisperseTime);

  const connectivity &Connectivity = Domain.LocalConnectivities_.at(DonorGridID).at(ReceiverGridID);

  const connectivity_r *Receivers;
  GetConnectivityReceiverSide(Connectivity, Receivers);

  std::map<int, domain::disperse_data> &DisperseDataRow = Domain.DisperseData_[DonorGridID];
  auto DisperseDataIter = DisperseDataRow.lower_bound(ReceiverGridID);
  if (DisperseDataIter == DisperseDataRow.end() || DisperseDataIter->first > ReceiverGridID) {
    DisperseDataIter = DisperseDataRow.emplace_hint(DisperseDataIter, ReceiverGridID,
      domain::disperse_data());
  }
  domain::disperse_data &DisperseData = DisperseDataIter->second;
  std::map<int, core::disperse> &Disperses = DisperseData.Disperses;

  auto Iter = Disperses.lower_bound(DisperseID);

  OVK_DEBUG_ASSERT(Iter == Disperses.end() || Iter->first > DisperseID, "Disperse %i already "
    "exists.", DisperseID);

  core::disperse Disperse;

  switch (DisperseOp) {
  case disperse_op::OVERWRITE:
    Disperse = core::MakeDisperseOverwrite(Receivers->Points_, ValueType, Count, GridValuesRange,
      GridValuesLayout, Domain.Profiler_);
    break;
  case disperse_op::APPEND:
    Disperse = core::MakeDisperseAppend(Receivers->Points_, ValueType, Count, GridValuesRange,
      GridValuesLayout, Domain.Profiler_);
    break;
  }

  Disperses.emplace_hint(Iter, DisperseID, std::move(Disperse));

  core::EndProfile(Domain.Profiler_, DisperseTime);

}

void DestroyDisperse(domain &Domain, int DonorGridID, int ReceiverGridID, int DisperseID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  int DisperseTime = core::GetProfilerTimerID(Domain.Profiler_, "Disperse");
  core::StartProfile(Domain.Profiler_, DisperseTime);

  auto DisperseDataRowIter = Domain.DisperseData_.find(DonorGridID);

  OVK_DEBUG_ASSERT(DisperseDataRowIter != Domain.DisperseData_.end(), "Disperse %i does not exist.",
    DisperseID);

  std::map<int, domain::disperse_data> &DisperseDataRow = DisperseDataRowIter->second;

  auto DisperseDataIter = DisperseDataRow.find(ReceiverGridID);

  OVK_DEBUG_ASSERT(DisperseDataIter != DisperseDataRow.end(), "Disperse %i does not exist.",
    DisperseID);

  domain::disperse_data &DisperseData = DisperseDataIter->second;
  std::map<int, core::disperse> &Disperses = DisperseData.Disperses;

  auto Iter = Disperses.find(DisperseID);

  OVK_DEBUG_ASSERT(Iter != Disperses.end(), "Disperse %i does not exist.", DisperseID);

  Disperses.erase(Iter);

  if (Disperses.empty()) {
    DisperseDataRow.erase(DisperseDataIter);
    if (DisperseDataRow.empty()) {
      Domain.DisperseData_.erase(DisperseDataRowIter);
    }
  }

  core::EndProfile(Domain.Profiler_, DisperseTime);

}

namespace {

const std::map<int, core::disperse> *FindDisperses(const domain &Domain, int DonorGridID, int
  ReceiverGridID) {

  const std::map<int, core::disperse> *Disperses = nullptr;

  auto DisperseDataRowIter = Domain.DisperseData_.find(DonorGridID);
  if (DisperseDataRowIter != Domain.DisperseData_.end()) {
    const std::map<int, domain::disperse_data> &DisperseDataRow = DisperseDataRowIter->second;
    auto DisperseDataIter = DisperseDataRow.find(ReceiverGridID);
    if (DisperseDataIter != DisperseDataRow.end()) {
      const domain::disperse_data &DisperseData = DisperseDataIter->second;
      Disperses = &DisperseData.Disperses;
    }
  }

  return Disperses;

}

std::map<int, core::disperse> *FindDisperses(domain &Domain, int DonorGridID, int ReceiverGridID) {

  std::map<int, core::disperse> *Disperses = nullptr;

  auto DisperseDataRowIter = Domain.DisperseData_.find(DonorGridID);
  if (DisperseDataRowIter != Domain.DisperseData_.end()) {
    std::map<int, domain::disperse_data> &DisperseDataRow = DisperseDataRowIter->second;
    auto DisperseDataIter = DisperseDataRow.find(ReceiverGridID);
    if (DisperseDataIter != DisperseDataRow.end()) {
      domain::disperse_data &DisperseData = DisperseDataIter->second;
      Disperses = &DisperseData.Disperses;
    }
  }

  return Disperses;

}

}

void Disperse(domain &Domain, int DonorGridID, int ReceiverGridID, int DisperseID, const void *
  const *ReceiverValues, void **GridValues) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  int DisperseTime = core::GetProfilerTimerID(Domain.Profiler_, "Disperse");
  core::StartProfile(Domain.Profiler_, DisperseTime);

  std::map<int, core::disperse> *DispersesPtr = FindDisperses(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(DispersesPtr, "Disperse %i does not exist.", DisperseID);
  auto Iter = DispersesPtr->find(DisperseID);
  OVK_DEBUG_ASSERT(Iter != DispersesPtr->end(), "Disperse %i does not exist.", DisperseID);
  core::disperse &Disperse = Iter->second;

  Disperse.Disperse(ReceiverValues, GridValues);

  core::EndProfile(Domain.Profiler_, DisperseTime);

}

bool DisperseExists(const domain &Domain, int DonorGridID, int ReceiverGridID, int DisperseID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(DisperseID >= 0, "Invalid disperse ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  bool Exists = false;

  const std::map<int, core::disperse> *DispersesPtr = FindDisperses(Domain, DonorGridID,
    ReceiverGridID);

  if (DispersesPtr) {
    auto Iter = DispersesPtr->find(DisperseID);
    if (Iter != DispersesPtr->end()) {
      Exists = true;
    }
  }

  return Exists;

}

void GetNextAvailableDisperseID(const domain &Domain, int DonorGridID, int ReceiverGridID, int
  &DisperseID) {

  OVK_DEBUG_ASSERT((Domain.Config_ & domain_config::EXCHANGE) != domain_config::NONE, "Domain %s "
    "is not configured for exchange.", Domain.Name_);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");

  OVK_DEBUG_ASSERT(ConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);

  if (OVK_DEBUG) {
    const grid_info &GridInfo = Domain.GridInfo_.at(ReceiverGridID);
    std::string GridName;
    GetGridInfoName(GridInfo, GridName);
    OVK_DEBUG_ASSERT(RankHasGrid(Domain, ReceiverGridID), "Grid %s does not have local data on "
      "rank @rank@.", GridName);
  }

  const std::map<int, core::disperse> *DispersesPtr = FindDisperses(Domain, DonorGridID,
    ReceiverGridID);

  if (DispersesPtr) {
    DisperseID = GetNextAvailableID(*DispersesPtr);
  } else {
    DisperseID = 0;
  }

}

namespace {

template <typename T> int GetNextAvailableID(const std::map<int, T> &Map) {

  auto Iter = Map.begin();

  if (Iter == Map.end() || Iter->first > 0) {
    return 0;
  }

  auto PrevIter = Iter;
  ++Iter;
  while (Iter != Map.end()) {
    if (Iter->first - PrevIter->first > 1) {
      return PrevIter->first + 1;
    }
    PrevIter = Iter;
    ++Iter;
  }

  return PrevIter->first + 1;

}

array<long long> GetSendRecvOrder(const array<int,2> &ReceiverPoints, const range
  &ReceiverGridGlobalRange) {

  long long NumReceivers = ReceiverPoints.Size(1);

  array<long long> Order({NumReceivers});

  using range_indexer = indexer<long long, int, MAX_DIMS, array_layout::GRID>;
  range_indexer ReceiverGridGlobalIndexer(ReceiverGridGlobalRange);

  array<long long> ReceiverIndices({NumReceivers});

  for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
    tuple<int> Point = {
      ReceiverPoints(0,iReceiver),
      ReceiverPoints(1,iReceiver),
      ReceiverPoints(2,iReceiver)
    };
    ReceiverIndices(iReceiver) = ReceiverGridGlobalIndexer.ToIndex(Point);
  }

  bool Sorted = true;

  long long PrevIndex = 0;
  for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
    if (ReceiverIndices(iReceiver) < PrevIndex) {
      Sorted = false;
      break;
    }
    PrevIndex = ReceiverIndices(iReceiver);
  }

  if (Sorted) {
    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      Order(iReceiver) = iReceiver;
    }
  } else {
    core::SortPermutation(ReceiverIndices, Order);
  }

  return Order;

}

void CreateDomainGridInfo(domain::grid_info &GridInfo, grid *Grid, core::comm_view Comm) {

  core::CreateGridInfo(GridInfo, Grid, Comm);
  GridInfo.EditRefCount_ = 0;

}

void DestroyDomainGridInfo(domain::grid_info &GridInfo) {

  core::DestroyGridInfo(GridInfo);

}

void CreateDomainConnectivityInfo(domain::connectivity_info &ConnectivityInfo, connectivity
  *Connectivity, core::comm_view Comm) {

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
