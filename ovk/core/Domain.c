// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Domain.h"

#include "ovk/core/AssemblyOptions.h"
#include "ovk/core/Connectivity.h"
#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Exchange.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/Logger.h"
#include "ovk/core/MPIUtils.h"
#include "ovk/core/OrderedMap.h"
#include "ovk/core/TextUtils.h"

static void CreateGridGlobal(t_domain_grid_container **GridContainer, int GridID,
  const ovk_grid_params *Params, MPI_Comm Comm, int CommRank, t_logger *Logger,
  t_error_handler *ErrorHandler);
static void DestroyGridGlobal(t_domain_grid_container **GridContainer);

static t_domain_grid_container *FindGrid(ovk_domain *Domain, int GridID);
static const t_domain_grid_container *FindGridC(const ovk_domain *Domain, int GridID);

static void EditGridGlobal(ovk_domain *Domain, int GridID, ovk_grid **Grid);
static void ReleaseGridGlobal(ovk_domain *Domain, int GridID, ovk_grid **Grid);

static void EnableConnectivityComponent(ovk_domain *Domain);
static void DisableConnectivityComponent(ovk_domain *Domain);

static void CreateConnectivitiesForGrid(ovk_domain *Domain, int GridID);
static void DestroyConnectivitiesForGrid(ovk_domain *Domain, int GridID);

static void CreateConnectivityGlobal(t_domain_connectivity_container **ConnectivityContainer,
  int NumDims, const t_domain_grid_container *DonorGridContainer,
  const t_domain_grid_container *ReceiverGridContainer, MPI_Comm Comm, int CommRank,
  t_logger *Logger, t_error_handler *ErrorHandler);
static void DestroyConnectivityGlobal(t_domain_connectivity_container **ConnectivityContainer);

static t_domain_connectivity_container *FindConnectivity(ovk_domain *Domain, int DonorGridID,
  int ReceiverGridID);
static const t_domain_connectivity_container *FindConnectivityC(const ovk_domain *Domain,
  int DonorGridID, int ReceiverGridID);

static void EditConnectivityGlobal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity);
static void ReleaseConnectivityGlobal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity);

static void EnableExchangeComponent(ovk_domain *Domain);
static void DisableExchangeComponent(ovk_domain *Domain);

static void CreateExchangesForGrid(ovk_domain *Domain, int GridID);
static void DestroyExchangesForGrid(ovk_domain *Domain, int GridID);

static void CreateExchangeGlobal(t_domain_exchange_container **ExchangeContainer,
  const t_domain_connectivity_container *ConnectivityContainer, MPI_Comm Comm, int CommRank,
  t_logger *Logger, t_error_handler *ErrorHandler);
static void DestroyExchangeGlobal(t_domain_exchange_container **ExchangeContainer);

// static t_domain_exchange_container *FindExchange(ovk_domain *Domain, int DonorGridID,
//   int ReceiverGridID);
static const t_domain_exchange_container *FindExchangeC(const ovk_domain *Domain, int DonorGridID,
  int ReceiverGridID);

static bool EditingProperties(const ovk_domain *Domain);
static bool EditingGrid(const ovk_domain *Domain, int GridID);
static bool EditingConnectivity(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID);

static void AssembleExchange(ovk_domain *Domain);

static void ResetAllConnectivityEdits(ovk_domain *Domain);

static void CreateGridContainer(t_domain_grid_container **Container, ovk_grid *Grid, MPI_Comm Comm,
  int CommRank);
static void DestroyGridContainer(t_domain_grid_container **Container);

static void CreateConnectivityContainer(t_domain_connectivity_container **Container,
  ovk_connectivity *Connectivity, MPI_Comm Comm, int CommRank);
static void DestroyConnectivityContainer(t_domain_connectivity_container **Container);

static void CreateExchangeContainer(t_domain_exchange_container **Container, ovk_exchange *Exchange,
  MPI_Comm Comm, int CommRank);
static void DestroyExchangeContainer(t_domain_exchange_container **Container);

static void DefaultProperties(ovk_domain_properties *Properties);

void PRIVATE(CreateDomain)(ovk_domain **Domain_, const ovk_domain_params *Params, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  OVK_DEBUG_ASSERT(Domain_, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  MPI_Comm Comm;
  MPI_Comm_dup(Params->comm, &Comm);

  MPI_Barrier(Comm);

  *Domain_ = malloc(sizeof(ovk_domain));
  ovk_domain *Domain = *Domain_;

  DefaultProperties(&Domain->properties);

  if (strlen(Params->name) > 0) {
    strncpy(Domain->properties.name, Params->name, OVK_NAME_LENGTH);
  } else {
    strcpy(Domain->properties.name, "Domain");
  }

  Domain->properties.num_dims = Params->num_dims;

  Domain->properties.comm = Comm;
  MPI_Comm_size(Domain->properties.comm, &Domain->properties.comm_size);
  MPI_Comm_rank(Domain->properties.comm, &Domain->properties.comm_rank);

  Domain->properties_edit_ref_count = 0;

  Domain->logger = Logger;
  Domain->error_handler = ErrorHandler;

  Domain->config = OVK_DOMAIN_CONFIG_NONE;

  OMCreate(&Domain->grids);
  Domain->grids_edit_ref_count = 0;

  Domain->connectivities_edit_ref_count = 0;

  if (Domain->properties.comm_rank == 0) {
    char ProcessesString[NUMBER_STRING_LENGTH+10];
    PluralizeLabel(Domain->properties.comm_size, "processes", "process", ProcessesString);
    LogStatus(Logger, true, 0, "Created %1iD domain %s on %s.", Domain->properties.num_dims,
      Domain->properties.name, ProcessesString);
  }

  MPI_Barrier(Domain->properties.comm);

}

void PRIVATE(DestroyDomain)(ovk_domain **Domain_) {

  OVK_DEBUG_ASSERT(Domain_, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(*Domain_, "Invalid domain pointer.");

  ovk_domain *Domain = *Domain_;

  MPI_Barrier(Domain->properties.comm);

  if (Domain->config & OVK_DOMAIN_CONFIG_EXCHANGE) {
    DisableExchangeComponent(Domain);
  }

  if (Domain->config & OVK_DOMAIN_CONFIG_CONNECTIVITY) {
    DisableConnectivityComponent(Domain);
  }

  t_ordered_map_entry *Entry = OMBegin(Domain->grids);
  while (Entry != OMEnd(Domain->grids)) {
    t_domain_grid_container *Container = OMData(Entry);
    if (Container->has_local_data) {
      DestroyGrid(&Container->grid);
    }
    DestroyGridContainer(&Container);
    Entry = OMNext(Entry);
  }
  OMDestroy(&Domain->grids);

  t_logger *Logger = Domain->logger;
  MPI_Comm Comm = Domain->properties.comm;
  bool IsRoot = Domain->properties.comm_rank == 0;
  char Name[OVK_NAME_LENGTH];
  strncpy(Name, Domain->properties.name, OVK_NAME_LENGTH);

  free_null(Domain_);

  MPI_Barrier(Comm);

  MPI_Comm_free(&Comm);

  LogStatus(Logger, IsRoot, 0, "Destroyed domain %s.", Name);

}

void ovkConfigureDomain(ovk_domain *Domain, ovk_domain_config Config) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(ValidDomainConfig(Config), "Invalid domain config.");
  OVK_DEBUG_ASSERT(!EditingProperties(Domain), "Cannot configure domain while editing properties.");
  OVK_DEBUG_ASSERT(!EditingGrid(Domain, OVK_ALL_GRIDS), "Cannot configure domain while editing "
    "grids.");
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot configure "
    "domain while editing connectivities.");

  MPI_Barrier(Domain->properties.comm);

  ovk_domain_config OldConfig = Domain->config;

  bool HasGeometry = Config & OVK_DOMAIN_CONFIG_GEOMETRY;
  bool HasOverlap = Config & OVK_DOMAIN_CONFIG_OVERLAP;
  bool HasConnectivity = Config & OVK_DOMAIN_CONFIG_CONNECTIVITY;
  bool HasExchange = Config & OVK_DOMAIN_CONFIG_EXCHANGE;

  OVK_DEBUG_ASSERT(!HasOverlap || HasGeometry, "Domain overlap component requires geometry "
    "component.");
  OVK_DEBUG_ASSERT(!HasExchange || HasConnectivity, "Domain exchange component requires "
    "connectivity component.");

  Domain->config = Config;

  ovk_domain_config ConfigAdded = ~OldConfig & Config;
  ovk_domain_config ConfigRemoved = OldConfig & ~Config;

  if (ConfigAdded & OVK_DOMAIN_CONFIG_CONNECTIVITY) {
    EnableConnectivityComponent(Domain);
  } else if (ConfigRemoved & OVK_DOMAIN_CONFIG_CONNECTIVITY) {
    DisableConnectivityComponent(Domain);
  }

  if (ConfigAdded & OVK_DOMAIN_CONFIG_EXCHANGE) {
    EnableExchangeComponent(Domain);
  } else if (ConfigRemoved & OVK_DOMAIN_CONFIG_EXCHANGE) {
    DisableExchangeComponent(Domain);
  }

  MPI_Barrier(Domain->properties.comm);

}

void ovkGetDomainConfiguration(ovk_domain *Domain, ovk_domain_config *Config) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Config, "Invalid config pointer.");

  *Config = Domain->config;

}

void ovkGetDomainProperties(const ovk_domain *Domain, const ovk_domain_properties **Properties) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");

  *Properties = &Domain->properties;

}

void ovkEditDomainProperties(ovk_domain *Domain, ovk_domain_properties **Properties) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(!EditingGrid(Domain, OVK_ALL_GRIDS), "Cannot edit properties while editing "
    "grids.");
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot edit "
    "grid while editing connectivities.");

  bool StartEdit = Domain->properties_edit_ref_count == 0;
  ++Domain->properties_edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Domain->properties.comm);
  }

  *Properties = &Domain->properties;

}

void ovkReleaseDomainProperties(ovk_domain *Domain, ovk_domain_properties **Properties) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(*Properties == &Domain->properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(EditingProperties(Domain), "Unable to release properties; not currently being "
    "edited.");

  --Domain->properties_edit_ref_count;
  bool EndEdit = Domain->properties_edit_ref_count == 0;

  *Properties = NULL;

  if (EndEdit) {
    MPI_Barrier(Domain->properties.comm);
  }

}

void ovkGetNextAvailableGridID(const ovk_domain *Domain, int *GridID) {

  *GridID = OMNextAvailableKey(Domain->grids);

}

void ovkCreateGridLocal(ovk_domain *Domain, int GridID, const ovk_grid_params *Params) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  OVK_DEBUG_ASSERT(!EditingProperties(Domain), "Cannot create grid while editing properties.");
  OVK_DEBUG_ASSERT(!EditingGrid(Domain, OVK_ALL_GRIDS), "Cannot create grid while editing other "
    "grids.");
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot create grid "
    "while editing connectivities.");

  OVK_DEBUG_ASSERT(!OMExists(Domain->grids, GridID), "Grid %i already exists.", GridID);

  MPI_Barrier(Domain->properties.comm);

  t_domain_grid_container *GridContainer;
  CreateGridGlobal(&GridContainer, GridID, Params, Domain->properties.comm,
    Domain->properties.comm_rank, Domain->logger, Domain->error_handler);

  OMInsert(Domain->grids, GridID, GridContainer);

  ++Domain->properties.num_grids;

  if (Domain->config & OVK_DOMAIN_CONFIG_CONNECTIVITY) {
    CreateConnectivitiesForGrid(Domain, GridID);
  }

  if (Domain->config & OVK_DOMAIN_CONFIG_EXCHANGE) {
    CreateExchangesForGrid(Domain, GridID);
  }

  MPI_Barrier(Domain->properties.comm);

}

void ovkCreateGridRemote(ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");

  MPI_Barrier(Domain->properties.comm);

  OVK_DEBUG_ASSERT(!EditingProperties(Domain), "Cannot create grid while editing properties.");
  OVK_DEBUG_ASSERT(!EditingGrid(Domain, OVK_ALL_GRIDS), "Cannot create grid while editing other "
    "grids.");
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot create grid "
    "while editing connectivities.");

  OVK_DEBUG_ASSERT(!OMExists(Domain->grids, GridID), "Grid %i already exists.", GridID);

  t_domain_grid_container *GridContainer;
  CreateGridGlobal(&GridContainer, GridID, NULL, Domain->properties.comm,
    Domain->properties.comm_rank, Domain->logger, Domain->error_handler);

  OMInsert(Domain->grids, GridID, GridContainer);

  ++Domain->properties.num_grids;

  if (Domain->config & OVK_DOMAIN_CONFIG_CONNECTIVITY) {
    CreateConnectivitiesForGrid(Domain, GridID);
  }

  if (Domain->config & OVK_DOMAIN_CONFIG_EXCHANGE) {
    CreateExchangesForGrid(Domain, GridID);
  }

  MPI_Barrier(Domain->properties.comm);

}

static void CreateGridGlobal(t_domain_grid_container **GridContainer, int GridID,
  const ovk_grid_params *Params, MPI_Comm Comm, int CommRank, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  bool IsLocal = Params != NULL;

  if (OVK_DEBUG) {
    int IsLocalInt = IsLocal ? 1 : 0;
    int AtLeastOneLocal;
    MPI_Allreduce(&IsLocalInt, &AtLeastOneLocal, 1, MPI_INT, MPI_LOR, Comm);
    OVK_DEBUG_ASSERT(AtLeastOneLocal, "Grid must exist on at least one process.");
  }

  ovk_grid *Grid = NULL;
  if (IsLocal) {
    CreateGrid(&Grid, GridID, Params, Logger, ErrorHandler);
  }

  CreateGridContainer(GridContainer, Grid, Comm, CommRank);

}

void ovkDestroyGrid(ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  MPI_Barrier(Domain->properties.comm);

  OVK_DEBUG_ASSERT(!EditingProperties(Domain), "Cannot destroy grid while editing properties.");
  OVK_DEBUG_ASSERT(!EditingGrid(Domain, OVK_ALL_GRIDS), "Cannot destroy grid while editing other "
    "grids.");
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot destroy grid "
    "while editing connectivities.");

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");

  OVK_DEBUG_ASSERT(OMExists(Domain->grids, GridID), "Grid %i does not exist.", GridID);

  if (Domain->config & OVK_DOMAIN_CONFIG_EXCHANGE) {
    DestroyExchangesForGrid(Domain, GridID);
  }

  if (Domain->config & OVK_DOMAIN_CONFIG_CONNECTIVITY) {
    DestroyConnectivitiesForGrid(Domain, GridID);
  }

  t_domain_grid_container *GridContainer = OMRemove(Domain->grids, GridID);

  DestroyGridGlobal(&GridContainer);

  --Domain->properties.num_grids;

  MPI_Barrier(Domain->properties.comm);

}

static void DestroyGridGlobal(t_domain_grid_container **GridContainer_) {

  t_domain_grid_container *GridContainer = *GridContainer_;

  bool IsLocal = GridContainer->has_local_data;

  if (IsLocal) {
    DestroyGrid(&GridContainer->grid);
  }

  DestroyGridContainer(GridContainer_);

}

static t_domain_grid_container *FindGrid(ovk_domain *Domain, int GridID) {

  t_ordered_map_entry *Entry = OMFind(Domain->grids, GridID);

  if (Entry == OMEnd(Domain->grids)) return NULL;

  return OMData(Entry);

}

static const t_domain_grid_container *FindGridC(const ovk_domain *Domain, int GridID) {

  const t_ordered_map_entry *Entry = OMFindC(Domain->grids, GridID);

  if (Entry == OMEndC(Domain->grids)) return NULL;

  return OMDataC(Entry);

}

bool ovkGridExists(const ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");

  return FindGridC(Domain, GridID) != NULL;

}

void ovkGetGridInfo(const ovk_domain *Domain, int GridID, const ovk_grid_info **GridInfo) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(GridInfo, "Invalid grid info pointer.");

  const t_domain_grid_container *Container = FindGridC(Domain, GridID);

  *GridInfo = Container->info;

}

bool ovkRankHasGrid(const ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");

  const t_domain_grid_container *Container = FindGridC(Domain, GridID);
  OVK_DEBUG_ASSERT(Container, "Grid %i does not exist.", GridID);

  return Container->has_local_data;

}

void ovkGetGrid(const ovk_domain *Domain, int GridID, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  const t_domain_grid_container *Container = FindGridC(Domain, GridID);
  OVK_DEBUG_ASSERT(Container, "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(Container->has_local_data, "Grid %s does not have local data on rank @rank@.",
    Container->info->name);

  *Grid = Container->grid;

}

void ovkEditGridLocal(ovk_domain *Domain, int GridID, ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  EditGridGlobal(Domain, GridID, Grid);

}

void ovkEditGridRemote(ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, GridID), "Grid %i does not exist.", GridID);

  EditGridGlobal(Domain, GridID, NULL);

}

static void EditGridGlobal(ovk_domain *Domain, int GridID, ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(!EditingProperties(Domain), "Cannot edit grid while editing properties.");
  OVK_DEBUG_ASSERT(!EditingConnectivity(Domain, OVK_ALL_GRIDS, OVK_ALL_GRIDS), "Cannot edit "
    "grid while editing connectivities.");

  bool IsLocal = Grid != NULL;

  t_domain_grid_container *Container = FindGrid(Domain, GridID);

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Container->has_local_data, "Grid %s does not have local data on rank @rank@.",
      Container->info->name);
  }

//   bool StartEditAll = Domain->grids_edit_ref_count == 0;
  ++Domain->grids_edit_ref_count;
  bool StartEdit = Container->edit_ref_count == 0;
  ++Container->edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Domain->properties.comm);
  }

  if (IsLocal) {
    *Grid = Container->grid;
  }

}

void ovkReleaseGridLocal(ovk_domain *Domain, int GridID, ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, GridID), "Grid %i does not exist.", GridID);
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(*Grid, "Invalid grid pointer.");

  ReleaseGridGlobal(Domain, GridID, Grid);

}

void ovkReleaseGridRemote(ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, GridID), "Grid %i does not exist.", GridID);

  ReleaseGridGlobal(Domain, GridID, NULL);

}

static void ReleaseGridGlobal(ovk_domain *Domain, int GridID, ovk_grid **Grid) {

  bool IsLocal = Grid != NULL;

  t_domain_grid_container *Container = FindGrid(Domain, GridID);

  OVK_DEBUG_ASSERT(EditingGrid(Domain, GridID), "Unable to release grid %s; not currently being "
    "edited.", Container->info->name);

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Container->has_local_data, "Grid %s does not have local data on rank @rank@.",
      Container->info->name);
    OVK_DEBUG_ASSERT(*Grid == Container->grid, "Invalid grid pointer.");
  }

  --Container->edit_ref_count;
  bool EndEdit = Container->edit_ref_count == 0;
  --Domain->grids_edit_ref_count;
//   bool EndEditAll = Domain->grids_edit_ref_count == 0;

  if (IsLocal) {
    *Grid = NULL;
  }

  if (EndEdit) {
    MPI_Barrier(Domain->properties.comm);
  }

}

static void EnableConnectivityComponent(ovk_domain *Domain) {

  OMCreate(&Domain->connectivities);

  const t_ordered_map_entry *DonorGridEntry = OMBeginC(Domain->grids);

  while (DonorGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *DonorGridContainer = OMDataC(DonorGridEntry);
    int DonorGridID = DonorGridContainer->info->id;

    t_ordered_map *ConnectivityRow;
    OMCreate(&ConnectivityRow);
    OMInsert(Domain->connectivities, DonorGridID, ConnectivityRow);

    const t_ordered_map_entry *ReceiverGridEntry = OMBeginC(Domain->grids);

    while (ReceiverGridEntry != OMEndC(Domain->grids)) {

      const t_domain_grid_container *ReceiverGridContainer = OMDataC(ReceiverGridEntry);
      int ReceiverGridID = ReceiverGridContainer->info->id;

      if (ReceiverGridID != DonorGridID) {
        t_domain_connectivity_container *ConnectivityContainer;
        CreateConnectivityGlobal(&ConnectivityContainer, Domain->properties.num_dims,
          DonorGridContainer, ReceiverGridContainer, Domain->properties.comm,
          Domain->properties.comm_rank, Domain->logger, Domain->error_handler);
        OMInsert(ConnectivityRow, ReceiverGridID, ConnectivityContainer);
      }

      ReceiverGridEntry = OMNextC(ReceiverGridEntry);

    }

    DonorGridEntry = OMNextC(DonorGridEntry);

  }

}

static void DisableConnectivityComponent(ovk_domain *Domain) {

  t_ordered_map_entry *RowEntry = OMBegin(Domain->connectivities);

  while (RowEntry != OMEnd(Domain->connectivities)) {

    t_ordered_map *ConnectivityRow = OMData(RowEntry);
    t_ordered_map_entry *Entry = OMBegin(ConnectivityRow);

    while (Entry != OMEnd(ConnectivityRow)) {
      t_domain_connectivity_container *ConnectivityContainer = OMData(Entry);
      DestroyConnectivityGlobal(&ConnectivityContainer);
      Entry = OMNext(Entry);
    }

    OMDestroy(&ConnectivityRow);

    RowEntry = OMNext(RowEntry);

  }

  OMDestroy(&Domain->connectivities);

}

static void CreateConnectivitiesForGrid(ovk_domain *Domain, int GridID) {

  const t_domain_grid_container *GridContainer = OMDataC(OMFindC(Domain->grids, GridID));

  t_ordered_map *ConnectivityRow;
  OMCreate(&ConnectivityRow);
  OMInsert(Domain->connectivities, GridID, ConnectivityRow);

  const t_ordered_map_entry *OtherGridEntry = OMBeginC(Domain->grids);

  while (OtherGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *OtherGridContainer = OMDataC(OtherGridEntry);
    int OtherGridID = OtherGridContainer->info->id;

    if (OtherGridID != GridID) {
      t_domain_connectivity_container *ConnectivityContainer;
      CreateConnectivityGlobal(&ConnectivityContainer, Domain->properties.num_dims, GridContainer,
        OtherGridContainer, Domain->properties.comm, Domain->properties.comm_rank, Domain->logger,
        Domain->error_handler);
      OMInsert(ConnectivityRow, OtherGridID, ConnectivityContainer);
    }

    OtherGridEntry = OMNextC(OtherGridEntry);

  }

  OtherGridEntry = OMBeginC(Domain->grids);

  while (OtherGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *OtherGridContainer = OMDataC(OtherGridEntry);
    int OtherGridID = OtherGridContainer->info->id;

    if (OtherGridID != GridID) {
      ConnectivityRow = OMData(OMFind(Domain->connectivities, OtherGridID));
      t_domain_connectivity_container *ConnectivityContainer;
      CreateConnectivityGlobal(&ConnectivityContainer, Domain->properties.num_dims,
        OtherGridContainer, GridContainer, Domain->properties.comm, Domain->properties.comm_rank,
        Domain->logger, Domain->error_handler);
      OMInsert(ConnectivityRow, GridID, ConnectivityContainer);
    }

    OtherGridEntry = OMNextC(OtherGridEntry);

  }

}

static void DestroyConnectivitiesForGrid(ovk_domain *Domain, int GridID) {

  t_ordered_map *ConnectivityRow = OMRemove(Domain->connectivities, GridID);

  const t_ordered_map_entry *OtherGridEntry = OMBeginC(Domain->grids);

  while (OtherGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *OtherGridContainer = OMDataC(OtherGridEntry);
    int OtherGridID = OtherGridContainer->info->id;

    if (OtherGridID != GridID) {
      t_domain_connectivity_container *ConnectivityContainer = OMData(OMFind(ConnectivityRow,
        OtherGridID));
      DestroyConnectivityGlobal(&ConnectivityContainer);
    }

    OtherGridEntry = OMNextC(OtherGridEntry);

  }

  OMDestroy(&ConnectivityRow);

  OtherGridEntry = OMBeginC(Domain->grids);

  while (OtherGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *OtherGridContainer = OMDataC(OtherGridEntry);
    int OtherGridID = OtherGridContainer->info->id;

    if (OtherGridID != GridID) {
      ConnectivityRow = OMData(OMFind(Domain->connectivities, OtherGridID));
      t_domain_connectivity_container *ConnectivityContainer = OMRemove(ConnectivityRow, GridID);
      DestroyConnectivityGlobal(&ConnectivityContainer);
    }

    OtherGridEntry = OMNextC(OtherGridEntry);

  }

}

static void CreateConnectivityGlobal(t_domain_connectivity_container **ConnectivityContainer,
  int NumDims, const t_domain_grid_container *DonorGridContainer,
  const t_domain_grid_container *ReceiverGridContainer, MPI_Comm Comm, int CommRank,
  t_logger *Logger, t_error_handler *ErrorHandler) {

  bool IsLocal = DonorGridContainer->has_local_data || ReceiverGridContainer->has_local_data;

  MPI_Comm ConnectivityComm;
  MPI_Comm_split(Comm, IsLocal, CommRank, &ConnectivityComm);

  ovk_connectivity *Connectivity = NULL;
  if (IsLocal) {
    const ovk_grid *DonorGrid = NULL;
    if (DonorGridContainer->has_local_data) DonorGrid = DonorGridContainer->grid;
    const ovk_grid *ReceiverGrid = NULL;
    if (ReceiverGridContainer->has_local_data) ReceiverGrid = ReceiverGridContainer->grid;
    CreateConnectivity(&Connectivity, NumDims, ConnectivityComm, DonorGrid, ReceiverGrid, Logger,
      ErrorHandler);
  }

  CreateConnectivityContainer(ConnectivityContainer, Connectivity, Comm, CommRank);

  MPI_Comm_free(&ConnectivityComm);

}

static void DestroyConnectivityGlobal(t_domain_connectivity_container **ConnectivityContainer_) {

  t_domain_connectivity_container *ConnectivityContainer = *ConnectivityContainer_;

  bool IsLocal = ConnectivityContainer->has_local_data;

  if (IsLocal) {
    DestroyConnectivity(&ConnectivityContainer->connectivity);
  }

  DestroyConnectivityContainer(ConnectivityContainer_);

}

static t_domain_connectivity_container *FindConnectivity(ovk_domain *Domain, int DonorGridID,
  int ReceiverGridID) {

  t_ordered_map_entry *RowEntry = OMFind(Domain->connectivities, DonorGridID);

  if (RowEntry == OMEnd(Domain->connectivities)) return NULL;

  t_ordered_map *ConnectivityRow = OMData(RowEntry);
  t_ordered_map_entry *Entry = OMFind(ConnectivityRow, ReceiverGridID);

  if (Entry == OMEnd(ConnectivityRow)) return NULL;

  return OMData(Entry);

}

static const t_domain_connectivity_container *FindConnectivityC(const ovk_domain *Domain,
  int DonorGridID, int ReceiverGridID) {

  const t_ordered_map_entry *RowEntry = OMFindC(Domain->connectivities, DonorGridID);

  if (RowEntry == OMEndC(Domain->connectivities)) return NULL;

  const t_ordered_map *ConnectivityRow = OMDataC(RowEntry);
  const t_ordered_map_entry *Entry = OMFindC(ConnectivityRow, ReceiverGridID);

  if (Entry == OMEndC(ConnectivityRow)) return NULL;

  return OMDataC(Entry);

}

bool ovkConnectivityExists(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);

  return FindConnectivityC(Domain, DonorGridID, ReceiverGridID) != NULL;

}

bool ovkRankHasConnectivity(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);

  const t_domain_connectivity_container *Container = FindConnectivityC(Domain, DonorGridID,
    ReceiverGridID);
  OVK_DEBUG_ASSERT(Container, "Connectivity (%i,%i) does not exist.", DonorGridID, ReceiverGridID);

  return Container->has_local_data;

}

void ovkGetConnectivity(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  const ovk_connectivity **Connectivity) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);

  const t_domain_connectivity_container *Container = FindConnectivityC(Domain, DonorGridID,
    ReceiverGridID);
  OVK_DEBUG_ASSERT(Container, "Connectivity (%i,%i) does not exist.", DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container->has_local_data, "Connectivity %s does not have local data on "
    "rank @rank@.", Container->info->name);

  *Connectivity = Container->connectivity;

}

void ovkEditConnectivityLocal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity "
    "(%i,%i) does not exist.", DonorGridID, ReceiverGridID);

  EditConnectivityGlobal(Domain, DonorGridID, ReceiverGridID, Connectivity);

}

void ovkEditConnectivityRemote(ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity "
    "(%i,%i) does not exist.", DonorGridID, ReceiverGridID);

  EditConnectivityGlobal(Domain, DonorGridID, ReceiverGridID, NULL);

}

static void EditConnectivityGlobal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity) {

  OVK_DEBUG_ASSERT(!EditingProperties(Domain), "Cannot edit connectivity while editing properties.");
  OVK_DEBUG_ASSERT(!EditingGrid(Domain, OVK_ALL_GRIDS), "Cannot edit connectivity while editing "
    "grids.");

  bool IsLocal = Connectivity != NULL;

  t_domain_connectivity_container *Container = FindConnectivity(Domain, DonorGridID,
    ReceiverGridID);

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Container->has_local_data, "Connectivity %s does not have local data on "
      "rank @rank@.", Container->info->name);
  }

//   bool StartEditAll = Domain->connectivities_edit_ref_count == 0;
  ++Domain->connectivities_edit_ref_count;
  bool StartEdit = Container->edit_ref_count == 0;
  ++Container->edit_ref_count;

  if (StartEdit) {
    MPI_Barrier(Domain->properties.comm);
  }

  if (IsLocal) {
    *Connectivity = Container->connectivity;
  }

}

void ovkReleaseConnectivityLocal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity "
    "(%i,%i) does not exist.", DonorGridID, ReceiverGridID);

  ReleaseConnectivityGlobal(Domain, DonorGridID, ReceiverGridID, Connectivity);

}

void ovkReleaseConnectivityRemote(ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID), "Connectivity "
    "(%i,%i) does not exist.", DonorGridID, ReceiverGridID);

  ReleaseConnectivityGlobal(Domain, DonorGridID, ReceiverGridID, NULL);

}

static void ReleaseConnectivityGlobal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity) {

  bool IsLocal = Connectivity != NULL;

  t_domain_connectivity_container *Container = FindConnectivity(Domain, DonorGridID,
    ReceiverGridID);

  OVK_DEBUG_ASSERT(EditingConnectivity(Domain, DonorGridID, ReceiverGridID), "Unable to release "
    "connectivity %s; not currently being edited.", Container->info->name);

  if (OVK_DEBUG && IsLocal) {
    OVK_DEBUG_ASSERT(Container->has_local_data, "Connectivity %s does not have local data on "
      "rank @rank@.", Container->info->name);
    OVK_DEBUG_ASSERT(*Connectivity == Container->connectivity, "Invalid connectivity pointer.");
  }

  --Container->edit_ref_count;
  bool EndEdit = Container->edit_ref_count == 0;
  --Domain->connectivities_edit_ref_count;
//   bool EndEditAll = Domain->connectivities_edit_ref_count == 0;

  if (IsLocal) {
    *Connectivity = NULL;
  }

  if (EndEdit) {
    MPI_Barrier(Domain->properties.comm);
  }

}

static void EnableExchangeComponent(ovk_domain *Domain) {

  OMCreate(&Domain->exchanges);

  const t_ordered_map_entry *DonorGridEntry = OMBeginC(Domain->grids);

  while (DonorGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *DonorGridContainer = OMDataC(DonorGridEntry);
    int DonorGridID = DonorGridContainer->info->id;

    t_ordered_map *ExchangeRow;
    OMCreate(&ExchangeRow);
    OMInsert(Domain->exchanges, DonorGridID, ExchangeRow);

    const t_ordered_map_entry *ReceiverGridEntry = OMBeginC(Domain->grids);

    while (ReceiverGridEntry != OMEndC(Domain->grids)) {

      const t_domain_grid_container *ReceiverGridContainer = OMDataC(ReceiverGridEntry);
      int ReceiverGridID = ReceiverGridContainer->info->id;

      if (ReceiverGridID != DonorGridID) {
        const t_domain_connectivity_container *ConnectivityContainer = FindConnectivityC(Domain,
          DonorGridID, ReceiverGridID);
        t_domain_exchange_container *ExchangeContainer;
        CreateExchangeGlobal(&ExchangeContainer, ConnectivityContainer, Domain->properties.comm,
          Domain->properties.comm_rank, Domain->logger, Domain->error_handler);
        OMInsert(ExchangeRow, ReceiverGridID, ExchangeContainer);
      }

      ReceiverGridEntry = OMNextC(ReceiverGridEntry);

    }

    DonorGridEntry = OMNextC(DonorGridEntry);

  }

}

static void DisableExchangeComponent(ovk_domain *Domain) {

  t_ordered_map_entry *RowEntry = OMBegin(Domain->exchanges);

  while (RowEntry != OMEnd(Domain->exchanges)) {

    t_ordered_map *ExchangeRow = OMData(RowEntry);
    t_ordered_map_entry *Entry = OMBegin(ExchangeRow);

    while (Entry != OMEnd(ExchangeRow)) {
      t_domain_exchange_container *ExchangeContainer = OMData(Entry);
      DestroyExchangeGlobal(&ExchangeContainer);
      Entry = OMNext(Entry);
    }

    OMDestroy(&ExchangeRow);

    RowEntry = OMNext(RowEntry);

  }

  OMDestroy(&Domain->exchanges);

}

static void CreateExchangesForGrid(ovk_domain *Domain, int GridID) {

  t_ordered_map *ExchangeRow;
  OMCreate(&ExchangeRow);
  OMInsert(Domain->exchanges, GridID, ExchangeRow);

  const t_ordered_map_entry *OtherGridEntry = OMBeginC(Domain->grids);

  while (OtherGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *OtherGridContainer = OMDataC(OtherGridEntry);
    int OtherGridID = OtherGridContainer->info->id;

    if (OtherGridID != GridID) {
      const t_domain_connectivity_container *ConnectivityContainer = FindConnectivityC(Domain,
        GridID, OtherGridID);
      t_domain_exchange_container *ExchangeContainer;
      CreateExchangeGlobal(&ExchangeContainer, ConnectivityContainer, Domain->properties.comm,
        Domain->properties.comm_rank, Domain->logger, Domain->error_handler);
      OMInsert(ExchangeRow, OtherGridID, ExchangeContainer);
    }

    OtherGridEntry = OMNextC(OtherGridEntry);

  }

  OtherGridEntry = OMBeginC(Domain->grids);

  while (OtherGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *OtherGridContainer = OMDataC(OtherGridEntry);
    int OtherGridID = OtherGridContainer->info->id;

    if (OtherGridID != GridID) {
      const t_domain_connectivity_container *ConnectivityContainer = FindConnectivityC(Domain,
        OtherGridID, GridID);
      ExchangeRow = OMData(OMFind(Domain->exchanges, OtherGridID));
      t_domain_exchange_container *ExchangeContainer;
      CreateExchangeGlobal(&ExchangeContainer, ConnectivityContainer, Domain->properties.comm,
        Domain->properties.comm_rank, Domain->logger, Domain->error_handler);
      OMInsert(ExchangeRow, GridID, ExchangeContainer);
    }

    OtherGridEntry = OMNextC(OtherGridEntry);

  }

}

static void DestroyExchangesForGrid(ovk_domain *Domain, int GridID) {

  t_ordered_map *ExchangeRow = OMRemove(Domain->exchanges, GridID);

  const t_ordered_map_entry *OtherGridEntry = OMBeginC(Domain->grids);

  while (OtherGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *OtherGridContainer = OMDataC(OtherGridEntry);
    int OtherGridID = OtherGridContainer->info->id;

    if (OtherGridID != GridID) {
      t_domain_exchange_container *ExchangeContainer = OMData(OMFind(ExchangeRow, OtherGridID));
      DestroyExchangeGlobal(&ExchangeContainer);
    }

    OtherGridEntry = OMNextC(OtherGridEntry);

  }

  OMDestroy(&ExchangeRow);

  OtherGridEntry = OMBeginC(Domain->grids);

  while (OtherGridEntry != OMEndC(Domain->grids)) {

    const t_domain_grid_container *OtherGridContainer = OMDataC(OtherGridEntry);
    int OtherGridID = OtherGridContainer->info->id;

    if (OtherGridID != GridID) {
      ExchangeRow = OMData(OMFind(Domain->exchanges, OtherGridID));
      t_domain_exchange_container *ExchangeContainer = OMRemove(ExchangeRow, GridID);
      DestroyExchangeGlobal(&ExchangeContainer);
    }

    OtherGridEntry = OMNextC(OtherGridEntry);

  }

}

static void CreateExchangeGlobal(t_domain_exchange_container **ExchangeContainer,
  const t_domain_connectivity_container *ConnectivityContainer, MPI_Comm Comm, int CommRank,
  t_logger *Logger, t_error_handler *ErrorHandler) {

  bool IsLocal = ConnectivityContainer->has_local_data;

  ovk_exchange *Exchange = NULL;
  if (IsLocal) {
    ovk_connectivity *Connectivity = ConnectivityContainer->connectivity;
    CreateExchange(&Exchange, Connectivity, Logger, ErrorHandler);
  }

  CreateExchangeContainer(ExchangeContainer, Exchange, Comm, CommRank);

}

static void DestroyExchangeGlobal(t_domain_exchange_container **ExchangeContainer_) {

  t_domain_exchange_container *ExchangeContainer = *ExchangeContainer_;

  bool IsLocal = ExchangeContainer->has_local_data;

  if (IsLocal) {
    DestroyExchange(&ExchangeContainer->exchange);
  }

  DestroyExchangeContainer(ExchangeContainer_);

}

// static t_domain_exchange_container *FindExchange(ovk_domain *Domain, int DonorGridID,
//   int ReceiverGridID) {

//   t_ordered_map_entry *RowEntry = OMFind(Domain->exchanges, DonorGridID);

//   if (RowEntry == OMEnd(Domain->exchanges)) return NULL;

//   t_ordered_map *ExchangeRow = OMData(RowEntry);
//   t_ordered_map_entry *Entry = OMFind(ExchangeRow, ReceiverGridID);

//   if (Entry == OMEnd(ExchangeRow)) return NULL;

//   return OMData(Entry);

// }

static const t_domain_exchange_container *FindExchangeC(const ovk_domain *Domain, int DonorGridID,
  int ReceiverGridID) {

  const t_ordered_map_entry *RowEntry = OMFindC(Domain->exchanges, DonorGridID);

  if (RowEntry == OMEndC(Domain->exchanges)) return NULL;

  const t_ordered_map *ExchangeRow = OMDataC(RowEntry);
  const t_ordered_map_entry *Entry = OMFindC(ExchangeRow, ReceiverGridID);

  if (Entry == OMEndC(ExchangeRow)) return NULL;

  return OMDataC(Entry);

}

bool ovkExchangeExists(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);

  return FindExchangeC(Domain, DonorGridID, ReceiverGridID) != NULL;

}

bool ovkRankHasExchange(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);

  const t_domain_exchange_container *Container = FindExchangeC(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container, "Exchange (%i,%i) does not exist.", DonorGridID, ReceiverGridID);

  return Container->has_local_data;

}

void ovkGetExchange(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  const ovk_exchange **Exchange) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(Exchange, "Invalid exchange pointer.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);

  const t_domain_exchange_container *Container = FindExchangeC(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container, "Exchange (%i,%i) does not exist.", DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container->has_local_data, "Exchange %s does not have local data on rank @rank@.",
    Container->info->name);

  *Exchange = Container->exchange;

}

static bool EditingProperties(const ovk_domain *Domain) {

  return Domain->properties_edit_ref_count > 0;

}

static bool EditingGrid(const ovk_domain *Domain, int GridID) {

  if (GridID == OVK_ALL_GRIDS) {
    return Domain->grids_edit_ref_count > 0;
  } else {
    const t_domain_grid_container *Container = FindGridC(Domain, GridID);
    return Container->edit_ref_count > 0;
  }

}

static bool EditingConnectivity(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  if (DonorGridID == OVK_ALL_GRIDS && ReceiverGridID == OVK_ALL_GRIDS) {
    return Domain->connectivities_edit_ref_count > 0;
  } else {
    const t_domain_connectivity_container *Container = FindConnectivityC(Domain, DonorGridID,
      ReceiverGridID);
    return Container->edit_ref_count > 0;
  }

}

void ovkGetLocalDonorCount(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  size_t *NumDonors) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(NumDonors, "Invalid num donors pointer.");

  const ovk_connectivity *Connectivity;
  ovkGetConnectivity(Domain, DonorGridID, ReceiverGridID, &Connectivity);

  const ovk_connectivity_d *Donors;
  ovkGetConnectivityDonorSide(Connectivity, &Donors);

  const ovk_connectivity_d_properties *DonorsProperties;
  ovkGetConnectivityDonorSideProperties(Donors, &DonorsProperties);

  ovkGetConnectivityDonorSidePropertyCount(DonorsProperties, NumDonors);

}

void ovkGetLocalReceiverCount(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  size_t *NumReceivers) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(NumReceivers, "Invalid num receivers pointer.");

  const ovk_connectivity *Connectivity;
  ovkGetConnectivity(Domain, DonorGridID, ReceiverGridID, &Connectivity);

  const ovk_connectivity_r *Receivers;
  ovkGetConnectivityReceiverSide(Connectivity, &Receivers);

  const ovk_connectivity_r_properties *ReceiversProperties;
  ovkGetConnectivityReceiverSideProperties(Receivers, &ReceiversProperties);

  ovkGetConnectivityReceiverSidePropertyCount(ReceiversProperties, NumReceivers);

}

void ovkAssemble(ovk_domain *Domain, const ovk_assembly_options *Options) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  bool IsDomainRoot = Domain->properties.comm_rank == 0;

  LogStatus(Domain->logger, IsDomainRoot, 0, "Beginning overset assembly on domain %s.",
    Domain->properties.name);

  bool HasOverlap = Domain->config & OVK_DOMAIN_CONFIG_OVERLAP;
  bool HasConnectivity = Domain->config & OVK_DOMAIN_CONFIG_CONNECTIVITY;
  bool HasExchange = Domain->config & OVK_DOMAIN_CONFIG_EXCHANGE;

  if (HasOverlap) {
//     AssembleOverlap(Domain);
  }

  if (HasOverlap && HasConnectivity) {
//     AssembleConnectivity(Domain);
  }

  if (HasExchange) {
    AssembleExchange(Domain);
  }

  if (HasConnectivity) {
    ResetAllConnectivityEdits(Domain);
  }

  LogStatus(Domain->logger, IsDomainRoot, 0, "Finished overset assembly on domain %s.",
    Domain->properties.name);

}

static void AssembleExchange(ovk_domain *Domain) {

  t_ordered_map_entry *RowEntry = OMBegin(Domain->exchanges);
  while (RowEntry != OMEnd(Domain->exchanges)) {
    t_ordered_map *ExchangeRow = OMData(RowEntry);
    t_ordered_map_entry *Entry = OMBegin(ExchangeRow);
    while (Entry != OMEnd(ExchangeRow)) {
      t_domain_exchange_container *Container = OMData(Entry);
      if (Container->has_local_data) {
        UpdateExchange(Container->exchange);
      }
      Entry = OMNext(Entry);
    }
    RowEntry = OMNext(RowEntry);
  }

}

static void ResetAllConnectivityEdits(ovk_domain *Domain) {

  t_ordered_map_entry *RowEntry = OMBegin(Domain->connectivities);
  while (RowEntry != OMEnd(Domain->connectivities)) {
    t_ordered_map *ConnectivityRow = OMData(RowEntry);
    t_ordered_map_entry *Entry = OMBegin(ConnectivityRow);
    while (Entry != OMEnd(ConnectivityRow)) {
      t_domain_connectivity_container *Container = OMData(Entry);
      if (Container->has_local_data) {
        ResetConnectivityEdits(Container->connectivity);
      }
      Entry = OMNext(Entry);
    }
    RowEntry = OMNext(RowEntry);
  }

}

void ovkCollect(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_data_type DataType, int Count, ovk_collect_op CollectOp, const ovk_range *GridDataRange,
  ovk_array_layout GridDataLayout, const void **GridData, void **DonorData) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Domain->config & OVK_DOMAIN_CONFIG_EXCHANGE, "Domain %s is not configured for "
    "exchange.", Domain->properties.name);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ValidDataType(DataType), "Invalid data type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(ValidCollectOp(CollectOp), "Invalid collect operation.");
  // Note: Checking GridDataRange, GridData, and DonorData inside ExchangeCollect
  OVK_DEBUG_ASSERT(ValidArrayLayout(GridDataLayout), "Invalid grid data layout.");

  const t_domain_exchange_container *Container = FindExchangeC(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container, "Exchange (%i,%i) does not exist.", DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container->has_local_data, "Exchange %s does not have local data on rank @rank@.",
    Container->info->name);

  const ovk_exchange *Exchange = Container->exchange;

  OVK_DEBUG_ASSERT(ovkRankHasExchangeDonorSide(Exchange), "Exchange %s does not have "
    "donor-side data on rank @rank@.", Container->info->name);

  ExchangeCollect(Exchange, DataType, Count, CollectOp, GridDataRange, GridDataLayout, GridData,
    DonorData);

}

void ovkSend(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID, ovk_data_type DataType,
  int Count, const void **DonorData, int Tag, ovk_request **Request) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Domain->config & OVK_DOMAIN_CONFIG_EXCHANGE, "Domain %s is not configured for "
    "exchange.", Domain->properties.name);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ValidDataType(DataType), "Invalid data type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  // Note: Checking DonorData inside ExchangeSend
  OVK_DEBUG_ASSERT(Tag >= 0, "Invalid tag.");
  OVK_DEBUG_ASSERT(Request, "Invalid request pointer.");

  const t_domain_exchange_container *Container = FindExchangeC(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container, "Exchange (%i,%i) does not exist.", DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container->has_local_data, "Exchange %s does not have local data on rank @rank@.",
    Container->info->name);

  const ovk_exchange *Exchange = Container->exchange;

  OVK_DEBUG_ASSERT(ovkRankHasExchangeDonorSide(Exchange), "Exchange %s does not have "
    "donor-side data on rank @rank@.", Container->info->name);

  ExchangeSend(Exchange, DataType, Count, DonorData, Tag, Request);

}

void ovkReceive(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_data_type DataType, int Count, void **ReceiverData, int Tag, ovk_request **Request) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Domain->config & OVK_DOMAIN_CONFIG_EXCHANGE, "Domain %s is not configured for "
    "exchange.", Domain->properties.name);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ValidDataType(DataType), "Invalid data type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  // Note: Checking ReceiverData inside ExchangeReceive
  OVK_DEBUG_ASSERT(Tag >= 0, "Invalid tag.");
  OVK_DEBUG_ASSERT(Request, "Invalid request pointer.");

  const t_domain_exchange_container *Container = FindExchangeC(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container, "Exchange (%i,%i) does not exist.", DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container->has_local_data, "Exchange %s does not have local data on rank @rank@.",
    Container->info->name);

  const ovk_exchange *Exchange = Container->exchange;

  OVK_DEBUG_ASSERT(ovkRankHasExchangeReceiverSide(Exchange), "Exchange %s does not have "
    "receiver-side data on rank @rank@.", Container->info->name);

  ExchangeReceive(Exchange, DataType, Count, ReceiverData, Tag, Request);

}

void ovkWaitAll(int NumRequests, ovk_request **Requests) {

  OVK_DEBUG_ASSERT(NumRequests >= 0, "Invalid request count.");
  OVK_DEBUG_ASSERT(NumRequests == 0 || Requests, "Invalid requests pointer.");
  // Note: Not checking Requests[i] here on purpose -- allowed to be NULL

  ExchangeWaitAll(NumRequests, Requests);

}

void ovkWaitAny(int NumRequests, ovk_request **Requests, int *Index) {

  OVK_DEBUG_ASSERT(NumRequests >= 0, "Invalid request count.");
  OVK_DEBUG_ASSERT(NumRequests == 0 || Requests, "Invalid requests pointer.");
  // Note: Not checking Requests[i] here on purpose -- allowed to be NULL
  OVK_DEBUG_ASSERT(Index, "Invalid index pointer.");

  ExchangeWaitAny(NumRequests, Requests, Index);

}

void ovkDisperse(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_data_type DataType, int Count, ovk_disperse_op DisperseOp, const void **ReceiverData,
  const ovk_range *GridDataRange, ovk_array_layout GridDataLayout, void **GridData) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Domain->config & OVK_DOMAIN_CONFIG_EXCHANGE, "Domain %s is not configured for "
    "exchange.", Domain->properties.name);
  OVK_DEBUG_ASSERT(DonorGridID >= 0, "Invalid donor grid ID.");
  OVK_DEBUG_ASSERT(ReceiverGridID >= 0, "Invalid receiver grid ID.");
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, DonorGridID), "Grid %i does not exist.", DonorGridID);
  OVK_DEBUG_ASSERT(ovkGridExists(Domain, ReceiverGridID), "Grid %i does not exist.", ReceiverGridID);
  OVK_DEBUG_ASSERT(ovkExchangeExists(Domain, DonorGridID, ReceiverGridID), "Exchange (%i,%i) "
    "does not exist.", DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(ValidDataType(DataType), "Invalid data type.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count.");
  OVK_DEBUG_ASSERT(ValidDisperseOp(DisperseOp), "Invalid disperse operation.");
  // Note: Checking GridDataRange, GridData, and DonorData inside ExchangeCollect
  OVK_DEBUG_ASSERT(ValidArrayLayout(GridDataLayout), "Invalid grid data layout.");

  const t_domain_exchange_container *Container = FindExchangeC(Domain, DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container, "Exchange (%i,%i) does not exist.", DonorGridID, ReceiverGridID);
  OVK_DEBUG_ASSERT(Container->has_local_data, "Exchange %s does not have local data on rank @rank@.",
    Container->info->name);

  const ovk_exchange *Exchange = Container->exchange;

  OVK_DEBUG_ASSERT(ovkRankHasExchangeReceiverSide(Exchange), "Exchange %s does not have "
    "receiver-side data on rank @rank@.", Container->info->name);

  ExchangeDisperse(Exchange, DataType, Count, DisperseOp, ReceiverData, GridDataRange,
    GridDataLayout, GridData);

}

static void CreateGridContainer(t_domain_grid_container **Container_, ovk_grid *Grid, MPI_Comm Comm,
  int CommRank) {

  bool IsLocal = Grid != NULL;

  *Container_ = malloc(sizeof(t_domain_grid_container));
  t_domain_grid_container *Container = *Container_;

  CreateGridInfo(&Container->info, Grid, Comm, CommRank);
  Container->edit_ref_count = 0;

  Container->has_local_data = IsLocal;

  if (IsLocal) {
    Container->grid = Grid;
  }

}

static void DestroyGridContainer(t_domain_grid_container **Container_) {

  t_domain_grid_container *Container = *Container_;

  DestroyGridInfo(&Container->info);

  free_null(Container_);

}

static void CreateConnectivityContainer(t_domain_connectivity_container **Container_,
  ovk_connectivity *Connectivity, MPI_Comm Comm, int CommRank) {

  bool IsLocal = Connectivity != NULL;

  *Container_ = malloc(sizeof(t_domain_connectivity_container));
  t_domain_connectivity_container *Container = *Container_;

  CreateConnectivityInfo(&Container->info, Connectivity, Comm, CommRank);
  Container->edit_ref_count = 0;

  Container->has_local_data = IsLocal;

  if (IsLocal) {
    Container->connectivity = Connectivity;
  }

}

static void DestroyConnectivityContainer(t_domain_connectivity_container **Container_) {

  t_domain_connectivity_container *Container = *Container_;

  DestroyConnectivityInfo(&Container->info);

  free_null(Container_);

}

static void CreateExchangeContainer(t_domain_exchange_container **Container_,
  ovk_exchange *Exchange, MPI_Comm Comm, int CommRank) {

  bool IsLocal = Exchange != NULL;

  *Container_ = malloc(sizeof(t_domain_exchange_container));
  t_domain_exchange_container *Container = *Container_;

  CreateExchangeInfo(&Container->info, Exchange, Comm, CommRank);

  Container->has_local_data = IsLocal;

  if (IsLocal) {
    Container->exchange = Exchange;
  }

}

static void DestroyExchangeContainer(t_domain_exchange_container **Container_) {

  t_domain_exchange_container *Container = *Container_;

  DestroyExchangeInfo(&Container->info);

  free_null(Container_);

}

void ovkGetDomainParamName(const ovk_domain_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  strcpy(Name, Params->name);

}

void ovkCreateDomainParams(ovk_domain_params **Params_, int NumDims) {

  OVK_DEBUG_ASSERT(Params_, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(NumDims == 2 || NumDims == 3, "Invalid dimension.");

  *Params_ = malloc(sizeof(ovk_domain_params));
  ovk_domain_params *Params = *Params_;

  Params->num_dims = NumDims;
  memset(Params->name, 0, OVK_NAME_LENGTH);
  Params->comm = MPI_COMM_NULL;

}

void ovkDestroyDomainParams(ovk_domain_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  free_null(Params);

}

void ovkSetDomainParamName(ovk_domain_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  strncpy(Params->name, Name, OVK_NAME_LENGTH);

}

void ovkGetDomainParamDimension(const ovk_domain_params *Params, int *NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Params->num_dims;

}

void ovkGetDomainParamComm(const ovk_domain_params *Params, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Params->comm;

}

void ovkSetDomainParamComm(ovk_domain_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Params->comm = Comm;

}

static void DefaultProperties(ovk_domain_properties *Properties) {

  memset(Properties->name, 0, OVK_NAME_LENGTH);

  Properties->num_dims = 2;
  Properties->comm = MPI_COMM_NULL;
  Properties->comm_size = 0;
  Properties->comm_rank = 0;
  Properties->num_grids = 0;

}

void ovkGetDomainPropertyName(const ovk_domain_properties *Properties, char *Name) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  strcpy(Name, Properties->name);

}

void ovkGetDomainPropertyDimension(const ovk_domain_properties *Properties, int *NumDims) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  *NumDims = Properties->num_dims;

}

void ovkGetDomainPropertyComm(const ovk_domain_properties *Properties, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Properties->comm;

}

void ovkGetDomainPropertyCommSize(const ovk_domain_properties *Properties, int *CommSize) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  *CommSize = Properties->comm_size;

}

void ovkGetDomainPropertyCommRank(const ovk_domain_properties *Properties, int *CommRank) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  *CommRank = Properties->comm_rank;

}

void ovkGetDomainPropertyGridCount(const ovk_domain_properties *Properties, int *NumGrids) {

  OVK_DEBUG_ASSERT(Properties, "Invalid properties pointer.");
  OVK_DEBUG_ASSERT(NumGrids, "Invalid num grids pointer.");

  *NumGrids = Properties->num_grids;

}

// void ovkGetDomainGridIDs(const ovk_domain *Domain, int *GridIDs) {



//   int iGrid;

//   iGrid = 0;
//   t_ordered_map_entry *Entry = OMBegin(Domain->grids);
//   while (Entry != OMEnd(Domain->grids)) {
//     GridIDs[iGrid] = OMKey(Entry);
//     ++iGrid;
//     Entry = OMNext(Entry);
//   }

// }
