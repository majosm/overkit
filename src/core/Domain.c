// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "Domain.h"

#include "Debug.h"
#include "ErrorHandler.h"
#include "Global.h"
#include "Logger.h"
#include "OrderedMap.h"
#include "ParallelUtils.h"
#include "TextUtils.h"

typedef struct {
  bool local;
  ovk_grid *local_data;
  int root_rank;
} t_grid_container;

static void DefaultDomainProperties(ovk_domain_properties *Properties);

void CreateDomain(ovk_domain **Domain_, const ovk_domain_params *Params, t_logger *Logger,
  t_error_handler *ErrorHandler) {

  OVK_DEBUG_ASSERT(Params, "Invalid domain params pointer.");

  MPI_Comm Comm;
  MPI_Comm_dup(Params->comm, &Comm);

  MPI_Barrier(Comm);

  *Domain_ = malloc(sizeof(ovk_domain));
  ovk_domain *Domain = *Domain_;

  DefaultDomainProperties(&Domain->properties);

  if (strlen(Params->name) > 0) {
    strncpy(Domain->properties.name, Params->name, OVK_NAME_LENGTH);
  } else {
    strcpy(Domain->properties.name, "Domain");
  }

  Domain->properties.num_dims = Params->num_dims;

  Domain->properties.comm = Comm;
  MPI_Comm_size(Domain->properties.comm, &Domain->properties.comm_size);
  MPI_Comm_rank(Domain->properties.comm, &Domain->properties.comm_rank);

  Domain->logger = Logger;
  Domain->error_handler = ErrorHandler;

  Domain->config = OVK_DOMAIN_CONFIG_NONE;

  OMCreate(&Domain->grids);

  if (Domain->properties.comm_rank == 0) {
    char ProcessesString[32];
    PluralizeLabel(Domain->properties.comm_size, "processes", "process", ProcessesString);
    LogStatus(Logger, true, 0, "Created %1iD domain '%s' on %s.", Domain->properties.num_dims,
      Domain->properties.name, ProcessesString);
  }

  MPI_Barrier(Domain->properties.comm);

}

void DestroyDomain(ovk_domain **Domain_) {

  ovk_domain *Domain = *Domain_;

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  MPI_Barrier(Domain->properties.comm);

  t_ordered_map_entry *Entry = OMBegin(Domain->grids);
  while (Entry != OMEnd(Domain->grids)) {
    t_grid_container *Container = OMData(Entry);
    if (Container->local) {
      ovk_grid *Grid = Container->local_data;
      DestroyGrid(&Grid);
    }
    free(Container);
    Entry = OMNext(Entry);
  }
  OMDestroy(&Domain->grids);

  t_logger *Logger = Domain->logger;
  MPI_Comm Comm = Domain->properties.comm;
  bool IsDomainRoot = Domain->properties.comm_rank == 0;
  char Name[OVK_NAME_LENGTH];
  strncpy(Name, Domain->properties.name, OVK_NAME_LENGTH);

  free(*Domain_);
  *Domain_ = NULL;

  MPI_Barrier(Comm);

  MPI_Comm_free(&Comm);

  LogStatus(Logger, IsDomainRoot, 0, "Destroyed domain '%s'.", Name);

}

void ovkConfigureDomain(ovk_domain *Domain, ovk_domain_config Config) {

  Domain->config = Config;

}

void ovkGetDomainProperties(const ovk_domain *Domain, const ovk_domain_properties **Properties) {

  *Properties = &Domain->properties;

}

void ovkCreateGridParams(ovk_domain *Domain, ovk_grid_params **Params) {

  CreateGridParams(Params, Domain->properties.num_dims, Domain->properties.comm);

}

void ovkDestroyGridParams(ovk_domain *Domain, ovk_grid_params **Params) {

  OVK_DEBUG_ASSERT(*Params, "Invalid grid params pointer.");

  DestroyGridParams(Params);

}

void ovkCreateGridLocal(ovk_domain *Domain, int *GridID_, const ovk_grid_params *Params) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid grid params pointer.");

  MPI_Barrier(Domain->properties.comm);

  MPI_Comm TempComm;
  MPI_Comm_dup(Params->comm, &TempComm);
  int GridCommRank;
  MPI_Comm_rank(TempComm, &GridCommRank);
  MPI_Comm_free(&TempComm);

  bool IsGridRoot = GridCommRank == 0;

  int Local = 1;
  int AtLeastOneLocal = 0;
  MPI_Allreduce(&Local, &AtLeastOneLocal, 1, MPI_INT, MPI_LOR, Domain->properties.comm);

  int GridID;
  if (IsGridRoot) GridID = Params->id;
  BroadcastAnySource(&GridID, 1, MPI_INT, IsGridRoot, Domain->properties.comm);

  if (GridID == OVK_GENERATE_ID) {
    GridID = OMNextAvailableKey(Domain->grids);
  }

  ovk_grid *Grid;
  CreateGrid(&Grid, GridID, Params, Domain->logger, Domain->error_handler);

  int GridRootRank;
  if (IsGridRoot) GridRootRank = Domain->properties.comm_rank;
  BroadcastAnySource(&GridRootRank, 1, MPI_INT, IsGridRoot, Domain->properties.comm);

  t_grid_container *Container = malloc(sizeof(t_grid_container));
  Container->local = true;
  Container->local_data = Grid;
  Container->root_rank = GridRootRank;

  OMInsert(Domain->grids, GridID, Container);

  ++Domain->properties.num_grids;

  if (GridID_) {
    *GridID_ = GridID;
  }

  MPI_Barrier(Domain->properties.comm);

}

void ovkCreateGridRemote(ovk_domain *Domain, int *GridID_) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  MPI_Barrier(Domain->properties.comm);

  int Local = 0;
  int AtLeastOneLocal = 0;
  MPI_Allreduce(&Local, &AtLeastOneLocal, 1, MPI_INT, MPI_LOR, Domain->properties.comm);
  OVK_DEBUG_ASSERT(AtLeastOneLocal, "Grid must exist on at least one process.");

  int GridID;
  BroadcastAnySource(&GridID, 1, MPI_INT, false, Domain->properties.comm);

  if (GridID == OVK_GENERATE_ID) {
    GridID = OMNextAvailableKey(Domain->grids);
  }

  int GridRootRank;
  BroadcastAnySource(&GridRootRank, 1, MPI_INT, false, Domain->properties.comm);

  t_grid_container *Container = malloc(sizeof(t_grid_container));
  Container->local = false;
  Container->local_data = NULL;
  Container->root_rank = GridRootRank;

  OMInsert(Domain->grids, GridID, Container);

  ++Domain->properties.num_grids;

  if (GridID_) {
    *GridID_ = GridID;
  }

  MPI_Barrier(Domain->properties.comm);

}

void ovkDestroyGrid(ovk_domain *Domain, int *GridID_) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID_, "Invalid grid ID pointer.");

  MPI_Barrier(Domain->properties.comm);

  int GridID = *GridID_;

  t_ordered_map_entry *Entry = OMFind(Domain->grids, GridID);

  OVK_DEBUG_ASSERT(Entry != OMEnd(Domain->grids), "Grid %i does not exist.", GridID);

  t_grid_container *Container = OMData(Entry);

  if (Container->local) {
    ovk_grid *Grid = Container->local_data;
    DestroyGrid(&Grid);
  }

  free(Container);

  OMRemove(Domain->grids, GridID);

  --Domain->properties.num_grids;

  *GridID_ = -1;

  MPI_Barrier(Domain->properties.comm);

}

void ovkGetGrid(const ovk_domain *Domain, int GridID, const ovk_grid **Grid) {

  const t_ordered_map_entry *Entry = OMFindC(Domain->grids, GridID);

  OVK_DEBUG_ASSERT(Entry != OMEnd(Domain->grids), "Grid %i does not exist.", GridID);

  const t_grid_container *Container = OMDataC(Entry);

  OVK_DEBUG_ASSERT(Container->local, "Grid %i does not have local data on rank %i.", GridID,
    Domain->properties.comm_rank);

  if (Container->local) {
    *Grid = Container->local_data;
  } else {
    *Grid = NULL;
  }

}

void ovkGetGrids(const ovk_domain *Domain, int Count, int *GridIDs, const ovk_grid **Grids) {

  int i;

  OVK_DEBUG_ASSERT(Count > 0 && Count < OMSize(Domain->grids), "Invalid count.");

  for (i = 0; i < Count; ++i) {
    ovkGetGrid(Domain, GridIDs[i], &Grids[i]);
  }

}

void CreateDomainParams(ovk_domain_params **Params_, MPI_Comm DefaultComm) {

  *Params_ = malloc(sizeof(ovk_domain_params));
  ovk_domain_params *Params = *Params_;

  memset(Params->name, 0, OVK_NAME_LENGTH);

  Params->num_dims = 2;
  Params->comm = DefaultComm;

}

void DestroyDomainParams(ovk_domain_params **Params) {

  free(*Params);
  *Params = NULL;

}

void ovkGetDomainParamName(const ovk_domain_params *Params, char *Name) {

  strcpy(Name, Params->name);

}

void ovkSetDomainParamName(ovk_domain_params *Params, const char *Name) {

  strncpy(Params->name, Name, OVK_NAME_LENGTH);

}

void ovkGetDomainParamDimension(const ovk_domain_params *Params, int *NumDims) {

  *NumDims = Params->num_dims;

}

void ovkSetDomainParamDimension(ovk_domain_params *Params, int NumDims) {

  OVK_DEBUG_ASSERT(NumDims == 2 || NumDims == 3, "Invalid dimension.");

  Params->num_dims = NumDims;

}

void ovkGetDomainParamComm(const ovk_domain_params *Params, MPI_Comm *Comm) {

  *Comm = Params->comm;

}

void ovkSetDomainParamComm(ovk_domain_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Params->comm = Comm;

}

static void DefaultDomainProperties(ovk_domain_properties *Properties) {

  memset(Properties->name, 0, OVK_NAME_LENGTH);

  Properties->num_dims = 2;
  Properties->comm = MPI_COMM_NULL;
  Properties->comm_size = 0;
  Properties->comm_rank = 0;
  Properties->num_grids = 0;

}

void ovkGetDomainPropertyName(const ovk_domain_properties *Properties, char *Name) {

  strcpy(Name, Properties->name);

}

void ovkGetDomainPropertyDimension(const ovk_domain_properties *Properties, int *NumDims) {

  *NumDims = Properties->num_dims;

}

void ovkGetDomainPropertyComm(const ovk_domain_properties *Properties, MPI_Comm *Comm) {

  *Comm = Properties->comm;

}

void ovkGetDomainPropertyGridCount(const ovk_domain_properties *Properties, int *NumGrids) {

  *NumGrids = Properties->num_grids;

}
