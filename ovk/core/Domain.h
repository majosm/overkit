// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DOMAIN_INCLUDED
#define OVK_CORE_DOMAIN_INCLUDED

#include "ovk/core/ovkDomain.h"

#include "ovk/core/Connectivity.h"
#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Exchange.h"
#include "ovk/core/Global.h"
#include "ovk/core/Grid.h"
#include "ovk/core/Logger.h"
#include "ovk/core/OrderedMap.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_domain_params {
  char name[OVK_NAME_LENGTH];
  int num_dims;
  MPI_Comm comm;
};

struct ovk_domain_properties {
  char name[OVK_NAME_LENGTH];
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  int num_grids;
};

typedef struct {
  ovk_grid_info *info;
  int edit_ref_count;
  bool has_local_data;
  ovk_grid *grid;
} t_domain_grid_container;

typedef struct {
  ovk_connectivity_info *info;
  int edit_ref_count;
  bool has_local_data;
  ovk_connectivity *connectivity;
} t_domain_connectivity_container;

typedef struct {
  ovk_exchange_info *info;
  bool has_local_data;
  ovk_exchange *exchange;
} t_domain_exchange_container;

struct ovk_domain {
  ovk_domain_properties properties;
  int properties_edit_ref_count;
  t_logger *logger;
  t_error_handler *error_handler;
  ovk_domain_config config;
  t_ordered_map *grids;
  int grids_edit_ref_count;
  t_ordered_map *connectivities;
  int connectivities_edit_ref_count;
  t_ordered_map *exchanges;
};

void PRIVATE(CreateDomain)(ovk_domain **Domain, const ovk_domain_params *Params,
  t_logger *Logger, t_error_handler *ErrorHandler);
#define CreateDomain(...) PRIVATE(CreateDomain)(__VA_ARGS__)
void PRIVATE(DestroyDomain)(ovk_domain **Domain);
#define DestroyDomain(...) PRIVATE(DestroyDomain)(__VA_ARGS__)

void PRIVATE(CreateDomainParams)(ovk_domain_params **Params, MPI_Comm DefaultComm);
#define CreateDomainParams(...) PRIVATE(CreateDomainParams)(__VA_ARGS__)
void PRIVATE(DestroyDomainParams)(ovk_domain_params **Params);
#define DestroyDomainParams(...) PRIVATE(DestroyDomainParams)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
