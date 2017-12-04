// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DOMAIN_INCLUDED
#define OVK_CORE_DOMAIN_INCLUDED

#include "ovkDomain.h"

#include "ErrorHandler.h"
#include "Global.h"
#include "Grid.h"
#include "Logger.h"
#include "OrderedMap.h"

struct ovk_domain_params {
  char name[OVK_NAME_LENGTH];
  int num_dims;
  MPI_Comm comm;
};

struct ovk_domain {
  ovk_domain_properties *properties;
  t_logger *logger;
  t_error_handler *error_handler;
  ovk_domain_config config;
  t_ordered_map *grids;
};

struct ovk_domain_properties {
  char name[OVK_NAME_LENGTH];
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  int num_grids;
};

void CreateDomainParams(ovk_domain_params **Params, MPI_Comm DefaultComm);
void DestroyDomainParams(ovk_domain_params **Params);

void CreateDomain(ovk_domain **Domain, const ovk_domain_params *Params, t_logger *Logger,
  t_error_handler *ErrorHandler);
void DestroyDomain(ovk_domain **Domain);

#endif
