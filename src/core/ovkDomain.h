// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_DOMAIN_INCLUDED
#define OVK_CORE_PUBLIC_DOMAIN_INCLUDED

#include <ovkGlobal.h>
#include <ovkGrid.h>

struct ovk_domain_params;
typedef struct ovk_domain_params ovk_domain_params;

struct ovk_domain;
typedef struct ovk_domain ovk_domain;

struct ovk_domain_properties;
typedef struct ovk_domain_properties ovk_domain_properties;

void ovkConfigureDomain(ovk_domain *Domain, ovk_domain_config Config);

void ovkGetDomainProperties(const ovk_domain *Domain, const ovk_domain_properties **Properties);
// void ovkEditDomainProperties(ovk_domain *Domain, ovk_domain_properties **Properties);
// void ovkReleaseDomainProperties(ovk_domain *Domain, ovk_domain_properties **Properties);

void ovkCreateGridParams(ovk_domain *Domain, ovk_grid_params **Params);
void ovkDestroyGridParams(ovk_domain *Domain, ovk_grid_params **Params);

void ovkCreateGridLocal(ovk_domain *Domain, int *GridID, const ovk_grid_params *Params);
void ovkCreateGridRemote(ovk_domain *Domain, int *GridID);
void ovkDestroyGrid(ovk_domain *Domain, int *GridID);

void ovkGetGrid(const ovk_domain *Domain, int GridID, const ovk_grid **Grid);
void ovkGetGrids(const ovk_domain *Domain, int Count, int *GridIDs, const ovk_grid **Grids);

void ovkGetDomainParamDimension(const ovk_domain_params *Params, int *NumDims);
void ovkSetDomainParamDimension(ovk_domain_params *Params, int NumDims);
void ovkGetDomainParamComm(const ovk_domain_params *Params, MPI_Comm *Comm);
void ovkSetDomainParamComm(ovk_domain_params *Params, MPI_Comm Comm);

void ovkGetDomainPropertyDimension(const ovk_domain_properties *Properties, int *NumDims);
void ovkGetDomainPropertyComm(const ovk_domain_properties *Properties, MPI_Comm *Comm);

#endif
