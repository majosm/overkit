// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_DOMAIN_H_INCLUDED
#define OVK_CORE_C_DOMAIN_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core-c/Grid.h>
#include <ovk/core-c/Request.h>
#include <ovk/core/Domain.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_COMPONENT_TYPE_GEOMETRY = 0,
  OVK_COMPONENT_TYPE_STATE,
  OVK_COMPONENT_TYPE_CONNECTIVITY
} ovk_component_type;

static inline bool ovkValidComponentType(ovk_component_type ComponentType) {

  switch (ComponentType) {
  case OVK_COMPONENT_TYPE_GEOMETRY:
  case OVK_COMPONENT_TYPE_STATE:
  case OVK_COMPONENT_TYPE_CONNECTIVITY:
    return true;
  default:
    return false;
  }

}

struct ovk_domain;
typedef struct ovk_domain ovk_domain;

struct ovk_domain_params;
typedef struct ovk_domain_params ovk_domain_params;

void ovkCreateDomain(ovk_domain **Domain, ovk_shared_context *Context, ovk_domain_params **Params);
void ovkDestroyDomain(ovk_domain **Domain);

void ovkGetDomainContextC(const ovk_domain *Domain, const ovk_context **Context);
void ovkGetDomainContext(ovk_domain *Domain, ovk_context **Context);
void ovkGetDomainSharedContext(ovk_domain *Domain, ovk_shared_context **Context);

void ovkGetDomainName(const ovk_domain *Domain, char *Name);
void ovkGetDomainDimension(const ovk_domain *Domain, int *NumDims);
void ovkGetDomainComm(const ovk_domain *Domain, MPI_Comm *Comm);
void ovkGetDomainCommSize(const ovk_domain *Domain, int *CommSize);
void ovkGetDomainCommRank(const ovk_domain *Domain, int *CommRank);

void ovkGetGridCount(const ovk_domain *Domain, int *NumGrids);

void ovkGetGridIDs(const ovk_domain *Domain, int *GridIDs);

bool ovkGridExists(const ovk_domain *Domain, int GridID);
void ovkGetNextAvailableGridID(const ovk_domain *Domain, int *GridID);

void ovkCreateGrid(ovk_domain *Domain, int GridID, ovk_grid_params **MaybeParams);
void ovkCreateGrids(ovk_domain *Domain, int Count, const int *GridIDs, ovk_grid_params
  **MaybeParams);

void ovkDestroyGrid(ovk_domain *Domain, int GridID);
void ovkDestroyGrids(ovk_domain *Domain, int Count, const int *GridIDs);

void ovkGetGridInfo(const ovk_domain *Domain, int GridID, const ovk_grid_info **GridInfo);

void ovkGetLocalGridCount(const ovk_domain *Domain, int *NumLocalGrids);
void ovkGetLocalGridIDs(const ovk_domain *Domain, int *LocalGridIDs);

bool ovkGridIsLocal(const ovk_domain *Domain, int GridID);
void ovkGetGrid(const ovk_domain *Domain, int GridID, const ovk_grid **Grid);

bool ovkComponentExists(const ovk_domain *Domain, int ComponentID);
void ovkGetNextAvailableComponentID(const ovk_domain *Domain, int *ComponentID);
bool ovkComponentHasType(const ovk_domain *Domain, int ComponentID, ovk_component_type
  ComponentType);

// "Params" actual type is <component-specific-params-type> **
void ovkCreateComponent(ovk_domain *Domain, int ComponentID, ovk_component_type ComponentType, void
  *Params);
void ovkDestroyComponent(ovk_domain *Domain, int ComponentID);

// "Component" actual type is const T **
void ovkGetComponent(const ovk_domain *Domain, int ComponentID, ovk_component_type ComponentType,
  void *Component);
// "Component" actual type is T **
void ovkEditComponent(ovk_domain *Domain, int ComponentID, ovk_component_type ComponentType,
  void *Component);
void ovkRestoreComponent(ovk_domain *Domain, int ComponentID, ovk_component_type ComponentType,
  void *Component);

void ovkCreateDomainParams(ovk_domain_params **Params);
void ovkDestroyDomainParams(ovk_domain_params **Params);
void ovkGetDomainParamName(const ovk_domain_params *Params, char *Name);
void ovkSetDomainParamName(ovk_domain_params *Params, const char *Name);
void ovkGetDomainParamDimension(const ovk_domain_params *Params, int *NumDims);
void ovkSetDomainParamDimension(ovk_domain_params *Params, int NumDims);
void ovkGetDomainParamComm(const ovk_domain_params *Params, MPI_Comm *Comm);
void ovkSetDomainParamComm(ovk_domain_params *Params, MPI_Comm Comm);

#ifdef __cplusplus
}
#endif

#endif
