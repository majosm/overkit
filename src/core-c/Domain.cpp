// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Domain.h"

#include "ovk/core-c/ConnectivityComponent.h"
#include "ovk/core-c/Context.h"
#include "ovk/core-c/GeometryComponent.h"
#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core-c/OverlapComponent.h"
#include "ovk/core-c/StateComponent.h"
#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/ConnectivityComponent.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/GeometryComponent.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/ID.hpp"
#include "ovk/core/Optional.hpp"
#include "ovk/core/OverlapComponent.hpp"
#include "ovk/core/StateComponent.hpp"

#include <mpi.h>

#include <cstring>
#include <memory>
#include <string>
#include <utility>

extern "C" {

void ovkCreateDomain(ovk_domain **Domain, ovk_shared_context *Context, ovk_domain_params **Params) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto &ContextCPP = *reinterpret_cast<std::shared_ptr<ovk::context> *>(Context);
  auto ParamsCPPPtr = reinterpret_cast<ovk::domain::params *>(*Params);

  auto DomainCPPPtr = new ovk::domain(ovk::CreateDomain(ContextCPP, std::move(*ParamsCPPPtr)));

  delete ParamsCPPPtr;

  *Domain = reinterpret_cast<ovk_domain *>(DomainCPPPtr);
  *Params = nullptr;

}

void ovkDestroyDomain(ovk_domain **Domain) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(*Domain, "Invalid domain pointer.");

  auto DomainCPPPtr = reinterpret_cast<ovk::domain *>(*Domain);

  delete DomainCPPPtr;

  *Domain = nullptr;

}

void ovkGetDomainContextC(const ovk_domain *Domain, const ovk_context **Context) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *Context = reinterpret_cast<const ovk_context *>(&DomainCPP.Context());

}

void ovkGetDomainContext(ovk_domain *Domain, ovk_context **Context) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  *Context = reinterpret_cast<ovk_context *>(&DomainCPP.Context());

}

void ovkGetDomainSharedContext(ovk_domain *Domain, ovk_shared_context **Context) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  auto &ContextCPP = DomainCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

void ovkGetDomainName(const ovk_domain *Domain, char *Name) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  std::strcpy(Name, DomainCPP.Name().c_str());

}

void ovkGetDomainDimension(const ovk_domain *Domain, int *NumDims) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *NumDims = DomainCPP.Dimension();

}

void ovkGetDomainComm(const ovk_domain *Domain, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *Comm = DomainCPP.Comm();

}

void ovkGetDomainCommSize(const ovk_domain *Domain, int *CommSize) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *CommSize = DomainCPP.Comm().Size();

}

void ovkGetDomainCommRank(const ovk_domain *Domain, int *CommRank) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *CommRank = DomainCPP.Comm().Rank();

}

void ovkGetGridCount(const ovk_domain *Domain, int *NumGrids) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(NumGrids, "Invalid num grids pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *NumGrids = DomainCPP.GridCount();

}

void ovkGetGridIDs(const ovk_domain *Domain, int *GridIDs) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridIDs, "Invalid grid IDs pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);

  auto &GridIDsCPP = DomainCPP.GridIDs();
  for (int iGrid = 0; iGrid < GridIDsCPP.Count(); ++iGrid) {
    GridIDs[iGrid] = GridIDsCPP[iGrid];
  }

}

bool ovkGridExists(const ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  return DomainCPP.GridExists(GridID);

}

void ovkGetNextAvailableGridID(const ovk_domain *Domain, int *GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *GridID = ovk::NextAvailableID(DomainCPP.GridIDs());

}

void ovkCreateGrid(ovk_domain *Domain, int GridID, ovk_grid_params **MaybeParams) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);

  ovk::grid::params *ParamsCPPPtr = nullptr;
  if (MaybeParams && *MaybeParams) {
    ParamsCPPPtr = reinterpret_cast<ovk::grid::params *>(*MaybeParams);
  }

  ovk::optional<ovk::grid::params> MaybeParamsCPP;
  if (ParamsCPPPtr) {
    MaybeParamsCPP = std::move(*ParamsCPPPtr);
  }

  DomainCPP.CreateGrid(GridID, std::move(MaybeParamsCPP));

  if (ParamsCPPPtr) {
    delete ParamsCPPPtr;
    *MaybeParams = nullptr;
  }

}

void ovkCreateGrids(ovk_domain *Domain, int Count, const int *GridIDs, ovk_grid_params
  **MaybeParams) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");
  OVK_DEBUG_ASSERT(GridIDs || Count == 0, "Invalid grid IDs pointer.");
  OVK_DEBUG_ASSERT(MaybeParams || Count == 0, "Invalid maybe params pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);

  ovk::array<ovk::grid::params *> ParamsCPPPtrs({Count}, nullptr);
  for (int iCreate = 0; iCreate < Count; ++iCreate) {
    if (MaybeParams[iCreate]) {
      ParamsCPPPtrs(iCreate) = reinterpret_cast<ovk::grid::params *>(MaybeParams[iCreate]);
    }
  }

  ovk::array<ovk::optional<ovk::grid::params>> MaybeParamsCPP({Count});
  for (int iCreate = 0; iCreate < Count; ++iCreate) {
    if (ParamsCPPPtrs(iCreate)) {
      MaybeParamsCPP(iCreate) = std::move(*ParamsCPPPtrs(iCreate));
    }
  }

  DomainCPP.CreateGrids({GridIDs, {Count}}, std::move(MaybeParamsCPP));

  for (int iCreate = 0; iCreate < Count; ++iCreate) {
    if (ParamsCPPPtrs(iCreate)) {
      delete ParamsCPPPtrs(iCreate);
      MaybeParams[iCreate] = nullptr;
    }
  }

}

void ovkDestroyGrid(ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  DomainCPP.DestroyGrid(GridID);

}

void ovkDestroyGrids(ovk_domain *Domain, int Count, const int *GridIDs) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");
  OVK_DEBUG_ASSERT(GridIDs || Count == 0, "Invalid grid IDs pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);

  DomainCPP.DestroyGrids({GridIDs, {Count}});

}

void ovkGetGridInfo(const ovk_domain *Domain, int GridID, const ovk_grid_info **GridInfo) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridInfo, "Invalid grid info pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *GridInfo = reinterpret_cast<const ovk_grid_info *>(&DomainCPP.GridInfo(GridID));

}

void ovkGetLocalGridCount(const ovk_domain *Domain, int *NumLocalGrids) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(NumLocalGrids, "Invalid num local grids pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *NumLocalGrids = DomainCPP.LocalGridCount();

}

void ovkGetLocalGridIDs(const ovk_domain *Domain, int *LocalGridIDs) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(LocalGridIDs, "Invalid grid IDs pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);

  auto &LocalGridIDsCPP = DomainCPP.LocalGridIDs();
  for (int iLocalGrid = 0; iLocalGrid < LocalGridIDsCPP.Count(); ++iLocalGrid) {
    LocalGridIDs[iLocalGrid] = LocalGridIDsCPP[iLocalGrid];
  }

}

bool ovkGridIsLocal(const ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  return DomainCPP.GridIsLocal(GridID);

}

void ovkGetGrid(const ovk_domain *Domain, int GridID, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *Grid = reinterpret_cast<const ovk_grid *>(&DomainCPP.Grid(GridID));

}

bool ovkComponentExists(const ovk_domain *Domain, int ComponentID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  return DomainCPP.ComponentExists(ComponentID);

}

void ovkGetNextAvailableComponentID(const ovk_domain *Domain, int *ComponentID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(ComponentID, "Invalid num components pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  *ComponentID = ovk::NextAvailableID(DomainCPP.ComponentIDs());

}

}

namespace {

template <typename ComponentType, typename ComponentTypeCPP, typename ParamsType, typename
  ParamsTypeCPP> void CreateComponent(ovk_domain *Domain, int ComponentID, void *ParamsVoid) {

  auto Params = static_cast<ParamsType **>(ParamsVoid);

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  ParamsTypeCPP *ParamsCPPPtr = nullptr;
  if (Params && *Params) {
    ParamsCPPPtr = reinterpret_cast<ParamsTypeCPP *>(*Params);
  }

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  if (ParamsCPPPtr) {
    DomainCPP.CreateComponent<ComponentTypeCPP>(ComponentID, std::move(*ParamsCPPPtr));
  } else {
    DomainCPP.CreateComponent<ComponentTypeCPP>(ComponentID);
  }

  if (ParamsCPPPtr) {
    delete ParamsCPPPtr;
    *Params = NULL;
  }

}

}

extern "C" {

void ovkCreateComponent(ovk_domain *Domain, int ComponentID, ovk_component_type ComponentType,
  void *Params) {

  OVK_DEBUG_ASSERT(ovkValidComponentType(ComponentType), "Invalid component type.");

  switch (ComponentType) {
  case OVK_COMPONENT_TYPE_GEOMETRY:
    CreateComponent<ovk_geometry_component, ovk::geometry_component, ovk_geometry_component_params,
      ovk::geometry_component::params>(Domain, ComponentID, Params);
    break;
  case OVK_COMPONENT_TYPE_STATE:
    CreateComponent<ovk_state_component, ovk::state_component, ovk_state_component_params,
      ovk::state_component::params>(Domain, ComponentID, Params);
    break;
  case OVK_COMPONENT_TYPE_OVERLAP:
    CreateComponent<ovk_overlap_component, ovk::overlap_component, ovk_overlap_component_params,
      ovk::overlap_component::params>(Domain, ComponentID, Params);
    break;
  case OVK_COMPONENT_TYPE_CONNECTIVITY:
    CreateComponent<ovk_connectivity_component, ovk::connectivity_component,
      ovk_connectivity_component_params, ovk::connectivity_component::params>(Domain, ComponentID,
        Params);
    break;
  }

}

void ovkDestroyComponent(ovk_domain *Domain, int ComponentID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  DomainCPP.DestroyComponent(ComponentID);

}

}

namespace {

template <typename ComponentType, typename ComponentTypeCPP> void GetComponent(const ovk_domain
  *Domain, int ComponentID, void *ComponentVoid) {

  auto Component = static_cast<const ComponentType **>(ComponentVoid);

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Component, "Invalid component pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  auto &ComponentCPP = DomainCPP.Component<ComponentTypeCPP>(ComponentID);

  *Component = reinterpret_cast<const ComponentType *>(&ComponentCPP);

}

}

extern "C" {

void ovkGetComponent(const ovk_domain *Domain, int ComponentID, ovk_component_type ComponentType,
  void *Component) {

  OVK_DEBUG_ASSERT(ovkValidComponentType(ComponentType), "Invalid component type.");

  switch (ComponentType) {
  case OVK_COMPONENT_TYPE_GEOMETRY:
    GetComponent<ovk_geometry_component, ovk::geometry_component>(Domain, ComponentID, Component);
    break;
  case OVK_COMPONENT_TYPE_STATE:
    GetComponent<ovk_state_component, ovk::state_component>(Domain, ComponentID, Component);
    break;
  case OVK_COMPONENT_TYPE_OVERLAP:
    GetComponent<ovk_overlap_component, ovk::overlap_component>(Domain, ComponentID, Component);
    break;
  case OVK_COMPONENT_TYPE_CONNECTIVITY:
    GetComponent<ovk_connectivity_component, ovk::connectivity_component>(Domain, ComponentID,
      Component);
    break;
  }

}

}

namespace {

template <typename ComponentType, typename ComponentTypeCPP> void EditComponent(ovk_domain *Domain,
  int ComponentID, void *ComponentVoid) {

  auto Component = static_cast<ComponentType **>(ComponentVoid);

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Component, "Invalid component pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);

  ovk::edit_handle<ComponentTypeCPP> EditHandle = DomainCPP.EditComponent<ComponentTypeCPP>(
    ComponentID);
  auto ComponentCPPPtr = EditHandle.Release();

  *Component = reinterpret_cast<ComponentType *>(ComponentCPPPtr);

}

}

extern "C" {

void ovkEditComponent(ovk_domain *Domain, int ComponentID, ovk_component_type ComponentType,
  void *Component) {

  OVK_DEBUG_ASSERT(ovkValidComponentType(ComponentType), "Invalid component type.");

  switch (ComponentType) {
  case OVK_COMPONENT_TYPE_GEOMETRY:
    EditComponent<ovk_geometry_component, ovk::geometry_component>(Domain, ComponentID, Component);
    break;
  case OVK_COMPONENT_TYPE_STATE:
    EditComponent<ovk_state_component, ovk::state_component>(Domain, ComponentID, Component);
    break;
  case OVK_COMPONENT_TYPE_OVERLAP:
    EditComponent<ovk_overlap_component, ovk::overlap_component>(Domain, ComponentID, Component);
    break;
  case OVK_COMPONENT_TYPE_CONNECTIVITY:
    EditComponent<ovk_connectivity_component, ovk::connectivity_component>(Domain, ComponentID,
      Component);
    break;
  }

}

}

namespace {

template <typename ComponentType, typename ComponentTypeCPP> void RestoreComponent(ovk_domain
  *Domain, int ComponentID, void *ComponentVoid) {

  auto Component = static_cast<ComponentType **>(ComponentVoid);

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Component, "Invalid component pointer.");
  OVK_DEBUG_ASSERT(*Component, "Invalid component pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  DomainCPP.RestoreComponent(ComponentID);

  *Component = nullptr;

}

}

extern "C" {

void ovkRestoreComponent(ovk_domain *Domain, int ComponentID, ovk_component_type ComponentType, void
  *Component) { 

  OVK_DEBUG_ASSERT(ovkValidComponentType(ComponentType), "Invalid component type.");

  switch (ComponentType) {
  case OVK_COMPONENT_TYPE_GEOMETRY:
    RestoreComponent<ovk_geometry_component, ovk::geometry_component>(Domain, ComponentID,
      Component);
    break;
  case OVK_COMPONENT_TYPE_STATE:
    RestoreComponent<ovk_state_component, ovk::state_component>(Domain, ComponentID, Component);
    break;
  case OVK_COMPONENT_TYPE_OVERLAP:
    RestoreComponent<ovk_overlap_component, ovk::overlap_component>(Domain, ComponentID, Component);
    break;
  case OVK_COMPONENT_TYPE_CONNECTIVITY:
    RestoreComponent<ovk_connectivity_component, ovk::connectivity_component>(Domain, ComponentID,
      Component);
    break;
  }

}

void ovkCreateDomainParams(ovk_domain_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::domain::params();

  *Params = reinterpret_cast<ovk_domain_params *>(ParamsCPPPtr);

}

void ovkDestroyDomainParams(ovk_domain_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::domain::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetDomainParamName(const ovk_domain_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::domain::params *>(Params);
  std::strcpy(Name, ParamsCPP.Name().c_str());

}

void ovkSetDomainParamName(ovk_domain_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::domain::params *>(Params);
  ParamsCPP.SetName(Name);

}

void ovkGetDomainParamDimension(const ovk_domain_params *Params, int *NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::domain::params *>(Params);
  *NumDims = ParamsCPP.Dimension();

}

void ovkSetDomainParamDimension(ovk_domain_params *Params, int NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::domain::params *>(Params);
  ParamsCPP.SetDimension(NumDims);

}

void ovkGetDomainParamComm(const ovk_domain_params *Params, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::domain::params *>(Params);
  *Comm = ParamsCPP.Comm();

}

void ovkSetDomainParamComm(ovk_domain_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::domain::params *>(Params);
  ParamsCPP.SetComm(Comm);

}

}
