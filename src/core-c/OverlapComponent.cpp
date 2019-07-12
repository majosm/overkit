// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/OverlapComponent.h"

#include "ovk/core-c/Global.h"
#include "ovk/core-c/OverlapM.h"
#include "ovk/core-c/OverlapN.h"
#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/OverlapComponent.hpp"
#include "ovk/core/OverlapM.hpp"
#include "ovk/core/OverlapN.hpp"

#include <mpi.h>

#include <utility>

extern "C" {

int ovkOverlapCount(const ovk_overlap_component *OverlapComponent) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(OverlapComponent);
  return OverlapComponentCPP.OverlapCount();

}

bool ovkOverlapExists(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(OverlapComponent);
  return OverlapComponentCPP.OverlapExists(MGridID, NGridID);

}

void ovkCreateOverlap(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<ovk::overlap_component *>(OverlapComponent);
  OverlapComponentCPP.CreateOverlap(MGridID, NGridID);

}

void ovkCreateOverlaps(ovk_overlap_component *OverlapComponent, int Count, const int *MGridIDs,
  const int *NGridIDs) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");

  auto &OverlapComponentCPP = *reinterpret_cast<ovk::overlap_component *>(OverlapComponent);
  OverlapComponentCPP.CreateOverlaps({MGridIDs, {Count}}, {NGridIDs, {Count}});

}

void ovkDestroyOverlap(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<ovk::overlap_component *>(OverlapComponent);
  OverlapComponentCPP.DestroyOverlap(MGridID, NGridID);

}

void ovkDestroyOverlaps(ovk_overlap_component *OverlapComponent, int Count, const int *MGridIDs,
  const int *NGridIDs) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");

  auto &OverlapComponentCPP = *reinterpret_cast<ovk::overlap_component *>(OverlapComponent);
  OverlapComponentCPP.DestroyOverlaps({MGridIDs, {Count}}, {NGridIDs, {Count}});

}

int ovkLocalOverlapMCount(const ovk_overlap_component *OverlapComponent) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(OverlapComponent);
  return OverlapComponentCPP.LocalOverlapMCount();

}

int ovkLocalOverlapMCountForGrid(const ovk_overlap_component *OverlapComponent, int MGridID) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(OverlapComponent);
  return OverlapComponentCPP.LocalOverlapMCountForGrid(MGridID);

}

void ovkGetOverlapM(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID, const
  ovk_overlap_m **OverlapM) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");
  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(OverlapComponent);
  *OverlapM = reinterpret_cast<const ovk_overlap_m *>(&OverlapComponentCPP.OverlapM(MGridID,
    NGridID));

}

bool ovkEditingOverlapM(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(OverlapComponent);
  return OverlapComponentCPP.EditingOverlapM(MGridID, NGridID);

}

void ovkEditOverlapM(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID,
  ovk_overlap_m **OverlapM) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");
  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<ovk::overlap_component *>(OverlapComponent);

  ovk::edit_handle<ovk::overlap_m> EditHandle = OverlapComponentCPP.EditOverlapM(MGridID, NGridID);
  auto OverlapMCPPPtr = EditHandle.Release();

  *OverlapM = reinterpret_cast<ovk_overlap_m *>(OverlapMCPPPtr);

}

void ovkRestoreOverlapM(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID,
  ovk_overlap_m **OverlapM) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");
  OVK_DEBUG_ASSERT(OverlapM, "Invalid overlap M pointer.");
  OVK_DEBUG_ASSERT(*OverlapM, "Invalid overlap M pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<ovk::overlap_component *>(OverlapComponent);
  OverlapComponentCPP.RestoreOverlapM(MGridID, NGridID);

  *OverlapM = nullptr;

}

int ovkLocalOverlapNCount(const ovk_overlap_component *OverlapComponent) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(OverlapComponent);
  return OverlapComponentCPP.LocalOverlapNCount();

}

int ovkLocalOverlapNCountForGrid(const ovk_overlap_component *OverlapComponent, int NGridID) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(
    OverlapComponent);
  return OverlapComponentCPP.LocalOverlapNCountForGrid(NGridID);

}

void ovkGetOverlapN(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID, const
  ovk_overlap_n **OverlapN) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");
  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(OverlapComponent);
  *OverlapN = reinterpret_cast<const ovk_overlap_n *>(&OverlapComponentCPP
    .OverlapN(MGridID, NGridID));

}

bool ovkEditingOverlapN(const ovk_overlap_component *OverlapComponent, int MGridID, int NGridID) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<const ovk::overlap_component *>(OverlapComponent);
  return OverlapComponentCPP.EditingOverlapN(MGridID, NGridID);

}

void ovkEditOverlapN(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID,
  ovk_overlap_n **OverlapN) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");
  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<ovk::overlap_component *>(OverlapComponent);

  ovk::edit_handle<ovk::overlap_n> EditHandle = OverlapComponentCPP.EditOverlapN(MGridID, NGridID);
  auto OverlapNCPPPtr = EditHandle.Release();

  *OverlapN = reinterpret_cast<ovk_overlap_n *>(OverlapNCPPPtr);

}

void ovkRestoreOverlapN(ovk_overlap_component *OverlapComponent, int MGridID, int NGridID,
  ovk_overlap_n **OverlapN) {

  OVK_DEBUG_ASSERT(OverlapComponent, "Invalid overlap component pointer.");
  OVK_DEBUG_ASSERT(OverlapN, "Invalid overlap N pointer.");
  OVK_DEBUG_ASSERT(*OverlapN, "Invalid overlap N pointer.");

  auto &OverlapComponentCPP = *reinterpret_cast<ovk::overlap_component *>(OverlapComponent);
  OverlapComponentCPP.RestoreOverlapN(MGridID, NGridID);

  *OverlapN = nullptr;

}

void ovkCreateOverlapComponentParams(ovk_overlap_component_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::overlap_component::params();

  *Params = reinterpret_cast<ovk_overlap_component_params *>(ParamsCPPPtr);

}

void ovkDestroyOverlapComponentParams(ovk_overlap_component_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::overlap_component::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetOverlapComponentParamName(const ovk_overlap_component_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::overlap_component::params *>(Params);
  std::strcpy(Name, ParamsCPP.Name().c_str());

}

void ovkSetOverlapComponentParamName(ovk_overlap_component_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::overlap_component::params *>(Params);
  ParamsCPP.SetName(Name);

}

}
