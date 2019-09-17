// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/StateComponent.h"

#include "ovk/core-c/Global.h"
#include "ovk/core-c/State.h"
#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/State.hpp"
#include "ovk/core/StateComponent.hpp"

#include <mpi.h>

#include <utility>

extern "C" {

int ovkStateCount(const ovk_state_component *StateComponent) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");

  auto &StateComponentCPP = *reinterpret_cast<const ovk::state_component *>(StateComponent);
  return StateComponentCPP.StateCount();

}

bool ovkStateExists(const ovk_state_component *StateComponent, int GridID) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");

  auto &StateComponentCPP = *reinterpret_cast<const ovk::state_component *>(StateComponent);
  return StateComponentCPP.StateExists(GridID);

}

void ovkCreateState(ovk_state_component *StateComponent, int GridID, ovk_state_params
  **MaybeParams) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");

  auto &StateComponentCPP = *reinterpret_cast<ovk::state_component *>(StateComponent);

  ovk::state::params *ParamsCPPPtr = nullptr;
  if (MaybeParams && *MaybeParams) {
    ParamsCPPPtr = reinterpret_cast<ovk::state::params *>(*MaybeParams);
  }

  ovk::optional<ovk::state::params> MaybeParamsCPP;
  if (ParamsCPPPtr) {
    MaybeParamsCPP = std::move(*ParamsCPPPtr);
  }

  StateComponentCPP.CreateState(GridID, std::move(MaybeParamsCPP));

  if (ParamsCPPPtr) {
    delete ParamsCPPPtr;
    *MaybeParams = nullptr;
  }

}

void ovkCreateStates(ovk_state_component *StateComponent, int Count, const int *GridIDs,
  ovk_state_params **MaybeParams) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");

  auto &StateComponentCPP = *reinterpret_cast<ovk::state_component *>(StateComponent);

  if (MaybeParams) {

    ovk::array<ovk::state::params *> ParamsCPPPtrs({Count}, nullptr);
    for (int iCreate = 0; iCreate < Count; ++iCreate) {
      if (MaybeParams[iCreate]) {
        ParamsCPPPtrs(iCreate) = reinterpret_cast<ovk::state::params *>(MaybeParams[iCreate]);
      }
    }

    ovk::array<ovk::optional<ovk::state::params>> MaybeParamsCPP({Count});
    for (int iCreate = 0; iCreate < Count; ++iCreate) {
      if (ParamsCPPPtrs(iCreate)) {
        MaybeParamsCPP(iCreate) = std::move(*ParamsCPPPtrs(iCreate));
      }
    }

    StateComponentCPP.CreateStates({GridIDs, {Count}}, std::move(MaybeParamsCPP));

    for (int iCreate = 0; iCreate < Count; ++iCreate) {
      if (ParamsCPPPtrs(iCreate)) {
        delete ParamsCPPPtrs(iCreate);
        MaybeParams[iCreate] = nullptr;
      }
    }

  } else {

    StateComponentCPP.CreateStates({GridIDs, {Count}});

  }

}

void ovkDestroyState(ovk_state_component *StateComponent, int GridID) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");

  auto &StateComponentCPP = *reinterpret_cast<ovk::state_component *>(StateComponent);
  StateComponentCPP.DestroyState(GridID);

}

void ovkDestroyStates(ovk_state_component *StateComponent, int Count, const int *GridIDs) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");
  OVK_DEBUG_ASSERT(Count >= 0, "Invalid count value.");

  auto &StateComponentCPP = *reinterpret_cast<ovk::state_component *>(StateComponent);
  StateComponentCPP.DestroyStates({GridIDs, {Count}});

}

int ovkLocalStateCount(const ovk_state_component *StateComponent) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");

  auto &StateComponentCPP = *reinterpret_cast<const ovk::state_component *>(StateComponent);
  return StateComponentCPP.LocalStateCount();

}

void ovkGetState(const ovk_state_component *StateComponent, int GridID, const ovk_state **State) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");
  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");

  auto &StateComponentCPP = *reinterpret_cast<const ovk::state_component *>(StateComponent);
  *State = reinterpret_cast<const ovk_state *>(&StateComponentCPP.State(GridID));

}

bool ovkEditingState(const ovk_state_component *StateComponent, int GridID) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");

  auto &StateComponentCPP = *reinterpret_cast<const ovk::state_component *>(StateComponent);
  return StateComponentCPP.EditingState(GridID);

}

void ovkEditState(ovk_state_component *StateComponent, int GridID, ovk_state **State) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");
  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");

  auto &StateComponentCPP = *reinterpret_cast<ovk::state_component *>(StateComponent);

  ovk::edit_handle<ovk::state> EditHandle = StateComponentCPP.EditState(GridID);
  auto StateCPPPtr = EditHandle.Release();

  *State = reinterpret_cast<ovk_state *>(StateCPPPtr);

}

void ovkRestoreState(ovk_state_component *StateComponent, int GridID, ovk_state **State) {

  OVK_DEBUG_ASSERT(StateComponent, "Invalid state component pointer.");
  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(*State, "Invalid state pointer.");

  auto &StateComponentCPP = *reinterpret_cast<ovk::state_component *>(StateComponent);
  StateComponentCPP.RestoreState(GridID);

  *State = nullptr;

}

void ovkCreateStateComponentParams(ovk_state_component_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::state_component::params();

  *Params = reinterpret_cast<ovk_state_component_params *>(ParamsCPPPtr);

}

void ovkDestroyStateComponentParams(ovk_state_component_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::state_component::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetStateComponentParamName(const ovk_state_component_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::state_component::params *>(Params);
  std::strcpy(Name, ParamsCPP.Name().c_str());

}

void ovkSetStateComponentParamName(ovk_state_component_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::state_component::params *>(Params);
  ParamsCPP.SetName(Name);

}

}
