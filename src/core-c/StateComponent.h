// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_STATE_COMPONENT_H_INCLUDED
#define OVK_CORE_C_STATE_COMPONENT_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core-c/State.h>
#include <ovk/core/StateComponent.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_state_component;
typedef struct ovk_state_component ovk_state_component;

struct ovk_state_component_params;
typedef struct ovk_state_component_params ovk_state_component_params;

int ovkStateCount(const ovk_state_component *StateComponent);

bool ovkStateExists(const ovk_state_component *StateComponent, int GridID);

void ovkCreateState(ovk_state_component *StateComponent, int GridID, ovk_state_params
  **MaybeStateParams);
void ovkCreateStates(ovk_state_component *StateComponent, int Count, const int *GridIDs,
  ovk_state_params **MaybeStateParams);

void ovkDestroyState(ovk_state_component *StateComponent, int GridID);
void ovkDestroyStates(ovk_state_component *StateComponent, int Count, const int *GridIDs);

int LocalStateCount(const ovk_state_component *StateComponent);

void ovkGetState(const ovk_state_component *StateComponent, int GridID, const ovk_state **State);
bool ovkEditingState(const ovk_state_component *StateComponent, int GridID);
void ovkEditState(ovk_state_component *StateComponent, int GridID, ovk_state **State);
void ovkRestoreState(ovk_state_component *StateComponent, int GridID, ovk_state **State);

void ovkCreateStateComponentParams(ovk_state_component_params **Params);
void ovkDestroyStateComponentParams(ovk_state_component_params **Params);
void ovkGetStateComponentParamName(const ovk_state_component_params *Params, char *Name);
void ovkSetStateComponentParamName(ovk_state_component_params *Params, const char *Name);

#ifdef __cplusplus
}
#endif

#endif
