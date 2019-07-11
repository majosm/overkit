// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_STATE_H_INCLUDED
#define OVK_CORE_C_STATE_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/Global.h>
#include <ovk/core-c/Grid.h>
#include <ovk/core-c/State.h>
#include <ovk/core/State.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_state;
typedef struct ovk_state ovk_state;

struct ovk_state_params;
typedef struct ovk_state_params ovk_state_params;

void ovkGetStateContextC(const ovk_state *State, const ovk_context **Context);
void ovkGetStateContext(ovk_state *State, ovk_context **Context);
void ovkGetStateSharedContext(ovk_state *State, ovk_shared_context **Context);

void ovkGetStateGrid(const ovk_state *State, const ovk_grid **Grid);

void ovkGetStateComm(const ovk_state *State, MPI_Comm *Comm);
void ovkGetStateCommSize(const ovk_state *State, int *CommSize);
void ovkGetStateCommRank(const ovk_state *State, int *CommRank);

void ovkGetStateFlags(const ovk_state *State, const ovk_state_flags **Flags);
bool ovkEditingStateFlags(const ovk_state *State);
void ovkEditStateFlags(ovk_state *State, ovk_state_flags **Flags);
void ovkRestoreStateFlags(ovk_state *State, ovk_state_flags **Flags);

void ovkCreateStateParams(ovk_state_params **Params);
void ovkDestroyStateParams(ovk_state_params **Params);

#ifdef __cplusplus
}
#endif

#endif
