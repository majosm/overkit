// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/State.h"

#include "ovk/core-c/Context.h"
#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/DistributedField.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/State.hpp"

#include <mpi.h>

#include <memory>

void ovkGetStateContextC(const ovk_state *State, const ovk_context **Context) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &StateCPP = *reinterpret_cast<const ovk::state *>(State);
  *Context = reinterpret_cast<const ovk_context *>(&StateCPP.Context());

}

void ovkGetStateContext(ovk_state *State, ovk_context **Context) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &StateCPP = *reinterpret_cast<ovk::state *>(State);
  *Context = reinterpret_cast<ovk_context *>(&StateCPP.Context());

}

void ovkGetStateSharedContext(ovk_state *State, ovk_shared_context **Context) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &StateCPP = *reinterpret_cast<ovk::state *>(State);
  auto &ContextCPP = StateCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

void ovkGetStateGrid(const ovk_state *State, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &StateCPP = *reinterpret_cast<const ovk::state *>(State);
  *Grid = reinterpret_cast<const ovk_grid *>(&StateCPP.Grid());

}

void ovkGetStateComm(const ovk_state *State, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &StateCPP = *reinterpret_cast<const ovk::state *>(State);
  *Comm = StateCPP.Comm();

}

void ovkGetStateCommSize(const ovk_state *State, int *CommSize) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &StateCPP = *reinterpret_cast<const ovk::state *>(State);
  *CommSize = StateCPP.Comm().Size();

}

void ovkGetStateCommRank(const ovk_state *State, int *CommRank) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &StateCPP = *reinterpret_cast<const ovk::state *>(State);
  *CommRank = StateCPP.Comm().Rank();

}

void ovkGetStateFlags(const ovk_state *State, const ovk_state_flags **Flags) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(Flags, "Invalid flags pointer.");

  auto &StateCPP = *reinterpret_cast<const ovk::state *>(State);
  *Flags = reinterpret_cast<const ovk_state_flags *>(StateCPP.Flags().Data());

}

bool ovkEditingStateFlags(const ovk_state *State) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");

  auto &StateCPP = *reinterpret_cast<const ovk::state *>(State);
  return StateCPP.EditingFlags();

}

void ovkEditStateFlags(ovk_state *State, ovk_state_flags **Flags) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(Flags, "Invalid flags pointer.");

  auto &StateCPP = *reinterpret_cast<ovk::state *>(State);

  ovk::edit_handle<ovk::distributed_field<ovk::state_flags>> EditHandle = StateCPP.EditFlags();
  auto &FlagsCPP = *EditHandle.Release();

  *Flags = reinterpret_cast<ovk_state_flags *>(FlagsCPP.Data());

}

void ovkRestoreStateFlags(ovk_state *State, ovk_state_flags **Flags) {

  OVK_DEBUG_ASSERT(State, "Invalid state pointer.");
  OVK_DEBUG_ASSERT(Flags, "Invalid flags pointer.");

  auto &StateCPP = *reinterpret_cast<ovk::state *>(State);
  StateCPP.RestoreFlags();

  *Flags = nullptr;

}

void ovkCreateStateParams(ovk_state_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::state::params();

  *Params = reinterpret_cast<ovk_state_params *>(ParamsCPPPtr);

}

void ovkDestroyStateParams(ovk_state_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::state::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}
