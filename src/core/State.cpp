// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/State.hpp"

#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Editor.hpp"
#include "ovk/core/Event.hpp"
#include "ovk/core/Field.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Partition.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

namespace state_internal {

state_base::state_base(std::shared_ptr<context> &&Context, const grid &Grid):
  Context_(std::move(Context)),
  Grid_(&Grid),
  Comm_(Grid.Comm())
{
  MPI_Barrier(Comm_);
}

state_base::~state_base() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Logger.LogDebug(Comm_.Rank() == 0, 0, "Destroyed state %s.", Grid_->Name());
  }

}

}

state::state(std::shared_ptr<context> &&Context, const grid &Grid, params &&Params):
  state_base(std::move(Context), Grid),
  Flags_(Grid.SharedPartition(), state_flags::ACTIVE)
{

  MPI_Barrier(Comm_);

  core::logger &Logger = Context_->core_Logger();
  Logger.LogDebug(Comm_.Rank() == 0, 0, "Created state %s.", Grid.Name());

}

state::~state() noexcept {

  if (Context_) {
    // Barrier before cleaning up
    MPI_Barrier(Comm_);
  }

}

state state::internal_Create(std::shared_ptr<context> &&Context, const grid &Grid, params
  &&Params) {

  return {std::move(Context), Grid, std::move(Params)};

}

namespace core {

state CreateState(std::shared_ptr<context> Context, const grid &Grid, state::params Params)
  {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return state::internal_Create(std::move(Context), Grid, std::move(Params));

}

}

bool state::EditingFlags() const {

  return FlagsEditor_.Active();

}

edit_handle<distributed_field<state_flags>> state::EditFlags() {

  if (!FlagsEditor_.Active()) {
    MPI_Barrier(Comm_);
    floating_ref<state> FloatingRef = FloatingRefGenerator_.Generate(*this);
    auto DeactivateFunc = [FloatingRef] {
      state &State = *FloatingRef;
      State.OnFlagsEndEdit_();
      MPI_Barrier(State.Comm_);
    };
    FlagsEditor_.Activate(std::move(DeactivateFunc));
  }

  return FlagsEditor_.Edit(Flags_);

}

void state::RestoreFlags() {

  OVK_DEBUG_ASSERT(FlagsEditor_.Active(), "Unable to restore flags; not currently being edited.");

  FlagsEditor_.Restore();

}

void state::OnFlagsEndEdit_() {

  const grid &Grid = *Grid_;

  Grid.Partition().Exchange(Flags_);

  MPI_Barrier(Comm_);

  FlagsEvent_.Trigger();

}

}
