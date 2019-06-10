// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Domain.hpp"

#include "ovk/core/Comm.hpp"
#include "ovk/core/Component.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/DataType.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/IDMap.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <memory>
#include <string>
#include <utility>

namespace ovk {

domain::domain(std::shared_ptr<context> &&Context, params &&Params):
  domain_base(std::move(Context), std::move(*Params.Name_), Params.NumDims_, Params.Comm_)
{

  const core::logger &Logger = Context_->core_Logger();

  MPI_Barrier(Comm_);

  if (Comm_.Rank() == 0 && Logger.LoggingDebug()) {
    std::string ProcessesString = core::FormatNumber(Comm_.Size(), "processes", "process");
    Logger.LogDebug(true, 0, "Created %1iD domain %s on %s.", NumDims_, *Name_, ProcessesString);
  }

  Logger.LogStatus(Comm_.Rank() == 0, 0, "Done creating domain %s.", *Name_);

}

domain::~domain() noexcept {

  if (Context_) {
    MPI_Barrier(Comm_);
    core::logger &Logger = Context_->core_Logger();
    Logger.LogStatus(Comm_.Rank() == 0, 0, "Destroying domain %s...", *Name_);
  }

}

domain domain::internal_Create(std::shared_ptr<context> &&Context, params &&Params) {

  return {std::move(Context), std::move(Params)};

}

domain CreateDomain(std::shared_ptr<context> Context, domain::params Params) {

  OVK_DEBUG_ASSERT(Context, "Invalid context.");

  return domain::internal_Create(std::move(Context), std::move(Params));

}

bool domain::ComponentExists(int ComponentID) const {

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");

  return Components_.Contains(ComponentID);

}

void domain::DestroyComponent(int ComponentID) {

  MPI_Barrier(Comm_);

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");
  OVK_DEBUG_ASSERT(ComponentExists(ComponentID), "Component %i does not exist.", ComponentID);

  ComponentEvent_.Trigger(ComponentID, component_event_flags::DESTROY);

  MPI_Barrier(Comm_);

  Components_.Erase(ComponentID);

  MPI_Barrier(Comm_);

}

bool domain::EditingComponent(int ComponentID) const {

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");
  OVK_DEBUG_ASSERT(ComponentExists(ComponentID), "Component %i does not exist.", ComponentID);

  const component_data &Data = Components_(ComponentID);
  const editor &Editor = Data.Editor;

  return Editor.Active();

}

void domain::RestoreComponent(int ComponentID) {

  OVK_DEBUG_ASSERT(ComponentID >= 0, "Invalid component ID.");
  OVK_DEBUG_ASSERT(ComponentExists(ComponentID), "Component %i does not exist.", ComponentID);

  component_data &Data = Components_(ComponentID);
  editor &Editor = Data.Editor;

  OVK_DEBUG_ASSERT(Editor.Active(), "Unable to restore component %i; not currently being edited.",
    ComponentID);

  Editor.Restore();

}

domain::params &domain::params::SetName(std::string Name) {

  Name_ = std::move(Name);

  return *this;

}

domain::params &domain::params::SetDimension(int NumDims) {

  NumDims_ = NumDims;

  return *this;

}

domain::params &domain::params::SetComm(MPI_Comm Comm) {

  Comm_ = Comm;

  return *this;

}

}
