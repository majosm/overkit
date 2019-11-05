// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Context.hpp"

#include "ovk/core/Comm.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Error.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/TextProcessing.hpp"
#include "ovk/core/Profiler.hpp"

#include <mpi.h>

#include <string>
#include <utility>

namespace ovk {

namespace context_internal {

context_base::context_base(MPI_Comm Comm, log_level LogLevel):
  Exists_(true),
  Comm_(DuplicateComm(Comm)),
  Logger_(Comm_.Rank())
{

  if ((LogLevel & log_level::ERRORS) != log_level::NONE) Logger_.EnableErrorLogging();
  if ((LogLevel & log_level::WARNINGS) != log_level::NONE) Logger_.EnableWarningLogging();
  if ((LogLevel & log_level::STATUS) != log_level::NONE) Logger_.EnableStatusLogging();
  if ((LogLevel & log_level::DEBUG) != log_level::NONE) Logger_.EnableDebugLogging();

  Logger_.LogStatus(Comm_.Rank() == 0, 0, "Creating context...");

}

context_base::~context_base() noexcept {

  if (Exists_) {
    MPI_Barrier(Comm_);
    Logger_.LogStatus(Comm_.Rank() == 0, 0, "Done destroying context.");
  }

}

}

context::context(params &&Params):
  context_base(Params.Comm_, Params.LogLevel_),
  Profiler_(Comm_)
{

  MPI_Comm_set_errhandler(Comm_, MPI_ERRORS_RETURN);

  if (Params.Profiling_) {
    Profiler_.Enable();
  }

  MPI_Barrier(Comm_);

  if (Comm_.Rank() == 0 && Logger_.LoggingDebug()) {
    std::string ProcessesString = core::FormatNumber(Comm_.Size(), "processes", "process");
    Logger_.LogDebug(true, 0, "Created context on %s.", ProcessesString);
  }

  Logger_.LogStatus(Comm_.Rank() == 0, 0, "Done creating context.");

}

context::~context() noexcept {

  if (Exists_) {

    // Barrier before cleaning up
    MPI_Barrier(Comm_);

    Logger_.LogStatus(Comm_.Rank() == 0, 0, "Destroying context...");

  }

}

context context::internal_Create(params &&Params) {

  if (OVK_DEBUG) {
    int MPIInitialized;
    MPI_Initialized(&MPIInitialized);
    // Can't use OVK_DEBUG_ASSERT here because it calls MPI_Abort
    if (!MPIInitialized) {
      std::fprintf(stderr, "ERROR: MPI not initialized.\n");
      std::fflush(stderr);
      exit(1);
    }
  }

  return {std::move(Params)};

}

context CreateContext(context::params Params) {

  return context::internal_Create(std::move(Params));

}

optional<context> CreateContext(context::params Params, error &Error) {

  Error = error::NONE;

  optional<context> MaybeContext;

  try {
    MaybeContext = CreateContext(std::move(Params));
  } catch (const exception &Exception) {
    Error = Exception.Error();
  }

  return MaybeContext;

}

log_level context::LogLevel() const {

  log_level LogLevel = log_level::NONE;

  if (Logger_.LoggingErrors()) LogLevel &= log_level::ERRORS;
  if (Logger_.LoggingWarnings()) LogLevel &= log_level::WARNINGS;
  if (Logger_.LoggingStatus()) LogLevel &= log_level::STATUS;
  if (Logger_.LoggingDebug()) LogLevel &= log_level::DEBUG;

  return LogLevel;

}

void context::SetLogLevel(log_level LogLevel) {

  if ((LogLevel & log_level::ERRORS) != log_level::NONE) {
    Logger_.EnableErrorLogging();
  } else {
    Logger_.DisableErrorLogging();
  }
  if ((LogLevel & log_level::WARNINGS) != log_level::NONE) {
    Logger_.EnableWarningLogging();
  } else {
    Logger_.DisableWarningLogging();
  }
  if ((LogLevel & log_level::STATUS) != log_level::NONE) {
    Logger_.EnableStatusLogging();
  } else {
    Logger_.DisableStatusLogging();
  }
  if ((LogLevel & log_level::DEBUG) != log_level::NONE) {
    Logger_.EnableDebugLogging();
  } else {
    Logger_.DisableDebugLogging();
  }

}

void context::EnableProfiling() {

  Profiler_.Enable();

}

void context::DisableProfiling() {

  Profiler_.Disable();

}

std::string context::WriteProfile() const {

  return Profiler_.WriteProfile();

}

context::params &context::params::SetComm(MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Comm_ = Comm;

  return *this;

}

context::params &context::params::SetLogLevel(log_level LogLevel) {

  OVK_DEBUG_ASSERT(ValidLogLevel(LogLevel), "Invalid log level.");

  LogLevel_ = LogLevel;

  return *this;

}

context::params &context::params::SetProfiling(bool Profiling) {

  Profiling_ = Profiling;

  return *this;

}

}
