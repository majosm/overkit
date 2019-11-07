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

context_base::context_base(MPI_Comm Comm, bool LoggingErrors, bool LoggingWarnings, int
  StatusLoggingThreshold):
  Exists_(true),
  Comm_(DuplicateComm(Comm)),
  Logger_(Comm_.Rank())
{

  if (LoggingErrors) Logger_.EnableErrors();
  if (LoggingWarnings) Logger_.EnableWarnings();
  Logger_.SetStatusThreshold(StatusLoggingThreshold);

  Logger_.LogStatus(Comm_.Rank() == 0, "Creating context...");
  Level1_ = Logger_.IncreaseStatusLevelAndIndent();

}

context_base::~context_base() noexcept {

  if (Exists_) {
    MPI_Barrier(Comm_);
    Level1_.Reset();
    Logger_.LogStatus(Comm_.Rank() == 0, "Done destroying context.");
  }

}

}

context::context(params &&Params):
  context_base(Params.Comm_, Params.ErrorLogging_, Params.WarningLogging_,
    Params.StatusLoggingThreshold_),
  Profiler_(Comm_)
{

  MPI_Comm_set_errhandler(Comm_, MPI_ERRORS_RETURN);

  if (Params.Profiling_) {
    Profiler_.Enable();
  }

  MPI_Barrier(Comm_);

  if (Comm_.Rank() == 0 && Logger_.LoggingStatus()) {
    std::string ProcessesString = core::FormatNumber(Comm_.Size(), "processes", "process");
    Logger_.LogStatus(true, "Created context on %s.", ProcessesString);
  }

  Level1_.Reset();
  Logger_.LogStatus(Comm_.Rank() == 0, "Done creating context.");

}

context::~context() noexcept {

  if (Exists_) {

    // Barrier before cleaning up
    MPI_Barrier(Comm_);

    Logger_.LogStatus(Comm_.Rank() == 0, "Destroying context...");
    Level1_ = Logger_.IncreaseStatusLevelAndIndent();

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

void context::EnableErrorLogging() {

  Logger_.EnableErrors();

}

void context::DisableErrorLogging() {

  Logger_.DisableErrors();

}

void context::EnableWarningLogging() {

  Logger_.EnableWarnings();

}

void context::DisableWarningLogging() {

  Logger_.DisableWarnings();

}

void context::SetStatusLoggingThreshold(int StatusLoggingThreshold) {

  OVK_DEBUG_ASSERT(StatusLoggingThreshold >= 0, "Invalid status logging threshold.");

  Logger_.SetStatusThreshold(StatusLoggingThreshold);

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

context::params &context::params::SetErrorLogging(bool ErrorLogging) {

  ErrorLogging_ = ErrorLogging;

  return *this;

}

context::params &context::params::SetWarningLogging(bool WarningLogging) {

  WarningLogging_ = WarningLogging;

  return *this;

}

context::params &context::params::SetStatusLoggingThreshold(int StatusLoggingThreshold) {

  OVK_DEBUG_ASSERT(StatusLoggingThreshold >= 0, "Invalid status logging threshold.");

  StatusLoggingThreshold_ = StatusLoggingThreshold;

  return *this;

}

context::params &context::params::SetProfiling(bool Profiling) {

  Profiling_ = Profiling;

  return *this;

}

}
