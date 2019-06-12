// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Context.hpp"

#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Error.hpp"
#include "ovk/core/FloatingRef.hpp"
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
  Comm_(Comm),
  Logger_(LogLevel, Comm_.Rank())
{
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
  FloatingRefGenerator_(*this),
  Profiler_(Comm_)
{

  MPI_Comm_set_errhandler(Comm_, MPI_ERRORS_RETURN);

  if (OVK_TIMERS) Profiler_.Enable();

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

    // May want to move this somewhere else
    std::string ProfileTimesString = Profiler_.WriteProfile();
    if (Comm_.Rank() == 0) {
      printf("%s", ProfileTimesString.c_str());
    }

    Logger_.LogStatus(Comm_.Rank() == 0, 0, "Destroying context...");

  }

}

context context::internal_Create(params &&Params) {

  if (OVK_DEBUG) {
    int MPIInitialized;
    MPI_Initialized(&MPIInitialized);
    // Can't use OVK_DEBUG_ASSERT here because it calls MPI_Abort
    if (!MPIInitialized) {
      fprintf(stderr, "ERROR: MPI not initialized.\n"); fflush(stderr);
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

void context::SetLogLevel(log_level LogLevel) {

  Logger_.SetLevel(LogLevel);

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

}
