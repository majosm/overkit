// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Context.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/ErrorHandler.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <list>
#include <string>
#include <utility>

namespace ovk {

error CreateContext(context &Context, const context_params &Params) {

  int MPIInitialized;
  MPI_Initialized(&MPIInitialized);
  if (!MPIInitialized) {
    fprintf(stderr, "ERROR: MPI not initialized.\n"); fflush(stderr);
    if (Params.ErrorHandlerType_ == error_handler_type::ABORT) {
      exit(int(error::MPI));
    } else {
      return error::MPI;
    }
  }

  MPI_Comm_dup(Params.Comm_, &Context.Comm_);

  MPI_Barrier(Context.Comm_);

  MPI_Comm_size(Context.Comm_, &Context.CommSize_);
  MPI_Comm_rank(Context.Comm_, &Context.CommRank_);

  core::CreateLogger(Context.Logger_, Params.LogLevel_, Context.CommRank_);
  core::CreateErrorHandler(Context.ErrorHandler_, Params.ErrorHandlerType_);

  switch (Params.ErrorHandlerType_) {
  case error_handler_type::ABORT:
    MPI_Comm_set_errhandler(Context.Comm_, MPI_ERRORS_ARE_FATAL);
    break;
  case error_handler_type::RETURN:
    MPI_Comm_set_errhandler(Context.Comm_, MPI_ERRORS_RETURN);
    break;
  }

  MPI_Barrier(Context.Comm_);

  if (Context.CommRank_ == 0 && LoggingStatus(Context.Logger_)) {
    std::string ProcessesString = core::FormatNumber(Context.CommSize_, "processes", "process");
    core::LogStatus(Context.Logger_, true, 0, "Created context on %s.", ProcessesString);
  }

  return error::NONE;

}

void DestroyContext(context &Context) {

  MPI_Barrier(Context.Comm_);

  for (auto &Domain : Context.Domains_) {
    core::DestroyDomain(Domain);
  }
  Context.Domains_.clear();

  core::DestroyErrorHandler(Context.ErrorHandler_);

  MPI_Barrier(Context.Comm_);

  core::LogStatus(Context.Logger_, Context.CommRank_ == 0, 0, "Destroyed context.");
  core::DestroyLogger(Context.Logger_);

  MPI_Comm_free(&Context.Comm_);

}

void GetContextComm(const context &Context, MPI_Comm &Comm) {

  Comm = Context.Comm_;

}

void GetContextLogLevel(const context &Context, log_level &LogLevel) {

  GetLogLevel(Context.Logger_, LogLevel);

}

void SetContextLogLevel(context &Context, log_level LogLevel) {

  SetLogLevel(Context.Logger_, LogLevel);

}

void GetContextErrorHandlerType(const context &Context, error_handler_type
  &ErrorHandlerType) {

  GetErrorHandlerType(Context.ErrorHandler_, ErrorHandlerType);

}

void SetContextErrorHandlerType(context &Context, error_handler_type ErrorHandlerType) {

  SetErrorHandlerType(Context.ErrorHandler_, ErrorHandlerType);

}

void CreateDomain(context &Context, domain *&DomainPtr, const domain_params &Params) {

  domain Domain;
  core::CreateDomain(Domain, Params, Context.Logger_, Context.ErrorHandler_);

  Context.Domains_.emplace_back(std::move(Domain));

  DomainPtr = &Context.Domains_.back();

}

void DestroyDomain(context &Context, domain *&DomainPtr) {

  OVK_DEBUG_ASSERT(DomainPtr, "Invalid domain pointer.");

  auto Iter = Context.Domains_.begin();

  while (Iter != Context.Domains_.end()) {
    if (&(*Iter) == DomainPtr) break;
    ++Iter;
  }

  OVK_DEBUG_ASSERT(Iter != Context.Domains_.end(), "Domain does not belong to this context.");

  core::DestroyDomain(*Iter);

  Context.Domains_.erase(Iter);

  DomainPtr = NULL;

}

void CreateContextParams(context_params &Params) {

  Params.Comm_ = MPI_COMM_NULL;
  Params.LogLevel_ = log_level::ALL;
  if (OVK_DEBUG) {
    Params.ErrorHandlerType_ = error_handler_type::ABORT;
  } else {
    Params.ErrorHandlerType_ = error_handler_type::RETURN;
  }

}

void DestroyContextParams(context_params &) {}

void GetContextParamComm(const context_params &Params, MPI_Comm &Comm) {

  Comm = Params.Comm_;

}

void SetContextParamComm(context_params &Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Params.Comm_ = Comm;

}

void GetContextParamLogLevel(const context_params &Params, log_level &LogLevel) {

  LogLevel = Params.LogLevel_;

}

void SetContextParamLogLevel(context_params &Params, log_level LogLevel) {

  OVK_DEBUG_ASSERT(ValidLogLevel(LogLevel), "Invalid log level.");

  Params.LogLevel_ = LogLevel;

}

void GetContextParamErrorHandlerType(const context_params &Params, error_handler_type
  &ErrorHandlerType) {

  ErrorHandlerType = Params.ErrorHandlerType_;

}

void SetContextParamErrorHandlerType(context_params &Params, error_handler_type
  ErrorHandlerType) {

  OVK_DEBUG_ASSERT(ValidErrorHandlerType(ErrorHandlerType), "Invalid error handler type.");

  Params.ErrorHandlerType_ = ErrorHandlerType;

}

}
