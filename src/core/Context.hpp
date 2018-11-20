// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONTEXT_HPP_INCLUDED
#define OVK_CORE_CONTEXT_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Domain.hpp>
#include <ovk/core/ErrorHandler.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Logger.hpp>

#include <mpi.h>

#include <list>

namespace ovk {

struct context_params {
  MPI_Comm Comm_;
  log_level LogLevel_;
  error_handler_type ErrorHandlerType_;
};

struct context {
  mutable core::logger Logger_;
  mutable core::error_handler ErrorHandler_;
  MPI_Comm Comm_;
  int CommSize_;
  int CommRank_;
  std::list<domain> Domains_;
};

void CreateContextParams(context_params &Params);
void DestroyContextParams(context_params &Params);
void GetContextParamComm(const context_params &Params, MPI_Comm &Comm);
void SetContextParamComm(context_params &Params, MPI_Comm Comm);
void GetContextParamLogLevel(const context_params &Params, log_level &LogLevel);
void SetContextParamLogLevel(context_params &Params, log_level LogLevel);
void GetContextParamErrorHandlerType(const context_params &Params, error_handler_type &
  ErrorHandlerType);
void SetContextParamErrorHandlerType(context_params &Params, error_handler_type ErrorHandlerType);

error CreateContext(context &Context, const context_params &Params);
void DestroyContext(context &Context);

void GetContextComm(const context &Context, MPI_Comm &Comm);
void GetContextLogLevel(const context &Context, log_level &LogLevel);
void SetContextLogLevel(context &Context, log_level LogLevel);
void GetContextErrorHandlerType(const context &Context, error_handler_type &ErrorHandlerType);
void SetContextErrorHandlerType(context &Context, error_handler_type ErrorHandlerType);

void CreateDomain(context &Context, domain *&Domain, const domain_params &Params);
void DestroyDomain(context &Context, domain *&Domain);

}

#endif
