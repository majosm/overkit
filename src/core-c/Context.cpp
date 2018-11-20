// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Context.h"

#include "ovk/core-c/Constants.h"
#include "ovk/core-c/Domain.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Context.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

ovk_error ovkCreateContext(ovk_context **Context, const ovk_context_params *Params) {

  int MPIInitialized;
  MPI_Initialized(&MPIInitialized);
  if (!MPIInitialized) {
    fprintf(stderr, "ERROR: MPI not initialized.\n"); fflush(stderr);
    bool Exit;
    if (Params) {
      auto &ParamsCPP = *reinterpret_cast<const ovk::context_params *>(Params);
      Exit = ParamsCPP.ErrorHandlerType_ == ovk::error_handler_type::ABORT;
    } else {
      Exit = OVK_DEBUG;
    }
    if (Exit) {
      exit(OVK_ERROR_MPI);
    } else {
      return OVK_ERROR_MPI;
    }
  }

  // Have to do this after checking for MPI because it potentially calls MPI_Abort
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ContextCPPPtr = new ovk::context();
  auto &ParamsCPP = *reinterpret_cast<const ovk::context_params *>(Params);

  ovk::error ErrorCPP = ovk::CreateContext(*ContextCPPPtr, ParamsCPP);

  *Context = reinterpret_cast<ovk_context *>(ContextCPPPtr);

  return ovk_error(ErrorCPP);

}

void ovkDestroyContext(ovk_context **Context) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(*Context, "Invalid context pointer.");

  auto ContextCPPPtr = reinterpret_cast<ovk::context *>(*Context);

  ovk::DestroyContext(*ContextCPPPtr);

  delete ContextCPPPtr;

  *Context = nullptr;

}

void ovkGetContextComm(const ovk_context *Context, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);
  ovk::GetContextComm(ContextCPP, *Comm);

}

void ovkGetContextLogLevel(const ovk_context *Context, ovk_log_level *LogLevel) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(LogLevel, "Invalid log level pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);

  ovk::log_level LogLevelCPP;
  ovk::GetContextLogLevel(ContextCPP, LogLevelCPP);

  *LogLevel = ovk_log_level(LogLevelCPP);

}

void ovkSetContextLogLevel(ovk_context *Context, ovk_log_level LogLevel) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(ovkValidLogLevel(LogLevel), "Invalid log level.");

  auto &ContextCPP = *reinterpret_cast<ovk::context *>(Context);
  ovk::SetContextLogLevel(ContextCPP, ovk::log_level(LogLevel));

}

void ovkGetContextErrorHandlerType(const ovk_context *Context, ovk_error_handler_type
  *ErrorHandlerType) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(ErrorHandlerType, "Invalid error handler type pointer.");

  auto &ContextCPP = *reinterpret_cast<const ovk::context *>(Context);

  ovk::error_handler_type ErrorHandlerTypeCPP;
  ovk::GetContextErrorHandlerType(ContextCPP, ErrorHandlerTypeCPP);

  *ErrorHandlerType = ovk_error_handler_type(ErrorHandlerTypeCPP);

}

void ovkSetContextErrorHandlerType(ovk_context *Context, ovk_error_handler_type ErrorHandlerType) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(ovkValidErrorHandlerType(ErrorHandlerType), "Invalid error handler type.");

  auto &ContextCPP = *reinterpret_cast<ovk::context *>(Context);
  ovk::SetContextErrorHandlerType(ContextCPP, ovk::error_handler_type(ErrorHandlerType));

}

void ovkCreateDomain(ovk_context *Context, ovk_domain **Domain, const ovk_domain_params *Params) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ContextCPP = *reinterpret_cast<ovk::context *>(Context);
  auto &ParamsCPP = *reinterpret_cast<const ovk::domain_params *>(Params);
  ovk::domain *DomainCPPPtr;

  ovk::CreateDomain(ContextCPP, DomainCPPPtr, ParamsCPP);

  *Domain = reinterpret_cast<ovk_domain *>(DomainCPPPtr);

}

void ovkDestroyDomain(ovk_context *Context, ovk_domain **Domain) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(*Domain, "Invalid domain pointer.");

  auto &ContextCPP = *reinterpret_cast<ovk::context *>(Context);
  auto DomainCPPPtr = reinterpret_cast<ovk::domain *>(*Domain);

  ovk::DestroyDomain(ContextCPP, DomainCPPPtr);

  *Domain = nullptr;

}

void ovkCreateContextParams(ovk_context_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::context_params();

  ovk::CreateContextParams(*ParamsCPPPtr);

  *Params = reinterpret_cast<ovk_context_params *>(ParamsCPPPtr);

}

void ovkDestroyContextParams(ovk_context_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::context_params *>(*Params);

  ovk::DestroyContextParams(*ParamsCPPPtr);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetContextParamComm(const ovk_context_params *Params, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::context_params *>(Params);
  ovk::GetContextParamComm(ParamsCPP, *Comm);

}

void ovkSetContextParamComm(ovk_context_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::context_params *>(Params);
  ovk::SetContextParamComm(ParamsCPP, Comm);

}

void ovkGetContextParamLogLevel(const ovk_context_params *Params, ovk_log_level *LogLevel) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LogLevel, "Invalid log level pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::context_params *>(Params);

  ovk::log_level LogLevelCPP;
  ovk::GetContextParamLogLevel(ParamsCPP, LogLevelCPP);

  *LogLevel = ovk_log_level(LogLevelCPP);

}

void ovkSetContextParamLogLevel(ovk_context_params *Params, ovk_log_level LogLevel) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::context_params *>(Params);
  ovk::SetContextParamLogLevel(ParamsCPP, ovk::log_level(LogLevel));

}

void ovkGetContextParamErrorHandlerType(const ovk_context_params *Params, ovk_error_handler_type
  *ErrorHandlerType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(ErrorHandlerType, "Invalid error handler type pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::context_params *>(Params);

  ovk::error_handler_type ErrorHandlerTypeCPP;
  ovk::GetContextParamErrorHandlerType(ParamsCPP, ErrorHandlerTypeCPP);

  *ErrorHandlerType = ovk_error_handler_type(ErrorHandlerTypeCPP);

}

void ovkSetContextParamErrorHandlerType(ovk_context_params *Params, ovk_error_handler_type
  ErrorHandlerType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::context_params *>(Params);
  ovk::SetContextParamErrorHandlerType(ParamsCPP, ovk::error_handler_type(ErrorHandlerType));

}
