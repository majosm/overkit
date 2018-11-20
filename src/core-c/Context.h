// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_CONTEXT_H_INCLUDED
#define OVK_CORE_C_CONTEXT_H_INCLUDED

#include <ovk/core-c/Constants.h>
#include <ovk/core-c/Domain.h>
#include <ovk/core-c/Global.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_context_params;
typedef struct ovk_context_params ovk_context_params;

struct ovk_context;
typedef struct ovk_context ovk_context;

void ovkCreateContextParams(ovk_context_params **Params);
void ovkDestroyContextParams(ovk_context_params **Params);
void ovkGetContextParamComm(const ovk_context_params *Params, MPI_Comm *Comm);
void ovkSetContextParamComm(ovk_context_params *Params, MPI_Comm Comm);
void ovkGetContextParamLogLevel(const ovk_context_params *Params, ovk_log_level *LogLevel);
void ovkSetContextParamLogLevel(ovk_context_params *Params, ovk_log_level LogLevel);
void ovkGetContextParamErrorHandlerType(const ovk_context_params *Params, ovk_error_handler_type *
  ErrorHandlerType);
void ovkSetContextParamErrorHandlerType(ovk_context_params *Params, ovk_error_handler_type
  ErrorHandlerType);

ovk_error ovkCreateContext(ovk_context **Context, const ovk_context_params *Params);
void ovkDestroyContext(ovk_context **Context);

void ovkGetContextComm(const ovk_context *Context, MPI_Comm *Comm);
void ovkGetContextLogLevel(const ovk_context *Context, ovk_log_level *LogLevel);
void ovkSetContextLogLevel(ovk_context *Context, ovk_log_level LogLevel);
void ovkGetContextErrorHandlerType(const ovk_context *Context, ovk_error_handler_type
  *ErrorHandlerType);
void ovkSetContextErrorHandlerType(ovk_context *Context, ovk_error_handler_type ErrorHandlerType);

void ovkCreateDomain(ovk_context *Context, ovk_domain **Domain, const ovk_domain_params *Params);
void ovkDestroyDomain(ovk_context *Context, ovk_domain **Domain);

#ifdef __cplusplus
}
#endif

#endif
