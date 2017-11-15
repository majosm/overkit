// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CONTEXT_INCLUDED
#define OVK_CORE_PUBLIC_CONTEXT_INCLUDED

#include <ovkGlobal.h>
#include <ovkDomain.h>

struct ovk_context_params;
typedef struct ovk_context_params ovk_context_params;

struct ovk_context;
typedef struct ovk_context ovk_context;

struct ovk_context_properties;
typedef struct ovk_context_properties ovk_context_properties;

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

void ovkGetContextProperties(const ovk_context *Context, const ovk_context_properties **Properties);
// void ovkEditContextProperties(ovk_context *Context, ovk_context_properties **Properties);
// void ovkReleaseContextProperties(ovk_context *Context, ovk_context_properties **Properties);

void ovkCreateDomainParams(ovk_context *Context, ovk_domain_params **Params);
void ovkDestroyDomainParams(ovk_context *Context, ovk_domain_params **Params);

void ovkCreateDomain(ovk_context *Context, ovk_domain **Domain, const ovk_domain_params *Params);
void ovkDestroyDomain(ovk_context *Context, ovk_domain **Domain);

void ovkGetContextPropertyComm(const ovk_context_properties *Properties, MPI_Comm *Comm);
void ovkGetContextPropertyLogLevel(const ovk_context_properties *Properties, ovk_log_level *LogLevel);
void ovkSetContextPropertyLogLevel(ovk_context_properties *Properties, ovk_log_level LogLevel);
void ovkGetContextPropertyErrorHandlerType(const ovk_context_properties *Properties,
  ovk_error_handler_type *ErrorHandlerType);
void ovkSetContextPropertyErrorHandlerType(ovk_context_properties *Properties,
  ovk_error_handler_type ErrorHandlerType);

#endif
