// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "Context.h"

#include "Debug.h"
#include "Domain.h"
#include "ErrorHandler.h"
#include "Global.h"
#include "List.h"
#include "Logger.h"

static void DefaultContextProperties(ovk_context_properties *Properties);

ovk_error ovkCreateContext(ovk_context **Context_, const ovk_context_params *Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid context params pointer.");

  MPI_Comm Comm;
  MPI_Comm_dup(Params->comm, &Comm);

  MPI_Barrier(Comm);

  *Context_ = malloc(sizeof(ovk_context));
  ovk_context *Context = *Context_;

  DefaultContextProperties(&Context->properties);

  Context->properties.comm = Comm;
  MPI_Comm_size(Comm, &Context->properties.comm_size);
  MPI_Comm_rank(Comm, &Context->properties.comm_rank);

  Context->properties.log_level = Params->log_level;
  CreateLogger(&Context->logger, Params->log_level);

  Context->properties.error_handler_type = Params->error_handler_type;
  CreateErrorHandler(&Context->error_handler, Params->error_handler_type);

  ListCreate(&Context->domains);

  MPI_Barrier(Context->properties.comm);

  LogStatus(Context->logger, Context->properties.comm_rank == 0, 0, "Context created.");

  return OVK_ERROR_NONE;

}

void ovkDestroyContext(ovk_context **Context_) {

  ovk_context *Context = *Context_;

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  MPI_Barrier(Context->properties.comm);

  t_list_entry *Entry = ListBegin(Context->domains);
  while (Entry != ListEnd(Context->domains)) {
    ovk_domain *Domain = ListData(Entry);
    DestroyDomain(&Domain);
    Entry = ListNext(Entry);
  }
  ListDestroy(&Context->domains);

  t_logger *Logger = Context->logger;

  DestroyErrorHandler(&Context->error_handler);

  MPI_Comm Comm = Context->properties.comm;
  bool IsCoreRank = Context->properties.comm_rank == 0;

  free(*Context_);
  *Context_ = NULL;

  MPI_Barrier(Comm);

  LogStatus(Logger, IsCoreRank, 0, "Context destroyed.");
  DestroyLogger(&Logger);

  MPI_Comm_free(&Comm);

}

void ovkGetContextProperties(const ovk_context *Context, const ovk_context_properties **Properties)
  {

  *Properties = &Context->properties;

}

void ovkCreateDomainParams(ovk_context *Context, ovk_domain_params **Params) {

  CreateDomainParams(Params, Context->properties.comm);

}

void ovkDestroyDomainParams(ovk_context *Context, ovk_domain_params **Params) {

  OVK_DEBUG_ASSERT(*Params, "Invalid domain params pointer.");

  DestroyDomainParams(Params);

}

void ovkCreateDomain(ovk_context *Context, ovk_domain **Domain_, const ovk_domain_params *Params) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Domain_, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid domain params pointer.");

  ovk_domain *Domain;
  CreateDomain(&Domain, Params, Context->logger, Context->error_handler);

  ListPushBack(Context->domains, Domain);

  *Domain_ = Domain;

}

void ovkDestroyDomain(ovk_context *Context, ovk_domain **Domain_) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Domain_, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(*Domain_, "Invalid domain pointer.");

  t_list_entry *Entry = ListBegin(Context->domains);
  while (Entry != ListEnd(Context->domains)) {
    ovk_domain *Domain = ListData(Entry);
    if (Domain == *Domain_) {
      break;
    }
    Entry = ListNext(Entry);
  }

  OVK_DEBUG_ASSERT(Entry != ListEnd(Context->domains), "Domain does not belong to this context.");

  ovk_domain *Domain = ListRemove(Context->domains, &Entry);

  DestroyDomain(&Domain);

  *Domain_ = NULL;

}

void ovkCreateContextParams(ovk_context_params **Params_) {

  *Params_ = malloc(sizeof(ovk_context_params));
  ovk_context_params *Params = *Params_;

  Params->comm = MPI_COMM_NULL;
  Params->log_level = OVK_LOG_ALL;
  if (OVK_DEBUG) {
    Params->error_handler_type = OVK_ERROR_HANDLER_ABORT;
  } else {
    Params->error_handler_type = OVK_ERROR_HANDLER_RETURN;
  }

}

void ovkDestroyContextParams(ovk_context_params **Params) {

  OVK_DEBUG_ASSERT(*Params, "Invalid context params pointer.");

  free(*Params);
  *Params = NULL;

}

void ovkGetContextParamComm(const ovk_context_params *Params, MPI_Comm *Comm) {

  *Comm = Params->comm;

}

void ovkSetContextParamComm(ovk_context_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Params->comm = Comm;

}

void ovkGetContextParamLogLevel(const ovk_context_params *Params, ovk_log_level *LogLevel) {

  *LogLevel = Params->log_level;

}

void ovkSetContextParamLogLevel(ovk_context_params *Params, ovk_log_level LogLevel) {

  OVK_DEBUG_ASSERT(ValidLogLevel(Params->log_level), "Invalid log level.");

  Params->log_level = LogLevel;

}

void ovkGetContextParamErrorHandlerType(const ovk_context_params *Params, ovk_error_handler_type *
  ErrorHandlerType) {

  *ErrorHandlerType = Params->error_handler_type;

}

void ovkSetContextParamErrorHandlerType(ovk_context_params *Params, ovk_error_handler_type
  ErrorHandlerType) {

  OVK_DEBUG_ASSERT(ValidErrorHandlerType(Params->error_handler_type), "Invalid error handler type.");

  Params->error_handler_type = ErrorHandlerType;

}

static void DefaultContextProperties(ovk_context_properties *Properties) {

  Properties->comm = MPI_COMM_NULL;
  Properties->comm_size = 0;
  Properties->comm_rank = -1;
  Properties->log_level = OVK_LOG_ALL;
  Properties->error_handler_type = OVK_ERROR_HANDLER_ABORT;

}

void ovkGetContextPropertyComm(const ovk_context_properties *Properties, MPI_Comm *Comm) {

  *Comm = Properties->comm;

}

void ovkGetContextPropertyLogLevel(const ovk_context_properties *Properties,
  ovk_log_level *LogLevel) {

  *LogLevel = Properties->log_level;

}

void ovkSetContextPropertyLogLevel(ovk_context_properties *Properties, ovk_log_level LogLevel) {

  Properties->log_level = LogLevel;

}

void ovkGetContextPropertyErrorHandlerType(const ovk_context_properties *Properties,
  ovk_error_handler_type *ErrorHandlerType) {

  *ErrorHandlerType = Properties->error_handler_type;

}

void ovkSetContextPropertyErrorHandlerType(ovk_context_properties *Properties,
  ovk_error_handler_type ErrorHandlerType) {

  Properties->error_handler_type = ErrorHandlerType;

}
