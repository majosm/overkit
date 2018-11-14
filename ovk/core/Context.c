// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Context.h"

#include "ovk/core/Domain.h"
#include "ovk/core/ErrorHandler.h"
#include "ovk/core/Global.h"
#include "ovk/core/List.h"
#include "ovk/core/Logger.h"
#include "ovk/core/TextUtils.h"

ovk_error ovkCreateContext(ovk_context **Context_, const ovk_context_params *Params) {

  int MPIInitialized;
  MPI_Initialized(&MPIInitialized);
  if (!MPIInitialized) {
    fprintf(stderr, "ERROR: MPI not initialized.\n"); fflush(stderr);
    bool Exit;
    if (Params) {
      Exit = Params->error_handler_type == OVK_ERROR_HANDLER_ABORT;
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
  OVK_DEBUG_ASSERT(Context_, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  int CommRank;
  MPI_Comm_rank(Params->comm, &CommRank);

  // Make sure size_t MPI type was detected successfully
  if (KMPI_UNSIGNED_SIZE == MPI_DATATYPE_NULL) {
    if (CommRank == 0) {
      fprintf(stderr, "ERROR: Failed to detect MPI datatype corresponding to size_t.\n");
      fflush(stderr);
    }
    if (Params->error_handler_type == OVK_ERROR_HANDLER_ABORT) {
      MPI_Abort(Params->comm, OVK_ERROR_MPI);
    } else {
      return OVK_ERROR_MPI;
    }
  }

  MPI_Comm Comm;
  MPI_Comm_dup(Params->comm, &Comm);

  MPI_Barrier(Comm);

  *Context_ = malloc(sizeof(ovk_context));
  ovk_context *Context = *Context_;

  Context->comm = Comm;
  MPI_Comm_size(Comm, &Context->comm_size);
  MPI_Comm_rank(Comm, &Context->comm_rank);

  CreateLogger(&Context->logger, Params->log_level, Context->comm_rank);
  CreateErrorHandler(&Context->error_handler, Params->error_handler_type);

  switch (Params->error_handler_type) {
  case OVK_ERROR_HANDLER_ABORT:
    MPI_Comm_set_errhandler(Context->comm, MPI_ERRORS_ARE_FATAL);
    break;
  case OVK_ERROR_HANDLER_RETURN:
    MPI_Comm_set_errhandler(Context->comm, MPI_ERRORS_RETURN);
    break;
  }

  ListCreate(&Context->domains);

  MPI_Barrier(Context->comm);

  if (Context->comm_rank == 0 && LoggingStatus(Context->logger)) {
    char ProcessesString[NUMBER_STRING_LENGTH+10];
    PluralizeLabel(Context->comm_size, "processes", "process", ProcessesString);
    LogStatus(Context->logger, true, 0, "Created context on %s.", ProcessesString);
  }

  return OVK_NO_ERROR;

}

void ovkDestroyContext(ovk_context **Context_) {

  OVK_DEBUG_ASSERT(Context_, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(*Context_, "Invalid context pointer.");

  ovk_context *Context = *Context_;

  MPI_Barrier(Context->comm);

  t_list_entry *Entry = ListBegin(Context->domains);
  while (Entry != ListEnd(Context->domains)) {
    ovk_domain *Domain = ListData(Entry);
    DestroyDomain(&Domain);
    Entry = ListNext(Entry);
  }
  ListDestroy(&Context->domains);

  t_logger *Logger = Context->logger;

  DestroyErrorHandler(&Context->error_handler);

  MPI_Comm Comm = Context->comm;
  bool IsRoot = Context->comm_rank == 0;

  free_null(Context_);

  MPI_Barrier(Comm);

  LogStatus(Logger, IsRoot, 0, "Destroyed context.");
  DestroyLogger(&Logger);

  MPI_Comm_free(&Comm);

}

void ovkGetContextComm(const ovk_context *Context, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Context->comm;

}

void ovkGetContextLogLevel(const ovk_context *Context, ovk_log_level *LogLevel) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(LogLevel, "Invalid log level pointer.");

  GetLogLevel(Context->logger, LogLevel);

}

void ovkSetContextLogLevel(ovk_context *Context, ovk_log_level LogLevel) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(ValidLogLevel(LogLevel), "Invalid log level.");

  SetLogLevel(Context->logger, LogLevel);

}

void ovkGetContextErrorHandlerType(const ovk_context *Context, ovk_error_handler_type
  *ErrorHandlerType) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(ErrorHandlerType, "Invalid error handler type pointer.");

  GetErrorHandlerType(Context->error_handler, ErrorHandlerType);

}

void ovkSetContextErrorHandlerType(ovk_context *Context, ovk_error_handler_type ErrorHandlerType) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(ValidErrorHandlerType(ErrorHandlerType), "Invalid error handler type.");

  SetErrorHandlerType(Context->error_handler, ErrorHandlerType);

}

void ovkCreateDomain(ovk_context *Context, ovk_domain **Domain_, const ovk_domain_params *Params) {

  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");
  OVK_DEBUG_ASSERT(Domain_, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

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

  OVK_DEBUG_ASSERT(Params_, "Invalid params pointer.");

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

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  free_null(Params);

}

void ovkGetContextParamComm(const ovk_context_params *Params, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  *Comm = Params->comm;

}

void ovkSetContextParamComm(ovk_context_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm != MPI_COMM_NULL, "Invalid MPI communicator.");

  Params->comm = Comm;

}

void ovkGetContextParamLogLevel(const ovk_context_params *Params, ovk_log_level *LogLevel) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(LogLevel, "Invalid log level pointer.");

  *LogLevel = Params->log_level;

}

void ovkSetContextParamLogLevel(ovk_context_params *Params, ovk_log_level LogLevel) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(ValidLogLevel(Params->log_level), "Invalid log level.");

  Params->log_level = LogLevel;

}

void ovkGetContextParamErrorHandlerType(const ovk_context_params *Params, ovk_error_handler_type *
  ErrorHandlerType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(ErrorHandlerType, "Invalid error handler type pointer.");

  *ErrorHandlerType = Params->error_handler_type;

}

void ovkSetContextParamErrorHandlerType(ovk_context_params *Params, ovk_error_handler_type
  ErrorHandlerType) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(ValidErrorHandlerType(Params->error_handler_type), "Invalid error handler type.");

  Params->error_handler_type = ErrorHandlerType;

}
