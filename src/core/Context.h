// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONTEXT_INCLUDED
#define OVK_CORE_CONTEXT_INCLUDED

#include "ovkContext.h"

#include "Domain.h"
#include "ErrorHandler.h"
#include "Global.h"
#include "List.h"
#include "Logger.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_context_params {
  MPI_Comm comm;
  ovk_log_level log_level;
  ovk_error_handler_type error_handler_type;
};

struct ovk_context_properties {
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  ovk_log_level log_level;
  ovk_error_handler_type error_handler_type;
};

struct ovk_context {
  ovk_context_properties properties;
  t_logger *logger;
  t_error_handler *error_handler;
  t_list *domains;
};

#ifdef __cplusplus
}
#endif

#endif
