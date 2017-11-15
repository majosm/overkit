// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_GLOBAL_INCLUDED
#define OVK_CORE_PUBLIC_GLOBAL_INCLUDED

// Headers that are needed for API
#include <mpi.h>
#include <stdbool.h>
#include <stdint.h>

#ifdef OVERKIT_DEBUG
  static const bool OVK_DEBUG = true;
#else
  static const bool OVK_DEBUG = false;
#endif

static const int OVK_MAX_DIMS = 3;

typedef enum {
  OVK_ERROR_HANDLER_ABORT,
  OVK_ERROR_HANDLER_RETURN
} ovk_error_handler_type;

typedef enum {
  OVK_ERROR_NONE = 0,
  OVK_MAX_ERROR
} ovk_error;

#endif
