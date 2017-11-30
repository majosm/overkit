// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_GLOBAL_INCLUDED
#define OVK_CORE_PUBLIC_GLOBAL_INCLUDED

// Headers that are needed for API
#include <mpi.h>
#include <stdbool.h>
#include <stdint.h>

#define ovk_min(a, b) ((a) < (b) ? (a) : (b))
#define ovk_max(a, b) ((a) > (b) ? (a) : (b))

#ifdef OVERKIT_DEBUG
  static const bool OVK_DEBUG = true;
#else
  static const bool OVK_DEBUG = false;
#endif

static const int OVK_MAX_DIMS = 3;
static const int OVK_NAME_LENGTH = 256;

typedef enum {
  OVK_LOG_NONE = 0,
  OVK_LOG_ERRORS = 1 << 0,
  OVK_LOG_WARNINGS = 1 << 1,
  OVK_LOG_STATUS = 1 << 2,
  OVK_LOG_ALL = OVK_LOG_ERRORS | OVK_LOG_WARNINGS | OVK_LOG_STATUS
} ovk_log_level;

typedef enum {
  OVK_ERROR_HANDLER_ABORT,
  OVK_ERROR_HANDLER_RETURN
} ovk_error_handler_type;

typedef enum {
  OVK_ERROR_NONE = 0,
  OVK_MAX_ERROR
} ovk_error;

typedef enum {
  OVK_DOMAIN_CONFIG_NONE = 0,
  OVK_DOMAIN_CONFIG_GEOMETRY = 1 << 0,
  OVK_DOMAIN_CONFIG_OVERLAP = 1 << 1,
  OVK_DOMAIN_CONFIG_CONNECTIVITY = 1 << 2,
  OVK_DOMAIN_CONFIG_EXCHANGE = 1 << 3,
} ovk_domain_config;

static const int OVK_GENERATE_ID = -1;

typedef enum {
  OVK_NO_OVERLAP_PERIODIC,
  OVK_OVERLAP_PERIODIC
} ovk_periodic_storage;

typedef enum {
  OVK_GEOMETRY_TYPE_CARTESIAN,
  OVK_GEOMETRY_TYPE_RECTILINEAR,
  OVK_GEOMETRY_TYPE_ORIENTED_CARTESIAN,
  OVK_GEOMETRY_TYPE_ORIENTED_RECTILINEAR,
  OVK_GEOMETRY_TYPE_CURVILINEAR
} ovk_geometry_type;

typedef enum {
  OVK_BOOL,
  OVK_INT,
  OVK_FLOAT,
  OVK_DOUBLE
} ovk_data_type;

typedef enum {
  OVK_COLLECT_INTERPOLATE
} ovk_collect_op;

typedef enum {
  OVK_DISPERSE_OVERWRITE
} ovk_disperse_op;

#endif
