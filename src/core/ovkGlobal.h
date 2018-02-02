// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_GLOBAL_INCLUDED
#define OVK_CORE_PUBLIC_GLOBAL_INCLUDED

// Headers that are needed for API
#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define ovk_min(a, b) ((a) < (b) ? (a) : (b))
#define ovk_max(a, b) ((a) > (b) ? (a) : (b))
#define ovk_clamp(a, b, c) (ovk_min(ovk_max((a), (b)), (c)))

// Apply prefix to names of internal functions that can't be defined as static due to being
// shared between multiple source files
#ifndef OVK_PRIVATE_PREFIX
#define OVK_PRIVATE_PREFIX ovk_INTERNAL_
#endif
// Hack to force evaluation of prefix before concatenating tokens
#define OVK_CONCAT_(Left, Right) Left ## Right
#define OVK_CONCAT(Left, Right) OVK_CONCAT_(Left, Right)
#define OVK_PRIVATE(Func) OVK_CONCAT(OVK_PRIVATE_PREFIX, Func)

#ifdef __cplusplus
extern "C" {
#endif

#ifdef OVERKIT_DEBUG
#define OVK_DEBUG true
void OVK_PRIVATE(DebugAssert)(const char *File, int Line, const char *Format, ...);
#define OVK_DEBUG_ASSERT(Condition, ...) \
  if (!(Condition)) OVK_PRIVATE(DebugAssert)(__FILE__, __LINE__, __VA_ARGS__);
#else
#define OVK_DEBUG false
#define OVK_DEBUG_ASSERT(...)
#endif

enum {
  OVK_MAX_DIMS = 3,
  OVK_NAME_LENGTH = 256,
  OVK_ALL_GRIDS = -1
};

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
  OVK_DOMAIN_CONFIG_ALL = OVK_DOMAIN_CONFIG_GEOMETRY | OVK_DOMAIN_CONFIG_OVERLAP |
    OVK_DOMAIN_CONFIG_CONNECTIVITY | OVK_DOMAIN_CONFIG_EXCHANGE
} ovk_domain_config;

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
  OVK_ROW_MAJOR,
  OVK_COLUMN_MAJOR,
  OVK_GRID_LAYOUT = OVK_COLUMN_MAJOR
} ovk_array_layout;

typedef enum {
  OVK_BOOL,
  OVK_BYTE,
  OVK_INT32,
  OVK_INT64,
  OVK_UINT32,
  OVK_UINT64,
  OVK_FLOAT,
  OVK_DOUBLE
} ovk_data_type;

typedef enum {
  OVK_COLLECT_INTERPOLATE
} ovk_collect_op;

typedef enum {
  OVK_DISPERSE_OVERWRITE
} ovk_disperse_op;

#ifdef __cplusplus
}
#endif

#endif
