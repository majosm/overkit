// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONSTANTS_BASE_H_INCLUDED
#define OVK_CORE_CONSTANTS_BASE_H_INCLUDED

#include <ovk/core/GlobalBase.h>

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  OVK_MAX_DIMS = 3,
  OVK_ALL_GRIDS = -1
};

typedef enum {
  OVK_LOG_NONE = 0,
  OVK_LOG_ERRORS = 1 << 0,
  OVK_LOG_WARNINGS = 1 << 1,
  OVK_LOG_STATUS = 1 << 2,
  OVK_LOG_DEBUG = 1 << 3,
  OVK_LOG_ALL =
    OVK_LOG_ERRORS |
    OVK_LOG_WARNINGS |
    OVK_LOG_STATUS |
    OVK_LOG_DEBUG
} ovk_log_level;

typedef enum {
  OVK_PERIODIC_STORAGE_UNIQUE,
  OVK_PERIODIC_STORAGE_DUPLICATED
} ovk_periodic_storage;

typedef enum {
  OVK_GEOMETRY_UNIFORM,
  OVK_GEOMETRY_ORIENTED_UNIFORM,
  OVK_GEOMETRY_RECTILINEAR,
  OVK_GEOMETRY_ORIENTED_RECTILINEAR,
  OVK_GEOMETRY_CURVILINEAR
} ovk_geometry_type;

typedef enum {
  OVK_OCCLUDES_NONE,
  OVK_OCCLUDES_ALL,
  OVK_OCCLUDES_COARSE
} ovk_occludes;

typedef enum {
  OVK_CONNECTION_NONE,
  OVK_CONNECTION_NEAREST,
  OVK_CONNECTION_LINEAR,
  OVK_CONNECTION_CUBIC
} ovk_connection_type;

typedef enum {
  OVK_ROW_MAJOR,
  OVK_COLUMN_MAJOR,
  OVK_GRID_LAYOUT = OVK_COLUMN_MAJOR
} ovk_array_layout;

typedef enum {
  OVK_BOOL,
  OVK_BYTE,
  OVK_INT,
  OVK_LONG,
  OVK_LONG_LONG,
  OVK_UNSIGNED_INT,
  OVK_UNSIGNED_LONG,
  OVK_UNSIGNED_LONG_LONG,
  OVK_FLOAT,
  OVK_DOUBLE,
#ifdef INT32_MAX
  OVK_INT32 = sizeof(int) == 4 ? OVK_INT : OVK_LONG,
#endif
#ifdef INT64_MAX
  OVK_INT64 =
    sizeof(int) == 8 ? OVK_INT :
    (sizeof(long) == 8 ? OVK_LONG :
    OVK_LONG_LONG),
#endif
#ifdef UINT32_MAX
  OVK_UNSIGNED_INT32 = sizeof(unsigned int) == 4 ? OVK_UNSIGNED_INT : OVK_UNSIGNED_LONG,
#endif
#ifdef UINT64_MAX
  OVK_UNSIGNED_INT64 =
    sizeof(unsigned int) == 8 ? OVK_UNSIGNED_INT :
    (sizeof(unsigned long) == 8 ? OVK_UNSIGNED_LONG :
    OVK_UNSIGNED_LONG_LONG)
#endif
} ovk_data_type;

typedef enum {
  OVK_COLLECT_NONE,
  OVK_COLLECT_ANY,
  OVK_COLLECT_NOT_ALL,
  OVK_COLLECT_ALL,
  OVK_COLLECT_INTERPOLATE
} ovk_collect_op;

typedef enum {
  OVK_DISPERSE_OVERWRITE,
  OVK_DISPERSE_APPEND
} ovk_disperse_op;

static inline bool ovkValidLogLevel(ovk_log_level LogLevel);
static inline bool ovkValidPeriodicStorage(ovk_periodic_storage PeriodicStorage);
static inline bool ovkValidGeometryType(ovk_geometry_type GeometryType);
static inline bool ovkValidOccludes(ovk_occludes Occludes);
static inline bool ovkValidConnectionType(ovk_connection_type ConnectionType);
static inline bool ovkValidArrayLayout(ovk_array_layout Layout);
static inline bool ovkValidDataType(ovk_data_type DataType);
static inline bool ovkValidCollectOp(ovk_collect_op CollectOp);
static inline bool ovkValidDisperseOp(ovk_disperse_op DisperseOp);

#ifdef __cplusplus
}
#endif

#include <ovk/core/ConstantsBase.inl>

#endif
