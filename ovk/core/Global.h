// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GLOBAL_INCLUDED
#define OVK_CORE_GLOBAL_INCLUDED

#include "ovk/core/ovkGlobal.h"

// Headers that are used by nearly every source file
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>

#ifndef __cplusplus
#undef min
#undef max
#define min(a, b) ovk_min(a, b)
#define max(a, b) ovk_max(a, b)
#endif
#define clamp(a, b, c) ovk_clamp(a, b, c)

// Use this in overkit.h to check whether internal headers were included into any public headers
#ifndef OVK_INTERNAL
#define OVK_INTERNAL
#endif

#define PRIVATE(Func) OVK_PRIVATE(Func)

#define free_null(Ptr) do { free(*Ptr); *Ptr = NULL; } while (false)

#ifdef __cplusplus
extern "C" {
#endif

enum {
  MAX_DIMS = OVK_MAX_DIMS,
  NUMBER_STRING_LENGTH = 32
};

// Define MPI datatype for size_t
static const MPI_Datatype KMPI_UNSIGNED_SIZE =
  sizeof(size_t) == sizeof(unsigned int) ? MPI_UNSIGNED :
  (sizeof(size_t) == sizeof(unsigned long) ? MPI_UNSIGNED_LONG :
  (sizeof(size_t) == sizeof(unsigned long long) ? MPI_UNSIGNED_LONG_LONG :
  MPI_DATATYPE_NULL));

static inline bool ValidLogLevel(ovk_log_level LogLevel) {

  return LogLevel >= OVK_LOG_NONE && LogLevel <= OVK_LOG_ALL;

}

static inline bool ValidErrorHandlerType(ovk_error_handler_type ErrorHandlerType) {

  switch (ErrorHandlerType) {
  case OVK_ERROR_HANDLER_ABORT:
  case OVK_ERROR_HANDLER_RETURN:
    return true;
  default:
    return false;
  }

}

static inline bool ValidError(ovk_error Error) {

  return Error >= OVK_NO_ERROR && Error < OVK_MAX_ERROR;

}

static inline bool ValidDomainConfig(ovk_domain_config DomainConfig) {

  return DomainConfig >= OVK_DOMAIN_CONFIG_NONE && DomainConfig <= OVK_DOMAIN_CONFIG_ALL;

}

static inline bool ValidPeriodicStorage(ovk_periodic_storage PeriodicStorage) {

  switch (PeriodicStorage) {
  case OVK_PERIODIC_STORAGE_UNIQUE:
  case OVK_PERIODIC_STORAGE_DUPLICATED:
    return true;
  default:
    return false;
  }

}

static inline bool ValidGeometryType(ovk_geometry_type GeometryType) {

  switch (GeometryType) {
  case OVK_GEOMETRY_UNIFORM:
  case OVK_GEOMETRY_ORIENTED_UNIFORM:
  case OVK_GEOMETRY_RECTILINEAR:
  case OVK_GEOMETRY_ORIENTED_RECTILINEAR:
  case OVK_GEOMETRY_CURVILINEAR:
    return true;
  default:
    return false;
  }

}

static inline bool ValidOccludes(ovk_occludes Occludes) {

  switch (Occludes) {
  case OVK_OCCLUDES_NONE:
  case OVK_OCCLUDES_ALL:
  case OVK_OCCLUDES_COARSE:
    return true;
  default:
    return false;
  }

}

static inline bool ValidConnectionType(ovk_connection_type ConnectionType) {

  switch (ConnectionType) {
  case OVK_CONNECTION_NONE:
  case OVK_CONNECTION_NEAREST:
  case OVK_CONNECTION_LINEAR:
  case OVK_CONNECTION_CUBIC:
    return true;
  default:
    return false;
  }

}

static inline bool ValidArrayLayout(ovk_array_layout Layout) {

  switch (Layout) {
  case OVK_ROW_MAJOR:
  case OVK_COLUMN_MAJOR:
    return true;
  default:
    return false;
  }

}

static inline bool ValidDataType(ovk_data_type DataType) {

  switch (DataType) {
  case OVK_BOOL:
  case OVK_BYTE:
  case OVK_INT32:
  case OVK_INT64:
  case OVK_UINT32:
  case OVK_UINT64:
  case OVK_FLOAT:
  case OVK_DOUBLE:
    return true;
  default:
    return false;
  }

}

static inline int DataTypeSize(ovk_data_type DataType) {

  switch (DataType) {
    case OVK_BOOL: return sizeof(unsigned char);
    case OVK_BYTE: return sizeof(unsigned char);
    case OVK_INT: return sizeof(int);
    case OVK_LONG: return sizeof(long);
    case OVK_LONG_LONG: return sizeof(long long);
    case OVK_INT32: return 4;
    case OVK_INT64: return 8;
    case OVK_UNSIGNED_INT: return sizeof(unsigned int);
    case OVK_UNSIGNED_LONG: return sizeof(unsigned long);
    case OVK_UNSIGNED_LONG_LONG: return sizeof(unsigned long long);
    case OVK_UNSIGNED_SIZE: return sizeof(size_t);
    case OVK_UINT32: return 4;
    case OVK_UINT64: return 8;
    case OVK_FLOAT: return sizeof(float);
    case OVK_DOUBLE: return sizeof(double);
    default: return 0;
  };

}

static inline MPI_Datatype DataTypeToMPI(ovk_data_type DataType) {

  switch (DataType) {
  case OVK_BOOL: return MPI_BYTE;
  case OVK_BYTE: return MPI_BYTE;
  case OVK_INT: return MPI_INT;
  case OVK_LONG: return MPI_LONG;
  case OVK_LONG_LONG: return MPI_LONG_LONG;
  case OVK_INT32: return MPI_INT;
  case OVK_INT64: return MPI_LONG_LONG;
  case OVK_UNSIGNED_INT: return MPI_UNSIGNED;
  case OVK_UNSIGNED_LONG: return MPI_UNSIGNED_LONG;
  case OVK_UNSIGNED_LONG_LONG: return MPI_UNSIGNED_LONG_LONG;
  case OVK_UNSIGNED_SIZE: return KMPI_UNSIGNED_SIZE;
  case OVK_UINT32: return MPI_UNSIGNED;
  case OVK_UINT64: return MPI_UNSIGNED_LONG_LONG;
  case OVK_FLOAT: return MPI_FLOAT;
  case OVK_DOUBLE: return MPI_DOUBLE;
  default: return MPI_DATATYPE_NULL;
  }

}

static inline bool ValidCollectOp(ovk_collect_op CollectOp) {

  switch (CollectOp) {
  case OVK_COLLECT_INTERPOLATE:
    return true;
  default:
    return false;
  }

}

static inline bool ValidDisperseOp(ovk_disperse_op DisperseOp) {

  switch (DisperseOp) {
  case OVK_DISPERSE_OVERWRITE:
    return true;
  default:
    return false;
  }

}

#ifdef __cplusplus
}
#endif

#endif
