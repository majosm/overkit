// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GLOBAL_INCLUDED
#define OVK_CORE_GLOBAL_INCLUDED

#include "ovkGlobal.h"

// Headers that are used by nearly every source file
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>

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

#define PRIVATE(func) OVK_PRIVATE(func)

#define free_null(ptr) do { free(*ptr); *ptr = NULL; } while (0)

#ifdef __cplusplus
extern "C" {
#endif

enum {
  MAX_DIMS = OVK_MAX_DIMS,
  NUMBER_STRING_LENGTH = 32
};

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

  return Error >= OVK_ERROR_NONE && Error < OVK_MAX_ERROR;

}

static inline bool ValidDomainConfig(ovk_domain_config DomainConfig) {

  return DomainConfig >= OVK_DOMAIN_CONFIG_NONE && DomainConfig <= OVK_DOMAIN_CONFIG_ALL;

}

static inline bool ValidPeriodicStorage(ovk_periodic_storage PeriodicStorage) {

  return PeriodicStorage == OVK_NO_OVERLAP_PERIODIC || PeriodicStorage == OVK_OVERLAP_PERIODIC;

}

static inline bool ValidGeometryType(ovk_geometry_type GeometryType) {

  switch (GeometryType) {
  case OVK_GEOMETRY_TYPE_CARTESIAN:
  case OVK_GEOMETRY_TYPE_RECTILINEAR:
  case OVK_GEOMETRY_TYPE_ORIENTED_CARTESIAN:
  case OVK_GEOMETRY_TYPE_ORIENTED_RECTILINEAR:
  case OVK_GEOMETRY_TYPE_CURVILINEAR:
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
    case OVK_BOOL: return 1;
    case OVK_BYTE: return 1;
    case OVK_INT32: return 4;
    case OVK_INT64: return 8;
    case OVK_UINT32: return 4;
    case OVK_UINT64: return 8;
    case OVK_FLOAT: return 4;
    case OVK_DOUBLE: return 8;
    default: return 0;
  };

}

static inline void DataTypeToMPI(ovk_data_type DataType, MPI_Datatype *MPIDataType,
  int *MPIDataSize) {

  switch (DataType) {
  case OVK_BOOL:
    *MPIDataType = MPI_BYTE;
    *MPIDataSize = 1;
    return;
  case OVK_BYTE:
    *MPIDataType = MPI_BYTE;
    *MPIDataSize = 1;
    return;
  case OVK_INT32:
    *MPIDataType = MPI_INT;
    *MPIDataSize = 4;
    return;
  case OVK_INT64:
    *MPIDataType = MPI_LONG_LONG;
    *MPIDataSize = 8;
    return;
  case OVK_UINT32:
    *MPIDataType = MPI_UNSIGNED;
    *MPIDataSize = 4;
    return;
  case OVK_UINT64:
  // TODO: Make this work (may need to create new MPI datatype in domain and pass it in)
//     *MPIDataType = MY_MPI_UNSIGNED_LONG_LONG;
//     *MPIDataSize = 8;
    OVK_DEBUG_ASSERT(false, "OVK_UINT64 not implemented yet.");
    return;
  case OVK_FLOAT:
    *MPIDataType = MPI_FLOAT;
    *MPIDataSize = 4;
    return;
  case OVK_DOUBLE:
    *MPIDataType = MPI_DOUBLE;
    *MPIDataSize = 8;
    return;
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
