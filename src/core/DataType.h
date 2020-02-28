// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DATA_TYPE_H_INCLUDED
#define OVK_CORE_DATA_TYPE_H_INCLUDED

#include <ovk/core/Global.h>

#include <mpi.h>

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_BOOL = 1,
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

static inline bool ovkValidDataType(ovk_data_type DataType) {

  switch (DataType) {
  case OVK_BOOL:
  case OVK_BYTE:
  case OVK_INT:
  case OVK_LONG:
  case OVK_LONG_LONG:
  case OVK_UNSIGNED_INT:
  case OVK_UNSIGNED_LONG:
  case OVK_UNSIGNED_LONG_LONG:
  case OVK_FLOAT:
  case OVK_DOUBLE:
    return true;
  default:
    return false;
  }

}

static inline int ovkDataTypeSize(ovk_data_type DataType) {

  switch (DataType) {
    case OVK_BOOL: return sizeof(byte);
    case OVK_BYTE: return sizeof(byte);
    case OVK_INT: return sizeof(int);
    case OVK_LONG: return sizeof(long);
    case OVK_LONG_LONG: return sizeof(long long);
    case OVK_UNSIGNED_INT: return sizeof(unsigned int);
    case OVK_UNSIGNED_LONG: return sizeof(unsigned long);
    case OVK_UNSIGNED_LONG_LONG: return sizeof(unsigned long long);
    case OVK_FLOAT: return sizeof(float);
    case OVK_DOUBLE: return sizeof(double);
    default: return 0;
  };

}

static inline bool ovkDataTypeIsIntegral(ovk_data_type DataType) {

  switch (DataType) {
  case OVK_INT:
  case OVK_LONG:
  case OVK_LONG_LONG:
  case OVK_BOOL:
  case OVK_BYTE:
  case OVK_UNSIGNED_INT:
  case OVK_UNSIGNED_LONG:
  case OVK_UNSIGNED_LONG_LONG:
    return true;
  default:
    return false;
  }

}

static inline bool ovkDataTypeIsFloatingPoint(ovk_data_type DataType) {

  switch (DataType) {
  case OVK_FLOAT:
  case OVK_DOUBLE:
    return true;
  default:
    return false;
  }

}

static inline bool ovkDataTypeIsSigned(ovk_data_type DataType) {

  switch (DataType) {
  case OVK_INT:
  case OVK_LONG:
  case OVK_LONG_LONG:
  case OVK_FLOAT:
  case OVK_DOUBLE:
    return true;
  default:
    return false;
  }

}

static inline bool ovkDataTypeIsUnsigned(ovk_data_type DataType) {

  switch (DataType) {
  case OVK_BOOL:
  case OVK_BYTE:
  case OVK_UNSIGNED_INT:
  case OVK_UNSIGNED_LONG:
  case OVK_UNSIGNED_LONG_LONG:
    return true;
  default:
    return false;
  }

}

static inline MPI_Datatype ovkDataTypeToMPI(ovk_data_type DataType) {

  switch (DataType) {
  case OVK_BYTE: return MPI_BYTE;
  case OVK_INT: return MPI_INT;
  case OVK_LONG: return MPI_LONG;
  case OVK_LONG_LONG: return MPI_LONG_LONG;
  case OVK_UNSIGNED_INT: return MPI_UNSIGNED;
  case OVK_UNSIGNED_LONG: return MPI_UNSIGNED_LONG;
  case OVK_UNSIGNED_LONG_LONG: return MPI_UNSIGNED_LONG_LONG;
  case OVK_FLOAT: return MPI_FLOAT;
  case OVK_DOUBLE: return MPI_DOUBLE;
  default: return MPI_DATATYPE_NULL;
  }

}

#ifdef __cplusplus
}
#endif

#endif
