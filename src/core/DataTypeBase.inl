// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifdef __cplusplus
extern "C" {
#endif

static inline int ovkDataTypeSize(ovk_data_type DataType) {

  switch (DataType) {
    case OVK_BOOL: return sizeof(unsigned char);
    case OVK_BYTE: return sizeof(unsigned char);
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
