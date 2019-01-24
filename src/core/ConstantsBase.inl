// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifdef __cplusplus
extern "C" {
#endif

static inline bool ovkValidLogLevel(ovk_log_level LogLevel) {

  return LogLevel >= OVK_LOG_NONE && LogLevel <= OVK_LOG_ALL;

}

static inline bool ovkValidErrorHandlerType(ovk_error_handler_type ErrorHandlerType) {

  switch (ErrorHandlerType) {
  case OVK_ERROR_HANDLER_ABORT:
  case OVK_ERROR_HANDLER_RETURN:
    return true;
  default:
    return false;
  }

}

static inline bool ovkValidError(ovk_error Error) {

  return Error >= OVK_NO_ERROR && Error < OVK_MAX_ERROR;

}

static inline bool ovkValidDomainConfig(ovk_domain_config DomainConfig) {

  return DomainConfig >= OVK_DOMAIN_CONFIG_NONE && DomainConfig <= OVK_DOMAIN_CONFIG_ALL;

}

static inline bool ovkValidPeriodicStorage(ovk_periodic_storage PeriodicStorage) {

  switch (PeriodicStorage) {
  case OVK_PERIODIC_STORAGE_UNIQUE:
  case OVK_PERIODIC_STORAGE_DUPLICATED:
    return true;
  default:
    return false;
  }

}

static inline bool ovkValidGeometryType(ovk_geometry_type GeometryType) {

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

static inline bool ovkValidOccludes(ovk_occludes Occludes) {

  switch (Occludes) {
  case OVK_OCCLUDES_NONE:
  case OVK_OCCLUDES_ALL:
  case OVK_OCCLUDES_COARSE:
    return true;
  default:
    return false;
  }

}

static inline bool ovkValidConnectionType(ovk_connection_type ConnectionType) {

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

static inline bool ovkValidArrayLayout(ovk_array_layout Layout) {

  switch (Layout) {
  case OVK_ROW_MAJOR:
  case OVK_COLUMN_MAJOR:
    return true;
  default:
    return false;
  }

}

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

static inline bool ovkValidCollectOp(ovk_collect_op CollectOp) {

  switch (CollectOp) {
  case OVK_COLLECT_NONE:
  case OVK_COLLECT_ANY:
  case OVK_COLLECT_NOT_ALL:
  case OVK_COLLECT_ALL:
  case OVK_COLLECT_INTERPOLATE:
    return true;
  default:
    return false;
  }

}

static inline bool ovkValidDisperseOp(ovk_disperse_op DisperseOp) {

  switch (DisperseOp) {
  case OVK_DISPERSE_OVERWRITE:
  case OVK_DISPERSE_APPEND:
    return true;
  default:
    return false;
  }

}

#ifdef __cplusplus
}
#endif
