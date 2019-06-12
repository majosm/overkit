// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONSTANTS_HPP_INCLUDED
#define OVK_CORE_CONSTANTS_HPP_INCLUDED

#include <ovk/core/ConstantsBase.h>
#include <ovk/core/Global.hpp>

#include <cstdint>

namespace ovk {

constexpr int MAX_DIMS = OVK_MAX_DIMS;
constexpr int ALL_GRIDS = OVK_ALL_GRIDS;

enum class log_level {
  NONE = OVK_LOG_NONE,
  ERRORS = OVK_LOG_ERRORS,
  WARNINGS = OVK_LOG_WARNINGS,
  STATUS = OVK_LOG_STATUS,
  DEBUG = OVK_LOG_DEBUG,
  ALL = OVK_LOG_ALL
};

constexpr inline log_level operator|(log_level Left, log_level Right) {
  return log_level(int(Left) | int(Right));
}
constexpr inline log_level operator&(log_level Left, log_level Right) {
  return log_level(int(Left) & int(Right));
}
constexpr inline log_level operator^(log_level Left, log_level Right) {
  return log_level(int(Left) ^ int(Right));
}
constexpr inline log_level operator~(log_level LogLevel) {
  return log_level(~int(LogLevel));
}
inline log_level operator|=(log_level &Left, log_level Right) {
  return Left = Left | Right;
}
inline log_level operator&=(log_level &Left, log_level Right) {
  return Left = Left & Right;
}
inline log_level operator^=(log_level &Left, log_level Right) {
  return Left = Left ^ Right;
}

enum class periodic_storage {
  UNIQUE = OVK_PERIODIC_STORAGE_UNIQUE,
  DUPLICATED = OVK_PERIODIC_STORAGE_DUPLICATED
};

enum class geometry_type {
  UNIFORM = OVK_GEOMETRY_UNIFORM,
  ORIENTED_UNIFORM = OVK_GEOMETRY_ORIENTED_UNIFORM,
  RECTILINEAR = OVK_GEOMETRY_RECTILINEAR,
  ORIENTED_RECTILINEAR = OVK_GEOMETRY_ORIENTED_RECTILINEAR,
  CURVILINEAR = OVK_GEOMETRY_CURVILINEAR
};

enum class occludes {
  NONE = OVK_OCCLUDES_NONE,
  ALL = OVK_OCCLUDES_ALL,
  COARSE = OVK_OCCLUDES_COARSE
};

enum class connection_type {
  NONE = OVK_CONNECTION_NONE,
  NEAREST = OVK_CONNECTION_NEAREST,
  LINEAR = OVK_CONNECTION_LINEAR,
  CUBIC = OVK_CONNECTION_CUBIC
};

enum class array_layout {
  ROW_MAJOR = OVK_ROW_MAJOR,
  COLUMN_MAJOR = OVK_COLUMN_MAJOR,
  GRID = OVK_GRID_LAYOUT
};

enum class data_type {
  BOOL = OVK_BOOL,
  BYTE = OVK_BYTE,
  INT = OVK_INT,
  LONG = OVK_LONG,
  LONG_LONG = OVK_LONG_LONG,
#ifdef INT32_MAX
  INT32 = OVK_INT32,
#endif
#ifdef INT64_MAX
  INT64 = OVK_INT64,
#endif
  UNSIGNED_INT = OVK_UNSIGNED_INT,
  UNSIGNED_LONG = OVK_UNSIGNED_LONG,
  UNSIGNED_LONG_LONG = OVK_UNSIGNED_LONG_LONG,
#ifdef UINT32_MAX
  UNSIGNED_INT32 = OVK_UNSIGNED_INT32,
#endif
#ifdef UINT64_MAX
  UNSIGNED_INT64 = OVK_UNSIGNED_INT64,
#endif
  FLOAT = OVK_FLOAT,
  DOUBLE = OVK_DOUBLE
};

enum class collect_op {
  NONE = OVK_COLLECT_NONE,
  ANY = OVK_COLLECT_ANY,
  NOT_ALL = OVK_COLLECT_NOT_ALL,
  ALL = OVK_COLLECT_ALL,
  INTERPOLATE = OVK_COLLECT_INTERPOLATE
};

enum class disperse_op {
  OVERWRITE = OVK_DISPERSE_OVERWRITE,
  APPEND = OVK_DISPERSE_APPEND
};

inline bool ValidLogLevel(log_level LogLevel) { return ovkValidLogLevel(ovk_log_level(LogLevel)); }
inline bool ValidPeriodicStorage(periodic_storage PeriodicStorage) { return ovkValidPeriodicStorage(ovk_periodic_storage(PeriodicStorage)); }
inline bool ValidGeometryType(geometry_type GeometryType) { return ovkValidGeometryType(ovk_geometry_type(GeometryType)); }
inline bool ValidOccludes(occludes Occludes) { return ovkValidOccludes(ovk_occludes(Occludes)); }
inline bool ValidConnectionType(connection_type ConnectionType) { return ovkValidConnectionType(ovk_connection_type(ConnectionType)); }
inline bool ValidArrayLayout(array_layout Layout) { return ovkValidArrayLayout(ovk_array_layout(Layout)); }
inline bool ValidDataType(data_type DataType) { return ovkValidDataType(ovk_data_type(DataType)); }
inline bool ValidCollectOp(collect_op CollectOp) { return ovkValidCollectOp(ovk_collect_op(CollectOp)); }
inline bool ValidDisperseOp(disperse_op DisperseOp) { return ovkValidDisperseOp(ovk_disperse_op(DisperseOp)); }

}

#endif
