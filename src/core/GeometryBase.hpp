// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GEOMETRY_BASE_HPP_INCLUDED
#define OVK_CORE_GEOMETRY_BASE_HPP_INCLUDED

#include <ovk/core/DataType.hpp>
#include <ovk/core/Geometry.h>
#include <ovk/core/Global.hpp>

#include <type_traits>

namespace ovk {

enum class geometry_type : typename std::underlying_type<ovk_geometry_type>::type {
  UNIFORM = OVK_GEOMETRY_TYPE_UNIFORM,
  ORIENTED_UNIFORM = OVK_GEOMETRY_TYPE_ORIENTED_UNIFORM,
  RECTILINEAR = OVK_GEOMETRY_TYPE_RECTILINEAR,
  ORIENTED_RECTILINEAR = OVK_GEOMETRY_TYPE_ORIENTED_RECTILINEAR,
  CURVILINEAR = OVK_GEOMETRY_TYPE_CURVILINEAR
};

inline bool ValidGeometryType(geometry_type GeometryType) {
  return ovkValidGeometryType(ovk_geometry_type(GeometryType));
}

namespace core {
template <> struct data_type_traits<geometry_type> : data_type_traits<typename std::underlying_type<
  geometry_type>::type> {};
}

}

#endif
