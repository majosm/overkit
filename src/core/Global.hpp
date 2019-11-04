// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GLOBAL_HPP_INCLUDED
#define OVK_CORE_GLOBAL_HPP_INCLUDED

#include <ovk/core/Global.h>

#include <type_traits>

namespace ovk {

// For detecting accidently-included C++ code in C headers
#ifndef OVK_CXX
#define OVK_CXX
#endif

// Test helper (allows access to class internals in tests)
namespace core {
template <typename T> class test_helper;
}

constexpr int MAX_DIMS = OVK_MAX_DIMS;
constexpr int ALL_GRIDS = OVK_ALL_GRIDS;

using byte = unsigned char;

enum class array_layout : typename std::underlying_type<ovk_array_layout>::type {
  ROW_MAJOR = OVK_ROW_MAJOR,
  COLUMN_MAJOR = OVK_COLUMN_MAJOR
};

inline bool ValidArrayLayout(array_layout Layout) {
  return ovkValidArrayLayout(ovk_array_layout(Layout));
}

enum class endian : typename std::underlying_type<ovk_endian>::type {
  LITTLE = OVK_LITTLE_ENDIAN,
  BIG = OVK_BIG_ENDIAN
};

inline bool ValidEndian(endian Endian) {
  return ovkValidEndian(ovk_endian(Endian));
}

}

#endif
