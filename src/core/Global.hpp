// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GLOBAL_HPP_INCLUDED
#define OVK_CORE_GLOBAL_HPP_INCLUDED

#include <ovk/core/Global.h>

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

enum class array_layout {
  ROW_MAJOR = OVK_ROW_MAJOR,
  COLUMN_MAJOR = OVK_COLUMN_MAJOR,
  GRID = OVK_GRID_LAYOUT
};

inline bool ValidArrayLayout(array_layout Layout) {
  return ovkValidArrayLayout(ovk_array_layout(Layout));
}

enum class endian {
  LITTLE = OVK_LITTLE_ENDIAN,
  BIG = OVK_BIG_ENDIAN
};

inline bool ValidEndian(endian Endian) {
  return ovkValidEndian(ovk_endian(Endian));
}

enum class periodic_storage {
  UNIQUE = OVK_PERIODIC_STORAGE_UNIQUE,
  DUPLICATED = OVK_PERIODIC_STORAGE_DUPLICATED
};

inline bool ValidPeriodicStorage(periodic_storage PeriodicStorage) {
  return ovkValidPeriodicStorage(ovk_periodic_storage(PeriodicStorage));
}

}

#endif
