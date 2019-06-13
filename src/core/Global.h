// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GLOBAL_H_INCLUDED
#define OVK_CORE_GLOBAL_H_INCLUDED

#include <ovk/core/Config.h>

#include <stdbool.h>

#define OVK_FORCE_INLINE __attribute__((always_inline)) inline

#define ovk_min(a, b) ((a) < (b) ? (a) : (b))
#define ovk_max(a, b) ((a) > (b) ? (a) : (b))

#ifdef __cplusplus
extern "C" {
#endif

enum {
  OVK_MAX_DIMS = 3,
  OVK_ALL_GRIDS = -1
};

typedef enum {
  OVK_ROW_MAJOR,
  OVK_COLUMN_MAJOR,
  OVK_GRID_LAYOUT = OVK_COLUMN_MAJOR
} ovk_array_layout;

static inline bool ovkValidArrayLayout(ovk_array_layout Layout) {

  switch (Layout) {
  case OVK_ROW_MAJOR:
  case OVK_COLUMN_MAJOR:
    return true;
  default:
    return false;
  }

}

typedef enum {
  OVK_LITTLE_ENDIAN,
  OVK_BIG_ENDIAN
} ovk_endian;

static inline bool ovkValidEndian(ovk_endian Endian) {

  switch (Endian) {
  case OVK_LITTLE_ENDIAN:
  case OVK_BIG_ENDIAN:
    return true;
  default:
    return false;
  }

}

typedef enum {
  OVK_PERIODIC_STORAGE_UNIQUE,
  OVK_PERIODIC_STORAGE_DUPLICATED
} ovk_periodic_storage;

static inline bool ovkValidPeriodicStorage(ovk_periodic_storage PeriodicStorage) {

  switch (PeriodicStorage) {
  case OVK_PERIODIC_STORAGE_UNIQUE:
  case OVK_PERIODIC_STORAGE_DUPLICATED:
    return true;
  default:
    return false;
  }

}

#ifdef __cplusplus
}
#endif

#endif
