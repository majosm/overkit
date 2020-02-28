// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CART_H_INCLUDED
#define OVK_CORE_CART_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

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
