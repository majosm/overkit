// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ERROR_BASE_H_INCLUDED
#define OVK_CORE_ERROR_BASE_H_INCLUDED

#include <ovk/core/GlobalBase.h>

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_ERROR_NONE = 0,
  OVK_ERROR_FILE_OPEN,
  OVK_ERROR_FILE_READ,
  OVK_ERROR_FILE_WRITE
} ovk_error;

static inline bool ovkValidError(ovk_error Error) {

  switch (Error) {
  case OVK_ERROR_NONE:
  case OVK_ERROR_FILE_OPEN:
  case OVK_ERROR_FILE_READ:
  case OVK_ERROR_FILE_WRITE:
    return true;
  default:
    return false;
  }

}

#ifdef __cplusplus
}
#endif

#endif
