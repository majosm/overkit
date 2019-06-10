// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_EXCHANGER_BASE_H_INCLUDED
#define OVK_CORE_EXCHANGER_BASE_H_INCLUDED

#include <ovk/core/GlobalBase.h>

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_COLLECT_NONE,
  OVK_COLLECT_ANY,
  OVK_COLLECT_NOT_ALL,
  OVK_COLLECT_ALL,
  OVK_COLLECT_INTERPOLATE
} ovk_collect_op;

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

typedef enum {
  OVK_DISPERSE_OVERWRITE,
  OVK_DISPERSE_APPEND
} ovk_disperse_op;

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

#endif
