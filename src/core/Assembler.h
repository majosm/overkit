// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ASSEMBLER_H_INCLUDED
#define OVK_CORE_ASSEMBLER_H_INCLUDED

#include <ovk/core/Global.h>

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_OCCLUDES_NONE,
  OVK_OCCLUDES_ALL,
  OVK_OCCLUDES_COARSE
} ovk_occludes;

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

typedef enum {
  OVK_CONNECTION_NONE,
  OVK_CONNECTION_NEAREST,
  OVK_CONNECTION_LINEAR,
  OVK_CONNECTION_CUBIC
} ovk_connection_type;

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

#ifdef __cplusplus
}
#endif

#endif
