// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_GLOBAL_INCLUDED
#define OVK_EXTRAS_GLOBAL_INCLUDED

#include "ovk/extras/ovkGlobal.h"

#include "ovk/core/Global.h"

// Headers that are used by nearly every source file
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline bool ValidEndian(ovk_ext_endian Endian) {

  switch (Endian) {
  case OVK_EXT_LITTLE_ENDIAN:
  case OVK_EXT_BIG_ENDIAN:
    return true;
  default:
    return false;
  }

}

static inline bool ValidXINTOUTFormat(ovk_ext_xintout_format Format) {

  switch (Format) {
  case OVK_EXT_XINTOUT_STANDARD:
  case OVK_EXT_XINTOUT_EXTENDED:
    return true;
  default:
    return false;
  }

}

#ifdef __cplusplus
}
#endif

#endif
