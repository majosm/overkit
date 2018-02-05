// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_PUBLIC_GLOBAL_INCLUDED
#define OVK_EXTRAS_PUBLIC_GLOBAL_INCLUDED

#include <ovk/core/ovkGlobal.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_EXT_LITTLE_ENDIAN,
  OVK_EXT_BIG_ENDIAN
} ovk_ext_endian;

typedef enum {
  OVK_EXT_XINTOUT_STANDARD,
  OVK_EXT_XINTOUT_EXTENDED
} ovk_ext_xintout_format;

#ifdef __cplusplus
}
#endif

#endif
