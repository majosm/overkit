// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_CONSTANTS_BASE_H_INCLUDED
#define OVK_EXTRAS_CONSTANTS_BASE_H_INCLUDED

#include <ovk/extras/GlobalBase.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  OVK_LITTLE_ENDIAN,
  OVK_BIG_ENDIAN
} ovk_endian;

#ifdef OVK_XPACC
typedef enum {
  OVK_XINTOUT_STANDARD,
  OVK_XINTOUT_EXTENDED
} ovk_xintout_format;
#endif

static inline bool ovkValidEndian(ovk_endian Endian);
#ifdef OVK_XPACC
static inline bool ovkValidXINTOUTFormat(ovk_xintout_format Format);
#endif

#ifdef __cplusplus
}
#endif

#include <ovk/extras/ConstantsBase.inl>

#endif
