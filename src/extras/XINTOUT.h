// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_XINTOUT_H_INCLUDED
#define OVK_EXTRAS_XINTOUT_H_INCLUDED

#include <ovk/extras/Global.h>

typedef enum {
  OVK_XINTOUT_STANDARD,
  OVK_XINTOUT_EXTENDED
} ovk_xintout_format;

static inline bool ovkValidXINTOUTFormat(ovk_xintout_format Format) {

  switch (Format) {
  case OVK_XINTOUT_STANDARD:
  case OVK_XINTOUT_EXTENDED:
    return true;
  default:
    return false;
  }

}

#endif
