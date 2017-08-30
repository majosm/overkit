// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GLOBAL_INCLUDED
#define OVK_CORE_GLOBAL_INCLUDED

#include "ovkGlobal.h"

// Headers that are used by nearly every source file
#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>

// Use this in overkit.h to check whether internal headers were included into any public headers
#ifndef OVK_INTERNAL
#define OVK_INTERNAL
#endif

static const int MAX_DIMS = OVK_MAX_DIMS;

static inline int min(a, b) {
  return a < b ? a : b;
}

static inline int max(a, b) {
  return a > b ? a : b;
}

#endif
