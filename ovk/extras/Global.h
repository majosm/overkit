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

// Use this in overkit-extras.h to check whether internal headers were included into any public
// headers
#ifndef OVK_EXTRAS_INTERNAL
#define OVK_EXTRAS_INTERNAL
#endif

#undef PRIVATE
#define PRIVATE(Func) OVK_EXTRAS_PRIVATE(Func)

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif
