// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_PUBLIC_GLOBAL_INCLUDED
#define OVK_EXTRAS_PUBLIC_GLOBAL_INCLUDED

#include <overkit.h>

#include <mpi.h>

// Apply prefix to names of internal functions that can't be defined as static due to being
// shared between multiple source files
#ifndef OVK_EXTRAS_PRIVATE_PREFIX
#define OVK_EXTRAS_PRIVATE_PREFIX ovkEXT_INTERNAL_
#endif
#define OVK_EXTRAS_PRIVATE(Func) OVK_CONCAT(OVK_EXTRAS_PRIVATE_PREFIX, Func)

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif
