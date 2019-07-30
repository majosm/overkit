// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_DEBUG_H_LOADED
#define OVK_SUPPORT_DEBUG_H_LOADED

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef OVK_POSIX_SYSTEM
void support_DebuggerAttachHelper();
#endif

#ifdef __cplusplus
}
#endif

#endif
