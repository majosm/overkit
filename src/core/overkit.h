// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_OVERKIT_INCLUDED
#define OVK_CORE_OVERKIT_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include <ovkContext.h>
#include <ovkGlobal.h>

// Make sure internal headers weren't accidentally included into any public headers
#ifdef OVK_INTERNAL
#error Overkit internal error: leaked internal header.
#endif

#ifdef __cplusplus
}
#endif

#endif
