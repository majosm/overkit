// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVERKIT_EXTRAS_INCLUDED
#define OVERKIT_EXTRAS_INCLUDED

#include <ovk/extras/ovkGlobal.h>

#ifdef OVERKIT_XPACC
#include <ovk/extras/ovkXINTOUT.h>
#endif

// Make sure internal headers weren't accidentally included into any public headers
#ifdef OVK_INTERNAL
#error Overkit internal error: leaked internal header.
#endif
#ifdef OVK_EXTRAS_INTERNAL
#error Overkit internal error: leaked internal header.
#endif

#endif
