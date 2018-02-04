// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVERKIT_INCLUDED
#define OVERKIT_INCLUDED

#include <ovk/core/ovkBox.h>
#include <ovk/core/ovkBox.inl>
#include <ovk/core/ovkCart.h>
#include <ovk/core/ovkCart.inl>
#include <ovk/core/ovkConnectivityD.h>
#include <ovk/core/ovkConnectivityR.h>
#include <ovk/core/ovkConnectivity.h>
#include <ovk/core/ovkContext.h>
#include <ovk/core/ovkDomain.h>
#include <ovk/core/ovkExchange.h>
#include <ovk/core/ovkGlobal.h>
#include <ovk/core/ovkGrid.h>
#include <ovk/core/ovkRange.h>
#include <ovk/core/ovkRange.inl>

// Make sure internal headers weren't accidentally included into any public headers
#ifdef OVK_INTERNAL
#error Overkit internal error: leaked internal header.
#endif

#endif
