// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_OVERKIT_INCLUDED
#define OVK_CORE_OVERKIT_INCLUDED

#include <ovkBox.h>
#include <ovkBox.inl>
#include <ovkCart.h>
#include <ovkCart.inl>
#include <ovkConnectivityD.h>
#include <ovkConnectivityR.h>
#include <ovkConnectivity.h>
#include <ovkContext.h>
#include <ovkDomain.h>
#include <ovkExchange.h>
#include <ovkGlobal.h>
#include <ovkGrid.h>
#include <ovkRange.h>
#include <ovkRange.inl>

// Make sure internal headers weren't accidentally included into any public headers
#ifdef OVK_INTERNAL
#error Overkit internal error: leaked internal header.
#endif

#endif
