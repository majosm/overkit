// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVERKIT_H_INCLUDED
#define OVERKIT_H_INCLUDED

#include <ovk/core-c/all.h>
#include <ovk/extras-c/all.h>

// Make sure C++ headers weren't accidentally included
#ifdef OVK_CXX
#error Overkit internal error: leaked C++ header.
#endif

#endif
