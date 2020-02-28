// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_DEBUG_HPP_LOADED
#define OVK_SUPPORT_DEBUG_HPP_LOADED

#include <support/Debug.h>

namespace support {

#ifdef OVK_POSIX_SYSTEM
void DebuggerAttachHelper();
#endif

}

#endif
