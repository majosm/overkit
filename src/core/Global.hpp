// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GLOBAL_HPP_INCLUDED
#define OVK_CORE_GLOBAL_HPP_INCLUDED

#include <ovk/core/GlobalBase.h>

namespace ovk {

// For detecting accidently-included C++ code in C headers
#ifndef OVK_CXX
#define OVK_CXX
#endif

// Test helper (allows access to class internals in tests)
namespace core {
template <typename T> class test_helper;
}

}

#endif
