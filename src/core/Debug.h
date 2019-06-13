// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DEBUG_H_INCLUDED
#define OVK_CORE_DEBUG_H_INCLUDED

#include <ovk/core/Global.h>

#ifdef __cplusplus
extern "C" {
#endif

#if OVK_DEBUG
void ovk_core_DebugExit(const char *File, int Line, const char *Format, ...);
#define OVK_DEBUG_ASSERT_C(Condition, ...) \
  if (!(Condition)) ovk_core_DebugExit(__FILE__, __LINE__, __VA_ARGS__)
#define OVK_DEBUG_ASSERT OVK_DEBUG_ASSERT_C
#else
#define OVK_DEBUG_ASSERT_C(...)
#define OVK_DEBUG_ASSERT(...)
#endif

#ifdef __cplusplus
}
#endif

#endif
