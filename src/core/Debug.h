// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DEBUG_H_INCLUDED
#define OVK_CORE_DEBUG_H_INCLUDED

#include <ovk/core/Global.h>

#ifdef __cplusplus
extern "C" {
#endif

// Global flag for debugging in hard-to-reach places
extern int ovk_core_DebugFlag;

#if OVK_DEBUG
void ovk_core_DebugExit(bool HasMPI, const char *File, int Line, const char *Format, ...);
#define OVK_DEBUG_ASSERT_C(Condition, ...) \
  if (!(Condition)) ovk_core_DebugExit(true, __FILE__, __LINE__, __VA_ARGS__)
#define OVK_DEBUG_ASSERT_NO_MPI_C(Condition, ...) \
  if (!(Condition)) ovk_core_DebugExit(false, __FILE__, __LINE__, __VA_ARGS__)
#define OVK_DEBUG_ASSERT OVK_DEBUG_ASSERT_C
#define OVK_DEBUG_ASSERT_NO_MPI OVK_DEBUG_ASSERT_NO_MPI_C
#else
#define OVK_DEBUG_ASSERT_C(...)
#define OVK_DEBUG_ASSERT_NO_MPI_C(...)
#define OVK_DEBUG_ASSERT(...)
#define OVK_DEBUG_ASSERT_NO_MPI(...)
#endif

#ifdef __cplusplus
}
#endif

#endif
