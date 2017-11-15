// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DEBUG_INCLUDED
#define OVK_CORE_DEBUG_INCLUDED

#include "Global.h"

#include <stdarg.h>
#include <stdio.h>

#ifdef OVERKIT_DEBUG
#define OVK_DEBUG_ASSERT(Condition, ...) \
  if (!(Condition)) DebugAssert(__FILE__, __LINE__, __VA_ARGS__);
#else
#define OVK_DEBUG_ASSERT(...)
#endif

static inline void DebugAssert(const char *File, int Line, const char *Format, ...) {

  va_list ArgList;

  fprintf(stderr, "DEBUG ERROR (line %i of file '%s'): ", Line, File);
  va_start(ArgList, Format);
  vfprintf(stderr, Format, ArgList);
  fprintf(stderr, "\n");
  va_end(ArgList);
  fflush(stderr);

  MPI_Abort(MPI_COMM_WORLD, 1);

}

#endif
