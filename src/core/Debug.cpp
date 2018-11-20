// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Debug.hpp"

#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <cstdarg>
#include <cstdio>

extern "C" {

void ovk_core_DebugAssert(const char *File, int Line, const char *Format, ...) {

  va_list ArgList;

  fprintf(stderr, "DEBUG ERROR (line %i of file '%s'): ", Line, File);
  va_start(ArgList, Format);
  vfprintf(stderr, Format, ArgList);
  fprintf(stderr, "\n");
  va_end(ArgList);
  fflush(stderr);

  MPI_Abort(MPI_COMM_WORLD, 1);

}

}
