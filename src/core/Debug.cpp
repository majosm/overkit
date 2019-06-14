// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Debug.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <cstdarg>
#include <cstdio>
#include <string>

extern "C" {

#if OVK_DEBUG

// Wrapper around ovk::core::DebugExit that can be called from C code
void ovk_core_DebugExit(const char *File, int Line, const char *Format, ...) {

  va_list ArgList1;
  va_start(ArgList1, Format);

  va_list ArgList2;
  va_copy(ArgList2, ArgList1);

  int NumChars = vsnprintf(nullptr, 0, Format, ArgList1);

  ovk::array<char> MessageChars({NumChars+1});

  va_end(ArgList1);

  vsnprintf(MessageChars.Data(), MessageChars.Count(), Format, ArgList2);

  va_end(ArgList2);

  std::string Message(MessageChars.LinearBegin(), MessageChars.LinearEnd());

  ovk::core::DebugExit(File, Line, Message);

}

#endif

}
