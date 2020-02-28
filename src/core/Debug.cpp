// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Debug.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <cstdarg>
#include <cstdio>
#include <string>

extern "C" {
int ovk_core_DebugFlag;
}

namespace ovk {
namespace core {
int &DebugFlag = ovk_core_DebugFlag;
}}

#if OVK_DEBUG

extern "C" {

// Wrapper around ovk::core::DebugExit that can be called from C code
void ovk_core_DebugExit(bool HasMPI, const char *File, int Line, const char *Format, ...) {

  std::va_list ArgList1;
  va_start(ArgList1, Format);

  std::va_list ArgList2;
  va_copy(ArgList2, ArgList1);

  int NumChars = std::vsnprintf(nullptr, 0, Format, ArgList1);

  ovk::array<char> MessageChars({NumChars+1});

  va_end(ArgList1);

  std::vsnprintf(MessageChars.Data(), MessageChars.Count(), Format, ArgList2);

  va_end(ArgList2);

  std::string Message(MessageChars.Begin(), MessageChars.End());

  ovk::core::DebugExit(HasMPI, File, Line, Message);

}

}

#endif
