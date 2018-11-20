// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DEBUG_HPP_INCLUDED
#define OVK_CORE_DEBUG_HPP_INCLUDED

#include <ovk/core/DebugBase.h>

#include <ovk/core/Global.hpp>
#include <ovk/core/TextProcessing.hpp>

#include <mpi.h>

#include <string>

#if OVK_DEBUG

namespace ovk {
namespace core {

template <typename... Ts> void DebugAssert(const char *File, int Line, const std::string &Format,
  const Ts &... Args) {

  std::string Message;
  Message += StringPrint("DEBUG ERROR (line %i of file '%s'): ", Line, File);
  Message += StringPrint(Format, Args...);

  fprintf(stderr, "%s\n", Message.c_str());
  fflush(stderr);

  MPI_Abort(MPI_COMM_WORLD, 1);

}

#define OVK_DEBUG_ASSERT_CPP(Condition, ...) \
  if (!(Condition)) ovk::core::DebugAssert(__FILE__, __LINE__, __VA_ARGS__)

#undef OVK_DEBUG_ASSERT
#define OVK_DEBUG_ASSERT OVK_DEBUG_ASSERT_CPP

}}

#else

#define OVK_DEBUG_ASSERT_CPP(...)

#endif

#endif
