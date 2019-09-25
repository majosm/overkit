// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DEBUG_HPP_INCLUDED
#define OVK_CORE_DEBUG_HPP_INCLUDED

#include <ovk/core/Debug.h>

#include <ovk/core/Global.hpp>
#include <ovk/core/TextProcessing.hpp>

#include <mpi.h>

#include <string>

namespace ovk {
namespace core {

extern int &DebugFlag;

#if OVK_DEBUG

template <typename... Ts> void DebugExit(const char *File, int Line, const std::string &Format,
  const Ts &... Args) {

  std::string Prefix = StringPrint("DEBUG ERROR (line %i of file '%s'): ", Line, File);

  int Rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

  std::string Message = StringPrint(Format, Args...);
  Message = StringReplace(Message, "@rank@", FormatNumber(Rank));

  std::fprintf(stderr, "%s%s\n", Prefix.c_str(), Message.c_str());
  std::fflush(stderr);

  MPI_Abort(MPI_COMM_WORLD, 1);

}

#define OVK_DEBUG_ASSERT_CPP(Condition, ...) \
  if (!(Condition)) ovk::core::DebugExit(__FILE__, __LINE__, __VA_ARGS__)

#undef OVK_DEBUG_ASSERT
#define OVK_DEBUG_ASSERT OVK_DEBUG_ASSERT_CPP

#else

#define OVK_DEBUG_ASSERT_CPP(...)

#endif

}}

#endif
