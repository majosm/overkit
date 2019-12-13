// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DEBUG_HPP_INCLUDED
#define OVK_CORE_DEBUG_HPP_INCLUDED

#include <ovk/core/Debug.h>

#include <ovk/core/Comm.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/TextProcessing.hpp>

#include <mpi.h>

#include <string>
#include <utility>

namespace ovk {
namespace core {

extern int &DebugFlag;

#if OVK_DEBUG

template <typename... Ts> void DebugExit(bool HasMPI, const char *File, int Line, const std::string
  &Format, const Ts &... Args) {

  std::string Prefix = StringPrint("DEBUG ERROR (line %i of file '%s'): ", Line, File);

  std::string Message = StringPrint(Format, Args...);

  if (HasMPI) {
    int Rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    Message = StringReplace(Message, "@rank@", FormatNumber(Rank));
  }

  std::fprintf(stderr, "%s%s\n", Prefix.c_str(), Message.c_str());
  std::fflush(stderr);

  if (HasMPI) {
    MPI_Abort(MPI_COMM_WORLD, 1);
  } else {
    exit(1);
  }

}

#define OVK_DEBUG_ASSERT_CPP(Condition, ...) \
  if (!(Condition)) ovk::core::DebugExit(true, __FILE__, __LINE__, __VA_ARGS__)

#define OVK_DEBUG_ASSERT_NO_MPI_CPP(Condition, ...) \
  if (!(Condition)) ovk::core::DebugExit(false, __FILE__, __LINE__, __VA_ARGS__)

#undef OVK_DEBUG_ASSERT
#define OVK_DEBUG_ASSERT OVK_DEBUG_ASSERT_CPP

#undef OVK_DEBUG_ASSERT_NO_MPI
#define OVK_DEBUG_ASSERT_NO_MPI OVK_DEBUG_ASSERT_NO_MPI_CPP

#else

#define OVK_DEBUG_ASSERT_CPP(...)
#define OVK_DEBUG_ASSERT_NO_MPI_CPP(...)

#endif

}}

#endif
