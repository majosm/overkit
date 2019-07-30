// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "support/Decomp.hpp"

#include <ovk/core/Comm.hpp>

#include <mpi.h>

#ifdef OVK_POSIX_SYSTEM
#include <sys/types.h>
#include <unistd.h>
#endif

#include <cstdio>

namespace support {

#ifdef OVK_POSIX_SYSTEM
void DebuggerAttachHelper() {

  ovk::comm_view Comm(MPI_COMM_WORLD);

  for (int Rank = 0; Rank < Comm.Size(); ++Rank) {
    if (Rank == Comm.Rank()) {
      std::printf("Rank %i has pid %i.\n", Rank, getpid()); std::fflush(stdout);
    }
  }

  sleep(10);

}
#endif

}

#ifdef __cplusplus
extern "C" {
#endif

#ifdef OVK_POSIX_SYSTEM
void support_DebuggerAttachHelper() {
  support::DebuggerAttachHelper();
}
#endif

#ifdef __cplusplus
}
#endif
