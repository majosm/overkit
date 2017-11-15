// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PARALLEL_UTILS_INCLUDED
#define OVK_CORE_PARALLEL_UTILS_INCLUDED

#include "Global.h"

static inline void BroadcastAnySource(void *Data, int Count, MPI_Datatype DataType, bool IsSource,
  MPI_Comm Comm) {

  int Rank;
  MPI_Comm_rank(Comm, &Rank);

  // Send to the root rank first if necessary
  if (Rank == 0 && !IsSource) {
    MPI_Recv(Data, Count, DataType, MPI_ANY_SOURCE, 0, Comm, MPI_STATUS_IGNORE);
  } else if (Rank != 0 && IsSource) {
    MPI_Send(Data, Count, DataType, 0, 0, Comm);
  }

  MPI_Bcast(Data, Count, DataType, 0, Comm);

}

#endif
