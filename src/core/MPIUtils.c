// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "MPIUtils.h"

#include "Global.h"

void PRIVATE(BroadcastAnySource)(void *Data, int Count, MPI_Datatype DataType, bool IsSource,
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

// Wrapper around MPI_Ibarrier (or crappy alternative if MPI_Ibarrier is not supported) 
// representing a global flag that gets set only after all processes call a routine
// (StartSignal) to set it
// Crappy alternative to Ibarrier:
// * StartSignal sends two messages to the root with different tags (note that the second
//   send is synchronous).
// * In CheckSignal, the root receives the first set of sends, and after they're all done
//   it receives the second set (the barrier is then 'done' when the synchronous sends are
//   matched with receives). The other ranks just test on the two send requests.
// This *seems* like it would be faster than having the root send messages back to each rank
// after the first set of receives, but in reality I have no idea. It would also be better
// to use some kind of tree-like communication instead of sending everything to the root,
// but who has time for that?
struct t_signal {
  MPI_Comm comm;
#ifdef HAVE_IBARRIER
  MPI_Request request;
#else
  char send_buffer[2], recv_buffer[1];
  int comm_size;
  int comm_rank;
  MPI_Request requests[2];
  int num_completed;
#endif
};

void PRIVATE(CreateSignal)(t_signal **Signal_, MPI_Comm Comm) {

  *Signal_ = malloc(sizeof(t_signal));
  t_signal *Signal = *Signal_;

#ifdef HAVE_IBARRIER
  Signal->comm = Comm;
  Signal->request = MPI_REQUEST_NULL;
#else
  MPI_Comm_dup(Comm, &Signal->comm);
  MPI_Comm_size(Signal->comm, &Signal->comm_size);
  MPI_Comm_rank(Signal->comm, &Signal->comm_rank);
  Signal->requests[0] = MPI_REQUEST_NULL;
  Signal->requests[1] = MPI_REQUEST_NULL;
  Signal->num_completed = 0;
#endif

}

void PRIVATE(StartSignal)(t_signal *Signal) {

#ifdef HAVE_IBARRIER
  MPI_Ibarrier(Signal->comm, &Signal->request);
#else
  if (Signal->comm_rank > 0) {
    MPI_Isend(Signal->send_buffer, 1, MPI_CHAR, 0, 0, Signal->comm, Signal->requests);
    MPI_Issend(Signal->send_buffer+1, 1, MPI_CHAR, 0, 1, Signal->comm, Signal->requests+1);
  }
#endif

}

void PRIVATE(CheckSignal)(t_signal *Signal, bool *Done_) {

  int Done;

#ifdef HAVE_IBARRIER
  MPI_Test(&Signal->request, &Done, MPI_STATUS_IGNORE);
#else
  if (Signal->comm_rank > 0) {
    MPI_Testall(2, Signal->requests, &Done, MPI_STATUSES_IGNORE);
  } else {
    Done = false;
    while (Signal->num_completed < Signal->comm_size-1) {
      int IncomingMessage;
      MPI_Status Status;
      MPI_Iprobe(MPI_ANY_SOURCE, 0, Signal->comm, &IncomingMessage, &Status);
      if (!IncomingMessage) break;
      MPI_Recv(Signal->recv_buffer, 1, MPI_CHAR, Status.MPI_SOURCE, 0, Signal->comm,
        MPI_STATUS_IGNORE);
      ++Signal->num_completed;
    }
    if (Signal->num_completed == Signal->comm_size-1) {
      int OtherRank;
      for (OtherRank = 1; OtherRank < Signal->comm_size; ++OtherRank) {
        MPI_Recv(Signal->recv_buffer, 1, MPI_CHAR, OtherRank, 1, Signal->comm, MPI_STATUS_IGNORE);
      }
      Done = true;
    }
  }
  if (Done) {
    Signal->requests[0] = MPI_REQUEST_NULL;
    Signal->requests[1] = MPI_REQUEST_NULL;
    Signal->num_completed = 0;
  }
#endif

  *Done_ = Done;

}

void PRIVATE(DestroySignal)(t_signal **Signal) {

#ifdef HAVE_IBARRIER
  // Do nothing
#else
  MPI_Comm_free(&(*Signal)->comm);
#endif

  free_null(Signal);

}
