// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Misc.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <set>

namespace ovk {
namespace core {

void BroadcastAnySource(void *Data, int Count, MPI_Datatype DataType, bool IsSource, MPI_Comm Comm)
  {

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

// Notes about crappy alternative to Ibarrier:
// * StartSignal sends two messages to the root with different tags (note that the second
//   send is synchronous).
// * In CheckSignal, the root receives the first set of sends, and after they're all done
//   it receives the second set (the barrier is then 'done' when the synchronous sends are
//   matched with receives). The other ranks just test on the two send requests.
// This *seems* like it would be faster than having the root send messages back to each rank
// after the first set of receives, but in reality I have no idea. It would also be better
// to use some kind of tree-like communication instead of sending everything to the root,
// but who has time for that?

void CreateSignal(signal &Signal, MPI_Comm Comm) {

#ifdef OVK_HAVE_MPI_IBARRIER
  Signal.Comm_ = Comm;
  Signal.Request_ = MPI_REQUEST_NULL;
#else
  MPI_Comm_dup(Comm, &Signal.Comm_);
  MPI_Comm_size(Signal.Comm_, &Signal.CommSize_);
  MPI_Comm_rank(Signal.Comm_, &Signal.CommRank_);
  Signal.Requests_[0] = MPI_REQUEST_NULL;
  Signal.Requests_[1] = MPI_REQUEST_NULL;
  Signal.NumCompleted_ = 0;
#endif

}

void StartSignal(signal &Signal) {

#ifdef OVK_HAVE_MPI_IBARRIER
  MPI_Ibarrier(Signal.Comm_, &Signal.Request_);
#else
  if (Signal.CommRank_ > 0) {
    MPI_Isend(Signal.SendBuffer_, 1, MPI_CHAR, 0, 0, Signal.Comm_, Signal.Requests_);
    MPI_Issend(Signal.SendBuffer_+1, 1, MPI_CHAR, 0, 1, Signal.Comm_, Signal.Requests_+1);
  }
#endif

}

void CheckSignal(signal &Signal, bool &Done) {

  int DoneInt;

#ifdef OVK_HAVE_MPI_IBARRIER
  MPI_Test(&Signal.Request_, &DoneInt, MPI_STATUS_IGNORE);
  Done = bool(DoneInt);
#else
  if (Signal.CommRank_ > 0) {
    MPI_Testall(2, Signal.Requests_, &DoneInt, MPI_STATUSES_IGNORE);
    Done = bool(DoneInt);
  } else {
    Done = false;
    while (Signal.NumCompleted_ < Signal.CommSize_-1) {
      int IncomingMessage;
      MPI_Status Status;
      MPI_Iprobe(MPI_ANY_SOURCE, 0, Signal.Comm_, &IncomingMessage, &Status);
      if (!IncomingMessage) break;
      MPI_Recv(Signal.RecvBuffer_, 1, MPI_CHAR, Status.MPI_SOURCE, 0, Signal.Comm_,
        MPI_STATUS_IGNORE);
      ++Signal.NumCompleted_;
    }
    if (Signal.NumCompleted_ == Signal.CommSize_-1) {
      int OtherRank;
      for (OtherRank = 1; OtherRank < Signal.CommSize_; ++OtherRank) {
        MPI_Recv(Signal.RecvBuffer_, 1, MPI_CHAR, OtherRank, 1, Signal.Comm_, MPI_STATUS_IGNORE);
      }
      Done = true;
    }
  }
  if (Done) {
    Signal.Requests_[0] = MPI_REQUEST_NULL;
    Signal.Requests_[1] = MPI_REQUEST_NULL;
    Signal.NumCompleted_ = 0;
  }
#endif

}

void DestroySignal(signal &Signal) {

#ifdef OVK_HAVE_MPI_IBARRIER
  // Do nothing
#else
  MPI_Comm_free(&Signal.Comm_);
#endif

}

void DynamicHandshake(MPI_Comm Comm, array_view<const int> DestRanks, array<int> &SourceRanks) {

  char SendBuffer[1], RecvBuffer[1];

  int NumDestRanks = DestRanks.Count();

  array<MPI_Request> SendRequests;
  SendRequests.Reserve(NumDestRanks);
  for (int iDestRank = 0; iDestRank < NumDestRanks; ++iDestRank) {
    MPI_Request &Request = SendRequests.Append();
    MPI_Issend(SendBuffer, 1, MPI_CHAR, DestRanks(iDestRank), 0, Comm, &Request);
  }

  signal AllSendsDoneSignal;
  CreateSignal(AllSendsDoneSignal, Comm);

  std::set<int> SourceRanksSet;

  bool Done = false;
  int SendsDone = false;
  while (!Done) {
    while (true) {
      int IncomingMessage;
      MPI_Status Status;
      MPI_Iprobe(MPI_ANY_SOURCE, 0, Comm, &IncomingMessage, &Status);
      if (!IncomingMessage) break;
      int SourceRank = Status.MPI_SOURCE;
      MPI_Recv(RecvBuffer, 1, MPI_CHAR, SourceRank, 0, Comm, MPI_STATUS_IGNORE);
      SourceRanksSet.insert(SourceRank);
    }
    if (SendsDone) {
      CheckSignal(AllSendsDoneSignal, Done);
    } else {
      MPI_Testall(NumDestRanks, SendRequests.Data(), &SendsDone, MPI_STATUSES_IGNORE);
      if (SendsDone) {
        StartSignal(AllSendsDoneSignal);
      }
    }
  }

  MPI_Barrier(Comm);

  DestroySignal(AllSendsDoneSignal);

  SourceRanks.Assign({int(SourceRanksSet.size())}, SourceRanksSet.begin());

}

}}
