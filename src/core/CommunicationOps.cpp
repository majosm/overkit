// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/CommunicationOps.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Set.hpp"

#include <mpi.h>

#include <string>

namespace ovk {
namespace core {

void BroadcastString(std::string &String, int Root, comm_view Comm) {

  bool IsRoot = Comm.Rank() == Root;

  int StringLength;
  if (IsRoot) StringLength = String.length();
  MPI_Bcast(&StringLength, 1, MPI_INT, Root, Comm);

  array<char> StringChars({StringLength});
  if (IsRoot) StringChars.Fill(String.begin());

  MPI_Bcast(StringChars.Data(), StringLength, MPI_CHAR, Root, Comm);

  String.assign(StringChars.Begin(), StringChars.End());

}

void BroadcastAnySource(void *Data, int Count, MPI_Datatype DataType, bool IsSource, comm_view Comm)
  {

  // Send to rank 0 first if necessary
  if (Comm.Rank() == 0 && !IsSource) {
    MPI_Recv(Data, Count, DataType, MPI_ANY_SOURCE, 0, Comm, MPI_STATUS_IGNORE);
  } else if (Comm.Rank() != 0 && IsSource) {
    MPI_Send(Data, Count, DataType, 0, 0, Comm);
  }

  MPI_Bcast(Data, Count, DataType, 0, Comm);

}

void BroadcastStringAnySource(std::string &String, bool IsSource, comm_view Comm) {

  // Send to rank 0 first if necessary
  if (Comm.Rank() == 0 && !IsSource) {
    int StringLength;
    MPI_Recv(&StringLength, 1, MPI_INT, MPI_ANY_SOURCE, 0, Comm, MPI_STATUS_IGNORE);
    array<char> StringChars({StringLength});
    MPI_Recv(StringChars.Data(), StringLength, MPI_CHAR, MPI_ANY_SOURCE, 0, Comm,
      MPI_STATUS_IGNORE);
    String.assign(StringChars.Begin(), StringChars.End());
  } else if (Comm.Rank() != 0 && IsSource) {
    int StringLength = String.length();
    MPI_Send(&StringLength, 1, MPI_INT, 0, 0, Comm);
    array<char> StringChars({StringLength}, String.begin());
    MPI_Send(StringChars.Data(), StringLength, MPI_CHAR, 0, 0, Comm);
  }

  BroadcastString(String, 0, Comm);

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

signal::signal(comm_view Comm):
  Comm_(DuplicateComm(Comm))
{}

void signal::Start() {

#ifdef OVK_HAVE_MPI_IBARRIER
  MPI_Ibarrier(Comm_, &Request_);
#else
  if (Comm_.Rank() > 0) {
    MPI_Isend(SendBuffer_, 1, MPI_UNSIGNED_CHAR, 0, 0, Comm_, Requests_);
    MPI_Issend(SendBuffer_+1, 1, MPI_UNSIGNED_CHAR, 0, 1, Comm_, Requests_+1);
  }
#endif

}

bool signal::Check() {

  bool Done;

#ifdef OVK_HAVE_MPI_IBARRIER
  int DoneInt;
  MPI_Test(&Request_, &DoneInt, MPI_STATUS_IGNORE);
  Done = bool(DoneInt);
#else
  if (Comm_.Rank() > 0) {
    int DoneInt;
    MPI_Testall(2, Requests_, &DoneInt, MPI_STATUSES_IGNORE);
    Done = bool(DoneInt);
  } else {
    Done = false;
    while (NumCompleted_ < Comm_.Size()-1) {
      int IncomingMessage;
      MPI_Status Status;
      MPI_Iprobe(MPI_ANY_SOURCE, 0, Comm_, &IncomingMessage, &Status);
      if (!IncomingMessage) break;
      MPI_Recv(RecvBuffer_, 1, MPI_UNSIGNED_CHAR, Status.MPI_SOURCE, 0, Comm_, MPI_STATUS_IGNORE);
      ++NumCompleted_;
    }
    if (NumCompleted_ == Comm_.Size()-1) {
      int OtherRank;
      for (OtherRank = 1; OtherRank < Comm_.Size(); ++OtherRank) {
        MPI_Recv(RecvBuffer_, 1, MPI_UNSIGNED_CHAR, OtherRank, 1, Comm_, MPI_STATUS_IGNORE);
      }
      Done = true;
    }
  }
  if (Done) {
    Requests_[0] = MPI_REQUEST_NULL;
    Requests_[1] = MPI_REQUEST_NULL;
    NumCompleted_ = 0;
  }
#endif

  return Done;

}

array<int> DynamicHandshake(comm_view Comm_, array_view<const int> Ranks) {

  // Duplicate comm to avoid matching with any sends/recvs before and after
  comm Comm = DuplicateComm(Comm_);

  byte SendData = 0;
  byte RecvData;

  int Rank = Comm.Rank();

  set<int> MatchedRanksSet;

  for (int OtherRank : Ranks) {
    if (OtherRank == Rank) {
      MatchedRanksSet.Insert(Rank);
    }
  }

  array<MPI_Request> SendRequests;
  SendRequests.Reserve(Ranks.Count());
  for (int OtherRank : Ranks) {
    if (OtherRank != Rank) {
      MPI_Issend(&SendData, 1, MPI_UNSIGNED_CHAR, OtherRank, 0, Comm, &SendRequests.Append());
    }
  }

  signal AllSendsDoneSignal(Comm);

  array<MPI_Request> RecvRequests;

  bool Done = false;
  int SendsDone = false;
  while (!Done) {
    while (true) {
      int IncomingMessage;
      MPI_Status Status;
      MPI_Iprobe(MPI_ANY_SOURCE, 0, Comm, &IncomingMessage, &Status);
      if (!IncomingMessage) break;
      int MatchedRank = Status.MPI_SOURCE;
      MPI_Irecv(&RecvData, 1, MPI_UNSIGNED_CHAR, MatchedRank, 0, Comm, &RecvRequests.Append());
      MatchedRanksSet.Insert(MatchedRank);
    }
    if (SendsDone) {
      Done = AllSendsDoneSignal.Check();
    } else {
      MPI_Testall(SendRequests.Count(), SendRequests.Data(), &SendsDone, MPI_STATUSES_IGNORE);
      if (SendsDone) {
        AllSendsDoneSignal.Start();
      }
    }
  }

  MPI_Waitall(RecvRequests.Count(), RecvRequests.Data(), MPI_STATUSES_IGNORE);

  return {MatchedRanksSet};

}

hang_detector::hang_detector(comm_view Comm, double Timeout):
  Comm_(DuplicateComm(Comm)),
  Signal_(Comm_),
  Timeout_(Timeout)
{}

void hang_detector::operator()() {

  timer Timer;

  Signal_.Start();
  Timer.Start();

  bool Done = false;
  while (!Done) {
    Done = Signal_.Check();
    if (Timer.Elapsed() > Timeout_) break;
  }

  if (!Done) {
    // Can't actually serialize printing, so we have to fake it
    int Rank, NumProcs;
    MPI_Comm_rank(Comm_, &Rank);
    MPI_Comm_size(Comm_, &NumProcs);
    Timer.Stop();
    Timer.Start();
    while (Timer.Elapsed() < 0.01*Rank) {}
    std::fprintf(stderr, "Hang detected: process %i reporting in.\n", Rank);
    std::fflush(stderr);
    while (Timer.Elapsed() < 0.01*NumProcs+1.) {}
    MPI_Abort(Comm_, 1);
  }

}

}}
