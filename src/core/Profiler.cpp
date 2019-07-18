// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Profiler.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <algorithm>
#include <iterator>
#include <string>
#include <utility>

namespace ovk {
namespace core {

const id_map<1,std::string,false> profiler::TimerNames_ = {
  {HALO_TIME, "Halo"},
  {HALO_SETUP_TIME, "Halo::Setup"},
  {HALO_EXCHANGE_TIME, "Halo::Exchange"},
  {HALO_EXCHANGE_PACK_TIME, "Halo::Exchange::Pack"},
  {HALO_EXCHANGE_MPI_TIME, "Halo::Exchange::MPI"},
  {HALO_EXCHANGE_UNPACK_TIME, "Halo::Exchange::Unpack"},
  {EXCHANGER_COLLECT_TIME, "Exchanger::Collect"},
  {EXCHANGER_COLLECT_PACK_TIME, "Exchanger::Collect::Pack"},
  {EXCHANGER_COLLECT_MPI_TIME, "Exchanger::Collect::MPI"},
  {EXCHANGER_COLLECT_REDUCE_TIME, "Exchanger::Collect::Reduce"},
  {EXCHANGER_SEND_RECV_TIME, "Exchanger::SendRecv"},
  {EXCHANGER_SEND_RECV_PACK_TIME, "Exchanger::SendRecv::Pack"},
  {EXCHANGER_SEND_RECV_MPI_TIME, "Exchanger::SendRecv::MPI"},
  {EXCHANGER_SEND_RECV_UNPACK_TIME, "Exchanger::SendRecv::Unpack"},
  {EXCHANGER_DISPERSE_TIME, "Exchanger::Disperse"},
  {XINTOUT_IMPORT_TIME, "XINTOUT::Import"},
  {XINTOUT_IMPORT_READ_TIME, "XINTOUT::Import::Read"},
  {XINTOUT_IMPORT_READ_MPI_IO_OPEN_TIME, "XINTOUT::Import::Read::MPI-IO::Open"},
  {XINTOUT_IMPORT_READ_MPI_IO_CLOSE_TIME, "XINTOUT::Import::Read::MPI-IO::Close"},
  {XINTOUT_IMPORT_READ_MPI_IO_READ_TIME, "XINTOUT::Import::Read::MPI-IO::Read"},
  {XINTOUT_IMPORT_READ_MPI_IO_OTHER_TIME, "XINTOUT::Import::Read::MPI-IO::Other"},
  {XINTOUT_IMPORT_MATCH_TIME, "XINTOUT::Import::Match"},
  {XINTOUT_IMPORT_MATCH_MAP_TO_BINS_TIME, "XINTOUT::Import::Match::MapToBins"},
  {XINTOUT_IMPORT_MATCH_HANDSHAKE_TIME, "XINTOUT::Import::Match::Handshake"},
  {XINTOUT_IMPORT_MATCH_SEND_TO_BINS_TIME, "XINTOUT::Import::Match::SendToBins"},
  {XINTOUT_IMPORT_MATCH_FILL_CONNECTION_DATA_TIME, "XINTOUT::Import::Match::FillConnectionData"},
  {XINTOUT_IMPORT_MATCH_RECV_FROM_BINS_TIME, "XINTOUT::Import::Match::RecvFromBins"},
  {XINTOUT_IMPORT_MATCH_UNPACK_TIME, "XINTOUT::Import::Match::Unpack"},
  {XINTOUT_IMPORT_DISTRIBUTE_TIME, "XINTOUT::Import::Distribute"},
  {XINTOUT_IMPORT_DISTRIBUTE_MAP_TO_BINS_TIME, "XINTOUT::Import::Distribute::MapToBins"},
  {XINTOUT_IMPORT_DISTRIBUTE_RETRIEVE_BINS_TIME, "XINTOUT::Import::Distribute::RetrieveBins"},
  {XINTOUT_IMPORT_DISTRIBUTE_FIND_RANKS_TIME, "XINTOUT::Import::Distribute::FindRanks"},
  {XINTOUT_IMPORT_DISTRIBUTE_HANDSHAKE_TIME, "XINTOUT::Import::Distribute::Handshake"},
  {XINTOUT_IMPORT_DISTRIBUTE_SEND_DATA_TIME, "XINTOUT::Import::Distribute::SendData"},
  {XINTOUT_IMPORT_SET_CONNECTIVITIES_TIME, "XINTOUT::Import::SetConnectivities"}
};

profiler::profiler(comm_view Comm):
  Comm_(Comm)
{
  OVK_DEBUG_ASSERT(TimerNames_.Count() == profiler_internal_TIMER_ID_COUNT, "Timer name map has "
    "incorrect size.");
}

void profiler::Enable() {

  Enabled_ = true;

}

void profiler::Disable() {

  Enabled_ = false;

}

void profiler::Start_(int TimerID) {

  OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < profiler_internal_TIMER_ID_COUNT, "Invalid timer ID.");

  timer_entry &Entry = Timers_.Fetch(TimerID);
  Entry.Timer.Start();
  ++Entry.ActiveCount;

}

void profiler::StartSync_(int TimerID, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < profiler_internal_TIMER_ID_COUNT, "Invalid timer ID.");

  MPI_Barrier(Comm);

  timer_entry &Entry = Timers_.Fetch(TimerID);
  Entry.Timer.Start();
  ++Entry.ActiveCount;

}

void profiler::Stop_(int TimerID) {

  OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < profiler_internal_TIMER_ID_COUNT, "Invalid timer ID.");

  OVK_DEBUG_ASSERT(Timers_.Contains(TimerID), "Timer %s is not active.", TimerNames_(TimerID));

  timer_entry &Entry = Timers_(TimerID);

  OVK_DEBUG_ASSERT(Entry.ActiveCount > 0, "Timer %s is not active.", TimerNames_(TimerID));

  --Entry.ActiveCount;

  if (Entry.ActiveCount == 0) {
    Entry.Timer.Stop();
  }

}

std::string profiler::WriteProfile() const {

  std::string ProfileString;

  if (Enabled_) {

    id_set<1> GlobalTimerIDs;

    array<int> TimerWasUsed({profiler_internal_TIMER_ID_COUNT}, 0);
    for (int TimerID : Timers_.Keys()) {
      TimerWasUsed(TimerID) = 1;
    }
    MPI_Allreduce(MPI_IN_PLACE, TimerWasUsed.Data(), TimerWasUsed.Count(), MPI_INT, MPI_LOR, Comm_);
    for (int TimerID = 0; TimerID < profiler_internal_TIMER_ID_COUNT; ++TimerID) {
      if (TimerWasUsed(TimerID)) {
        GlobalTimerIDs.Insert(TimerID);
      }
    }

    int NumGlobalTimers = GlobalTimerIDs.Count();

    array<double> Times;
    Times.Reserve(NumGlobalTimers);
    for (int iTimer = 0; iTimer < NumGlobalTimers; ++iTimer) {
      int TimerID = GlobalTimerIDs[iTimer];
      auto Iter = Timers_.Find(TimerID);
      if (Iter != Timers_.End()) {
        const timer &Timer = Iter->Value().Timer;
        Times.Append(Timer.Accumulated());
      } else {
        Times.Append(0.);
      }
    }

    array<double> MinTimes({NumGlobalTimers});
    array<double> MaxTimes({NumGlobalTimers});
    array<double> AvgTimes({NumGlobalTimers});

    MPI_Allreduce(Times.Data(), MinTimes.Data(), NumGlobalTimers, MPI_DOUBLE, MPI_MIN, Comm_);
    MPI_Allreduce(Times.Data(), MaxTimes.Data(), NumGlobalTimers, MPI_DOUBLE, MPI_MAX, Comm_);
    MPI_Allreduce(Times.Data(), AvgTimes.Data(), NumGlobalTimers, MPI_DOUBLE, MPI_SUM, Comm_);

    for (int iTimer = 0; iTimer < NumGlobalTimers; ++iTimer) {
      AvgTimes(iTimer) /= double(Comm_.Size());
    }

    for (int iTimer = 0; iTimer < NumGlobalTimers; ++iTimer) {
      int TimerID = GlobalTimerIDs[iTimer];
      const std::string &TimerName = TimerNames_(TimerID);
      ProfileString += StringPrint("%s: %f %f %f\n", TimerName, MinTimes(iTimer), MaxTimes(iTimer),
        AvgTimes(iTimer));
    }

  }

  return ProfileString;

}

}}
