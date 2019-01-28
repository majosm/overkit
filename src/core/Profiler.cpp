// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Profiler.hpp"

#include "ovk/core/Comm.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/TextProcessing.hpp"

#include <mpi.h>

#include <algorithm>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

namespace ovk {
namespace core {

void CreateProfiler(profiler &Profiler, const comm &Comm) {

  Profiler.Comm_ = DuplicateComm(Comm);

  Profiler.Enabled_ = false;

}

void DestroyProfiler(profiler &Profiler) {

  Profiler.Timers_.clear();
  Profiler.Comm_.Reset();

}

void EnableProfiler(profiler &Profiler) {

  Profiler.Enabled_ = true;

}

void DisableProfiler(profiler &Profiler) {

  Profiler.Enabled_ = false;

}

namespace {

inline std::vector<profiler::timer_entry>::iterator FindTimer(profiler &Profiler, const std::string
  &TimerName) {

  auto NameMatches = [&TimerName](const profiler::timer_entry &TimerEntry) -> bool {
    return TimerEntry.first == TimerName;
  };

  return std::find_if(Profiler.Timers_.begin(), Profiler.Timers_.end(), NameMatches);

}

inline std::vector<profiler::timer_entry>::const_iterator FindTimer(const profiler &Profiler, const
  std::string &TimerName) {

  auto NameMatches = [&TimerName](const profiler::timer_entry &TimerEntry) -> bool {
    return TimerEntry.first == TimerName;
  };

  return std::find_if(Profiler.Timers_.begin(), Profiler.Timers_.end(), NameMatches);

}

}

bool TimerExists(const profiler &Profiler, const std::string &TimerName) {

  return FindTimer(Profiler, TimerName) != Profiler.Timers_.end();

}

int AddProfilerTimer(profiler &Profiler, const std::string &TimerName) {

  auto Iter = FindTimer(Profiler, TimerName);

  if (Iter == Profiler.Timers_.end()) {

    profiler_internal::timer Timer;
    ResetTimer(Timer);

    Profiler.Timers_.emplace_back(TimerName, Timer);

    return Profiler.Timers_.size()-1;

  } else {

    return std::distance(Profiler.Timers_.begin(), Iter);

  }

}

void RemoveProfilerTimer(profiler &Profiler, const std::string &TimerName) {

  auto Iter = FindTimer(Profiler, TimerName);

  OVK_DEBUG_ASSERT(Iter != Profiler.Timers_.end(), "Timer does not exist.");

  Profiler.Timers_.erase(Iter);

}

int GetProfilerTimerID(const profiler &Profiler, const std::string &TimerName) {

  auto Iter = FindTimer(Profiler, TimerName);

  OVK_DEBUG_ASSERT(Iter != Profiler.Timers_.end(), "Timer does not exist.");

  return std::distance(Profiler.Timers_.begin(), Iter);

}

std::string WriteProfileTimes(const profiler &Profiler) {

  std::string ProfileTimesString;

  if (Profiler.Enabled_) {

    std::vector<double> Times(Profiler.Timers_.size());
    for (int iTimer = 0; iTimer < int(Profiler.Timers_.size()); ++iTimer) {
      const profiler::timer &Timer = Profiler.Timers_[iTimer].second;
      Times[iTimer] = GetAccumulatedTime(Timer);
    }

    std::vector<double> MinTimes(Profiler.Timers_.size());
    std::vector<double> MaxTimes(Profiler.Timers_.size());
    std::vector<double> AvgTimes(Profiler.Timers_.size());

    MPI_Allreduce(Times.data(), MinTimes.data(), Profiler.Timers_.size(), MPI_DOUBLE, MPI_MIN,
      Profiler.Comm_);
    MPI_Allreduce(Times.data(), MaxTimes.data(), Profiler.Timers_.size(), MPI_DOUBLE, MPI_MAX,
      Profiler.Comm_);
    MPI_Allreduce(Times.data(), AvgTimes.data(), Profiler.Timers_.size(), MPI_DOUBLE, MPI_SUM,
      Profiler.Comm_);

    for (int iTimer = 0; iTimer < int(Profiler.Timers_.size()); ++iTimer) {
      AvgTimes[iTimer] /= double(Profiler.Comm_.Size());
    }

    for (int iTimer = 0; iTimer < int(Profiler.Timers_.size()); ++iTimer) {
      const std::string &TimerName = Profiler.Timers_[iTimer].first;
      ProfileTimesString += StringPrint("%s: %f %f %f\n", TimerName, MinTimes[iTimer],
        MaxTimes[iTimer], AvgTimes[iTimer]);
    }

  }

  return ProfileTimesString;

}

}}
