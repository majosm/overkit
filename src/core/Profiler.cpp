// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
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

void CreateProfiler(profiler &Profiler, const comm &Comm) {

  Profiler.Comm_ = DuplicateComm(Comm);

  Profiler.Enabled_ = false;

}

void DestroyProfiler(profiler &Profiler) {

  Profiler.Timers_.Clear();
  Profiler.Comm_.Reset();

}

void EnableProfiler(profiler &Profiler) {

  Profiler.Enabled_ = true;

}

void DisableProfiler(profiler &Profiler) {

  Profiler.Enabled_ = false;

}

namespace {

inline array<profiler::timer_entry>::iterator FindTimer(profiler &Profiler, const std::string
  &TimerName) {

  auto NameMatches = [&TimerName](const profiler::timer_entry &TimerEntry) -> bool {
    return TimerEntry.first == TimerName;
  };

  return std::find_if(Profiler.Timers_.LinearBegin(), Profiler.Timers_.LinearEnd(), NameMatches);

}

inline array<profiler::timer_entry>::const_iterator FindTimer(const profiler &Profiler, const
  std::string &TimerName) {

  auto NameMatches = [&TimerName](const profiler::timer_entry &TimerEntry) -> bool {
    return TimerEntry.first == TimerName;
  };

  return std::find_if(Profiler.Timers_.LinearBegin(), Profiler.Timers_.LinearEnd(), NameMatches);

}

}

bool TimerExists(const profiler &Profiler, const char *TimerName) {

  if (Profiler.Enabled_) {
    return FindTimer(Profiler, std::string(TimerName)) != Profiler.Timers_.LinearEnd();
  } else {
    return false;
  }

}

int AddProfilerTimer(profiler &Profiler, const char *TimerName) {

  if (Profiler.Enabled_) {

    std::string TimerNameString(TimerName);

    auto Iter = FindTimer(Profiler, TimerNameString);

    if (Iter == Profiler.Timers_.LinearEnd()) {

      profiler_internal::timer Timer;
      ResetTimer(Timer);

      Profiler.Timers_.Append(TimerNameString, Timer);

      return Profiler.Timers_.Count()-1;

    } else {

      return std::distance(Profiler.Timers_.LinearBegin(), Iter);

    }

  } else {

    return -1;

  }

}

void RemoveProfilerTimer(profiler &Profiler, const char *TimerName) {

  if (Profiler.Enabled_) {

    auto Iter = FindTimer(Profiler, std::string(TimerName));

    if (Iter != Profiler.Timers_.LinearEnd()) {

      int Index = std::distance(Profiler.Timers_.LinearBegin(), Iter);

      int NumTimers = Profiler.Timers_.Count();

      array<profiler::timer_entry> NewTimers({NumTimers-1});

      for (int iTimer = 0; iTimer < Index; ++iTimer) {
        NewTimers(iTimer) = std::move(Profiler.Timers_(iTimer));
      }
      for (int iTimer = Index+1; iTimer < NumTimers; ++iTimer) {
        NewTimers(iTimer-1) = std::move(Profiler.Timers_(iTimer));
      }

      Profiler.Timers_ = std::move(NewTimers);

    }

  }

}

int GetProfilerTimerID(const profiler &Profiler, const char *TimerName) {

  if (Profiler.Enabled_) {

    std::string TimerNameString(TimerName);

    auto Iter = FindTimer(Profiler, TimerNameString);

    OVK_DEBUG_ASSERT(Iter != Profiler.Timers_.end(), "Timer '%s' does not exist.", TimerNameString);
    OVK_DEBUG_ASSERT(Iter != Profiler.Timers_.LinearEnd(), "Timer '%s' does not exist.",
      TimerNameString);

    return std::distance(Profiler.Timers_.LinearBegin(), Iter);

  } else {

    return -1;

  }

}

std::string WriteProfileTimes(const profiler &Profiler) {

  std::string ProfileTimesString;

  if (Profiler.Enabled_) {

    int NumTimers = Profiler.Timers_.Count();

    array<double> Times({NumTimers});
    for (int iTimer = 0; iTimer < int(NumTimers); ++iTimer) {
      const profiler::timer &Timer = Profiler.Timers_(iTimer).second;
      Times(iTimer) = GetAccumulatedTime(Timer);
    }

    array<double> MinTimes({NumTimers});
    array<double> MaxTimes({NumTimers});
    array<double> AvgTimes({NumTimers});

    MPI_Allreduce(Times.Data(), MinTimes.Data(), NumTimers, MPI_DOUBLE, MPI_MIN, Profiler.Comm_);
    MPI_Allreduce(Times.Data(), MaxTimes.Data(), NumTimers, MPI_DOUBLE, MPI_MAX, Profiler.Comm_);
    MPI_Allreduce(Times.Data(), AvgTimes.Data(), NumTimers, MPI_DOUBLE, MPI_SUM, Profiler.Comm_);

    for (int iTimer = 0; iTimer < NumTimers; ++iTimer) {
      AvgTimes(iTimer) /= double(Profiler.Comm_.Size());
    }

    for (int iTimer = 0; iTimer < NumTimers; ++iTimer) {
      const std::string &TimerName = Profiler.Timers_(iTimer).first;
      ProfileTimesString += StringPrint("%s: %f %f %f\n", TimerName, MinTimes(iTimer),
        MaxTimes(iTimer), AvgTimes(iTimer));
    }

  }

  return ProfileTimesString;

}

}}
