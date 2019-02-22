// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PROFILER_HPP_INCLUDED
#define OVK_CORE_PROFILER_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>

#include <mpi.h>

#include <string>
#include <utility>

namespace ovk {
namespace core {

namespace profiler_internal {

struct timer {
  bool Active_;
  double LastStartTime_;
  double LastEndTime_;
  double AccumulatedTime_;
};

OVK_FORCE_INLINE void ResetTimer(timer &Timer);
OVK_FORCE_INLINE void StartTimer(timer &Timer);
OVK_FORCE_INLINE void StopTimer(timer &Timer);
OVK_FORCE_INLINE double GetElapsedTime(const timer &Timer);
OVK_FORCE_INLINE double GetAccumulatedTime(const timer &Timer);

}

struct profiler {

  using timer = profiler_internal::timer;
  using timer_entry = std::pair<std::string, timer>;

  core::comm Comm_;
  bool Enabled_;
  array<timer_entry> Timers_;

};

void CreateProfiler(profiler &Profiler, const comm &Comm);
void DestroyProfiler(profiler &Profiler);

void EnableProfiler(profiler &Profiler);
void DisableProfiler(profiler &Profiler);

// Use const char * instead of std::string so calls have minimal overhead when profiler is disabled
bool TimerExists(const profiler &Profiler, const char *TimerName);
int AddProfilerTimer(profiler &Profiler, const char *TimerName);
void RemoveProfilerTimer(profiler &Profiler, const char *TimerName);
int GetProfilerTimerID(const profiler &Profiler, const char *TimerName);

OVK_FORCE_INLINE void StartProfile(profiler &Profiler, int TimerID);
OVK_FORCE_INLINE void StartProfileSync(profiler &Profiler, int TimerID, const comm &Comm);
OVK_FORCE_INLINE void EndProfile(profiler &Profiler, int TimerID);

std::string WriteProfileTimes(const profiler &Profiler);

}}

#include <ovk/core/Profiler.inl>

#endif
