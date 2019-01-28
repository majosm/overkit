// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PROFILER_HPP_INCLUDED
#define OVK_CORE_PROFILER_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>

#include <mpi.h>

#include <string>
#include <utility>
#include <vector>

namespace ovk {
namespace core {

namespace profiler_internal {

struct timer {
  bool Active_;
  double LastStartTime_;
  double LastEndTime_;
  double AccumulatedTime_;
};

inline void ResetTimer(timer &Timer);
inline void StartTimer(timer &Timer);
inline void StopTimer(timer &Timer);
inline double GetElapsedTime(const timer &Timer);
inline double GetAccumulatedTime(const timer &Timer);

}

struct profiler {

  using timer = profiler_internal::timer;
  using timer_entry = std::pair<std::string, timer>;

  core::comm Comm_;
  bool Enabled_;
  std::vector<timer_entry> Timers_;

};

void CreateProfiler(profiler &Profiler, const comm &Comm);
void DestroyProfiler(profiler &Profiler);

void EnableProfiler(profiler &Profiler);
void DisableProfiler(profiler &Profiler);

bool TimerExists(const profiler &Profiler, const std::string &TimerName);
int AddProfilerTimer(profiler &Profiler, const std::string &TimerName);
void RemoveProfilerTimer(profiler &Profiler, const std::string &TimerName);

int GetProfilerTimerID(const profiler &Profiler, const std::string &TimerName);

inline void StartProfile(profiler &Profiler, int TimerID);
inline void StartProfileSync(profiler &Profiler, int TimerID, const comm &Comm);
inline void EndProfile(profiler &Profiler, int TimerID);

std::string WriteProfileTimes(const profiler &Profiler);

}}

#include <ovk/core/Profiler.inl>

#endif
