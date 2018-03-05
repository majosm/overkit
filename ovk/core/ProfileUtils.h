// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PROFILE_UTILS_INCLUDED
#define OVK_CORE_PROFILE_UTILS_INCLUDED

#include "ovk/core/Global.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  bool active;
  double last_start_time;
  double last_end_time;
  double accumulated_time;
} t_timer;

typedef struct {
  char name[OVK_NAME_LENGTH];
  t_timer timer;
} t_profiler_timer_info;

typedef struct {
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  bool enabled;
  int num_timers;
  t_profiler_timer_info **timers;
} t_profiler;

static inline void ResetTimer(t_timer *Timer) {

  Timer->active = false;
  Timer->last_start_time = 0.;
  Timer->last_end_time = 0.;
  Timer->accumulated_time = 0.;

}

static inline void StartTimer(t_timer *Timer) {

  if (!Timer->active) {
    Timer->last_start_time = MPI_Wtime();
  }
  Timer->active = true;

}

static inline void StopTimer(t_timer *Timer) {

  if (Timer->active) {
    Timer->last_end_time = MPI_Wtime();
    Timer->accumulated_time += Timer->last_end_time - Timer->last_start_time;
  }
  Timer->active = false;

}

static inline double GetElapsedTime(t_timer *Timer) {

  if (Timer->active) {
    return MPI_Wtime() - Timer->last_start_time;
  } else {
    return Timer->last_end_time - Timer->last_start_time;
  }

}

static inline double GetAccumulatedTime(t_timer *Timer) {

  if (Timer->active) {
    return Timer->accumulated_time + (MPI_Wtime() - Timer->last_start_time);
  } else {
    return Timer->accumulated_time;
  }

}

void PRIVATE(CreateProfiler)(t_profiler **Profiler_, MPI_Comm Comm_);
#define CreateProfiler(...) PRIVATE(CreateProfiler)(__VA_ARGS__)

void PRIVATE(DestroyProfiler)(t_profiler **Profiler_);
#define DestroyProfiler(...) PRIVATE(DestroyProfiler)(__VA_ARGS__)

void PRIVATE(EnableProfiler)(t_profiler *Profiler);
#define EnableProfiler(...) PRIVATE(EnableProfiler)(__VA_ARGS__)

void PRIVATE(DisableProfiler)(t_profiler *Profiler);
#define DisableProfiler(...) PRIVATE(DisableProfiler)(__VA_ARGS__)

int PRIVATE(AddProfilerTimer)(t_profiler *Profiler, const char *TimerName);
#define AddProfilerTimer(...) PRIVATE(AddProfilerTimer)(__VA_ARGS__)

void PRIVATE(RemoveProfilerTimer)(t_profiler *Profiler, const char *TimerName);
#define RemoveProfilerTimer(...) PRIVATE(RemoveProfilerTimer)(__VA_ARGS__)

int PRIVATE(GetProfilerTimerID)(t_profiler *Profiler, const char *TimerName);
#define GetProfilerTimerID(...) PRIVATE(GetProfilerTimerID)(__VA_ARGS__)

void PRIVATE(WriteProfileTimes)(const t_profiler *Profiler, FILE *File);
#define WriteProfileTimes(...) PRIVATE(WriteProfileTimes)(__VA_ARGS__)

static inline void StartProfile(t_profiler *Profiler, int TimerID) {

  if (Profiler->enabled) {
    OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < Profiler->num_timers, "Invalid timer ID.");
    t_profiler_timer_info *TimerInfo = Profiler->timers[TimerID];
    StartTimer(&TimerInfo->timer);
  }

}

static inline void EndProfile(t_profiler *Profiler, int TimerID) {

  if (Profiler->enabled) {
    OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < Profiler->num_timers, "Invalid timer ID.");
    t_profiler_timer_info *TimerInfo = Profiler->timers[TimerID];
    StopTimer(&TimerInfo->timer);
  }

}

static inline void StartProfileSync(t_profiler *Profiler, int TimerID, MPI_Comm Comm) {

  if (Profiler->enabled) {
    OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < Profiler->num_timers, "Invalid timer ID.");
    t_profiler_timer_info *TimerInfo = Profiler->timers[TimerID];
    MPI_Barrier(Comm);
    StartTimer(&TimerInfo->timer);
  }

}

static inline void EndProfileSync(t_profiler *Profiler, int TimerID, MPI_Comm Comm) {

  if (Profiler->enabled) {
    OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < Profiler->num_timers, "Invalid timer ID.");
    t_profiler_timer_info *TimerInfo = Profiler->timers[TimerID];
    MPI_Barrier(Comm);
    StopTimer(&TimerInfo->timer);
  }

}

#ifdef __cplusplus
}
#endif

#endif
