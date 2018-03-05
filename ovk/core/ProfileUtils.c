// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ProfileUtils.h"

#include "ovk/core/Global.h"

#include <stdio.h>
#include <string.h>

void PRIVATE(CreateProfiler)(t_profiler **Profiler_, MPI_Comm Comm_) {

  *Profiler_ = malloc(sizeof(t_profiler));
  t_profiler *Profiler = *Profiler_;

  MPI_Comm Comm;
  int CommSize, CommRank;
  MPI_Comm_dup(Comm_, &Comm);
  MPI_Comm_size(Comm, &CommSize);
  MPI_Comm_rank(Comm, &CommRank);

  Profiler->comm = Comm;
  Profiler->comm_size = CommSize;
  Profiler->comm_rank = CommRank;

  Profiler->enabled = false;

  Profiler->num_timers = 0;
  Profiler->timers = NULL;

}

void PRIVATE(DestroyProfiler)(t_profiler **Profiler_) {

  int iTimer;

  t_profiler *Profiler = *Profiler_;

  for (iTimer = 0; iTimer < Profiler->num_timers; ++iTimer) {
    t_profiler_timer_info *TimerInfo = Profiler->timers[iTimer];
    free(TimerInfo);
  }

  free(Profiler->timers);

  MPI_Comm_free(&Profiler->comm);

  free_null(Profiler_);

}

void PRIVATE(EnableProfiler)(t_profiler *Profiler) {

  Profiler->enabled = true; 

}

void PRIVATE(DisableProfiler)(t_profiler *Profiler) {

  Profiler->enabled = false; 

}

int PRIVATE(AddProfilerTimer)(t_profiler *Profiler, const char *TimerName) {

  if (OVK_DEBUG) {
    int iTimer;
    for (iTimer = 0; iTimer < Profiler->num_timers; ++iTimer) {
      OVK_DEBUG_ASSERT(strcmp(Profiler->timers[iTimer]->name, TimerName) != 0,
        "Timer already exists.");
    }
  }

  ++Profiler->num_timers;
  Profiler->timers = realloc(Profiler->timers, Profiler->num_timers*sizeof(t_profiler_timer_info *));

  t_profiler_timer_info *TimerInfo = malloc(sizeof(t_profiler_timer_info));
  strncpy(TimerInfo->name, TimerName, OVK_NAME_LENGTH);
  ResetTimer(&TimerInfo->timer);

  Profiler->timers[Profiler->num_timers-1] = TimerInfo;

  return Profiler->num_timers-1;

}

void PRIVATE(RemoveProfilerTimer)(t_profiler *Profiler, const char *TimerName) {

  int iTimer, jTimer;

  for (iTimer = 0; iTimer < Profiler->num_timers; ++iTimer) {
    if (strcmp(Profiler->timers[iTimer]->name, TimerName) == 0) break;
  }

  OVK_DEBUG_ASSERT(iTimer < Profiler->num_timers, "Timer does not exist.");

  t_profiler_timer_info *TimerInfo = Profiler->timers[iTimer];
  free(TimerInfo);

  for (jTimer = iTimer+1; jTimer < Profiler->num_timers; ++jTimer) {
    Profiler->timers[jTimer-1] = Profiler->timers[jTimer];
  }
  
  --Profiler->num_timers;
  Profiler->timers = realloc(Profiler->timers, Profiler->num_timers*sizeof(t_profiler_timer_info *));

}

int PRIVATE(GetProfilerTimerID)(t_profiler *Profiler, const char *TimerName) {

  int iTimer;

  for (iTimer = 0; iTimer < Profiler->num_timers; ++iTimer) {
    if (strcmp(Profiler->timers[iTimer]->name, TimerName) == 0) break;
  }

  OVK_DEBUG_ASSERT(iTimer < Profiler->num_timers, "Timer does not exist.");

  return iTimer;

}

void PRIVATE(WriteProfileTimes)(const t_profiler *Profiler, FILE *File) {

  int iTimer;

  if (Profiler->enabled) {

    double *Times = malloc(Profiler->num_timers*sizeof(double));
    for (iTimer = 0; iTimer < Profiler->num_timers; ++iTimer) {
      t_profiler_timer_info *TimerInfo = Profiler->timers[iTimer];
      Times[iTimer] = GetAccumulatedTime(&TimerInfo->timer);
    }

    double *MinTimes = NULL, *MaxTimes = NULL, *AvgTimes = NULL;
    if (Profiler->comm_rank == 0) {
      MinTimes = malloc(Profiler->num_timers*sizeof(double));
      MaxTimes = malloc(Profiler->num_timers*sizeof(double));
      AvgTimes = malloc(Profiler->num_timers*sizeof(double));
    }
    MPI_Reduce(Times, MinTimes, Profiler->num_timers, MPI_DOUBLE, MPI_MIN, 0, Profiler->comm);
    MPI_Reduce(Times, MaxTimes, Profiler->num_timers, MPI_DOUBLE, MPI_MAX, 0, Profiler->comm);
    MPI_Reduce(Times, AvgTimes, Profiler->num_timers, MPI_DOUBLE, MPI_SUM, 0, Profiler->comm);

    if (Profiler->comm_rank == 0) {

      for (iTimer = 0; iTimer < Profiler->num_timers; ++iTimer) {
        AvgTimes[iTimer] /= (double)Profiler->comm_size;
      }

      for (iTimer = 0; iTimer < Profiler->num_timers; ++iTimer) {
        const t_profiler_timer_info *TimerInfo = Profiler->timers[iTimer];
        fprintf(File, "%s: %f %f %f\n", TimerInfo->name, MinTimes[iTimer], MaxTimes[iTimer],
          AvgTimes[iTimer]);
      }

      fflush(File);

      free(MinTimes);
      free(MaxTimes);
      free(AvgTimes);

    }

    free(Times);

  }

}
