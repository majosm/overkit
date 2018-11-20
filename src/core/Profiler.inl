// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

namespace profiler_internal {

inline void ResetTimer(timer &Timer) {

  Timer.Active_ = false;
  Timer.LastStartTime_ = 0.;
  Timer.LastEndTime_ = 0.;
  Timer.AccumulatedTime_ = 0.;

}

inline void StartTimer(timer &Timer) {

  if (!Timer.Active_) {
    Timer.LastStartTime_ = MPI_Wtime();
  }
  Timer.Active_ = true;

}

inline void StopTimer(timer &Timer) {

  if (Timer.Active_) {
    Timer.LastEndTime_ = MPI_Wtime();
    Timer.AccumulatedTime_ += Timer.LastEndTime_ - Timer.LastStartTime_;
  }
  Timer.Active_ = false;

}

inline double GetElapsedTime(const timer &Timer) {

  if (Timer.Active_) {
    return MPI_Wtime() - Timer.LastStartTime_;
  } else {
    return Timer.LastEndTime_ - Timer.LastStartTime_;
  }

}

inline double GetAccumulatedTime(const timer &Timer) {

  if (Timer.Active_) {
    return Timer.AccumulatedTime_ + (MPI_Wtime() - Timer.LastStartTime_);
  } else {
    return Timer.AccumulatedTime_;
  }

}

}

inline void StartProfile(profiler &Profiler, int TimerID) {

  if (Profiler.Enabled_) {
    OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < int(Profiler.Timers_.size()), "Invalid timer ID.");
    StartTimer(Profiler.Timers_[TimerID].second);
  }

}

inline void StartProfileSync(profiler &Profiler, int TimerID, MPI_Comm Comm) {

  if (Profiler.Enabled_) {
    OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < int(Profiler.Timers_.size()), "Invalid timer ID.");
    MPI_Barrier(Comm);
    StartTimer(Profiler.Timers_[TimerID].second);
  }

}

inline void EndProfile(profiler &Profiler, int TimerID) {

  if (Profiler.Enabled_) {
    OVK_DEBUG_ASSERT(TimerID >= 0 && TimerID < int(Profiler.Timers_.size()), "Invalid timer ID.");
    StopTimer(Profiler.Timers_[TimerID].second);
  }

}

}}
