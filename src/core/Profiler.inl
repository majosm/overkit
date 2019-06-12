// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

OVK_FORCE_INLINE void timer::Start() {

  if (!Active_) {
    LastStartTime_ = MPI_Wtime();
  }

  Active_ = true;

}

OVK_FORCE_INLINE void timer::Stop() {

  if (Active_) {
    LastEndTime_ = MPI_Wtime();
    AccumulatedTime_ += LastEndTime_ - LastStartTime_;
  }

  Active_ = false;

}

inline void timer::Reset() {

  *this = timer();

}

OVK_FORCE_INLINE double timer::Elapsed() const {

  if (Active_) {
    return MPI_Wtime() - LastStartTime_;
  } else {
    return LastEndTime_ - LastStartTime_;
  }

}

OVK_FORCE_INLINE double timer::Accumulated() const {

  if (Active_) {
    return AccumulatedTime_ + (MPI_Wtime() - LastStartTime_);
  } else {
    return AccumulatedTime_;
  }

}

OVK_FORCE_INLINE void profiler::Start(int TimerID) {

  if (Enabled_) {
    // Don't want to force everything inline, just the if statement
    Start_(TimerID);
  }

}

OVK_FORCE_INLINE void profiler::StartSync(int TimerID, MPI_Comm Comm) {

  if (Enabled_) {
    // Don't want to force everything inline, just the if statement
    StartSync_(TimerID, Comm);
  }

}

OVK_FORCE_INLINE void profiler::Stop(int TimerID) {

  if (Enabled_) {
    // Don't want to force everything inline, just the if statement
    Stop_(TimerID);
  }

}

}}
