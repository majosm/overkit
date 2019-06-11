// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PROFILER_HPP_INCLUDED
#define OVK_CORE_PROFILER_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IDMap.hpp>

#include <mpi.h>

#include <string>
#include <utility>

namespace ovk {
namespace core {

class timer {

public:

  timer() = default;

  OVK_FORCE_INLINE void Start();
  OVK_FORCE_INLINE void Stop();

  inline void Reset();

  OVK_FORCE_INLINE double Elapsed() const;

  OVK_FORCE_INLINE double Accumulated() const;

private:

  bool Active_ = false;
  double AccumulatedTime_ = 0.;
  double LastStartTime_;
  double LastEndTime_;

};

class profiler {

public:

  enum : int {
    HALO_TIME = 0,
    HALO_SETUP_TIME,
    HALO_EXCHANGE_TIME,
    HALO_EXCHANGE_PACK_TIME,
    HALO_EXCHANGE_MPI_TIME,
    HALO_EXCHANGE_UNPACK_TIME,
    EXCHANGER_COLLECT_TIME,
    EXCHANGER_COLLECT_MPI_TIME,
    EXCHANGER_COLLECT_PACK_TIME,
    EXCHANGER_COLLECT_REDUCE_TIME,
    EXCHANGER_SEND_RECV_TIME,
    EXCHANGER_SEND_RECV_PACK_TIME,
    EXCHANGER_SEND_RECV_MPI_TIME,
    EXCHANGER_SEND_RECV_UNPACK_TIME,
    EXCHANGER_DISPERSE_TIME,
    XINTOUT_IMPORT_TIME,
    XINTOUT_IMPORT_READ_TIME,
    XINTOUT_IMPORT_READ_MPI_IO_OPEN_TIME,
    XINTOUT_IMPORT_READ_MPI_IO_CLOSE_TIME,
    XINTOUT_IMPORT_READ_MPI_IO_READ_TIME,
    XINTOUT_IMPORT_READ_MPI_IO_OTHER_TIME,
    XINTOUT_IMPORT_MATCH_TIME,
    XINTOUT_IMPORT_MATCH_MAP_TO_BINS_TIME,
    XINTOUT_IMPORT_MATCH_HANDSHAKE_TIME,
    XINTOUT_IMPORT_MATCH_SEND_TO_BINS_TIME,
    XINTOUT_IMPORT_MATCH_FILL_CONNECTION_DATA_TIME,
    XINTOUT_IMPORT_MATCH_RECV_FROM_BINS_TIME,
    XINTOUT_IMPORT_MATCH_UNPACK_TIME,
    XINTOUT_IMPORT_DISTRIBUTE_TIME,
    XINTOUT_IMPORT_DISTRIBUTE_MAP_TO_BINS_TIME,
    XINTOUT_IMPORT_DISTRIBUTE_RETRIEVE_BINS_TIME,
    XINTOUT_IMPORT_DISTRIBUTE_FIND_RANKS_TIME,
    XINTOUT_IMPORT_DISTRIBUTE_HANDSHAKE_TIME,
    XINTOUT_IMPORT_DISTRIBUTE_SEND_DATA_TIME,
    XINTOUT_IMPORT_SET_CONNECTIVITIES_TIME,
    profiler_internal_TIMER_ID_COUNT
  };

  profiler() = default;
  explicit profiler(comm_view Comm);

  bool Enabled() const { return Enabled_; }

  void Enable();
  void Disable();

  OVK_FORCE_INLINE void Start(int TimerID);
  // Could use comm_view here, but want minimal overhead
  OVK_FORCE_INLINE void StartSync(int TimerID, MPI_Comm Comm);
  OVK_FORCE_INLINE void Stop(int TimerID);

  std::string WriteProfile() const;

private:

  struct timer_entry {
    timer Timer;
    int ActiveCount = 0;
    // XL doesn't think timer_entry is default-constructible without this (' = default' doesn't
    // work either)
    timer_entry() {}
  };

  comm_view Comm_ = MPI_COMM_SELF;
  bool Enabled_ = false;
  id_map<1,timer_entry> Timers_;

  void Start_(int TimerID);
  void StartSync_(int TimerID, MPI_Comm Comm);
  void Stop_(int TimerID);

  // Set non-contiguous because std::string is not noexcept movable until C++17
  static const id_map<1,std::string,false> TimerNames_;

};

}}

#include <ovk/core/Profiler.inl>

#endif
