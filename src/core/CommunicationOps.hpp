// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COMMUNICATION_OPS_HPP_INCLUDED
#define OVK_CORE_COMMUNICATION_OPS_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/ScopeGuard.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <string>
#include <utility>

namespace ovk {
namespace core {

void BroadcastString(std::string &String, int Root, comm_view Comm);

// Like MPI_Bcast, but can be used when source rank is not known globally
void BroadcastAnySource(void *Data, int Count, MPI_Datatype DataType, bool IsSource, comm_view
  Comm);
void BroadcastStringAnySource(std::string &String, bool IsSource, comm_view Comm);

// Wrapper around MPI_Ibarrier (or crappy alternative if MPI_Ibarrier is not supported) 
// representing a global flag that gets set only after all processes call a function
// (StartSignal) to set it
class signal {

public:

  signal(comm_view Comm);

  void Start();
  bool Check();

private:

  comm Comm_;
#ifdef OVK_HAVE_MPI_IBARRIER
  MPI_Request Request_ = MPI_REQUEST_NULL;
#else
  byte SendBuffer_[2] = {0, 0};
  byte RecvBuffer_[1] = {0};
  MPI_Request Requests_[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
  int NumCompleted_ = 0;
#endif

};

// Given known list of ranks on one end of communication, generate list of ranks on other end
array<int> DynamicHandshake(comm_view Comm, array_view<const int> Ranks);

// Run a section of code sequentially over each rank
template <typename F, OVK_FUNCDECL_REQUIRES(IsCallableWith<F &&>())> auto Serialize(comm_view Comm,
  F &&Func) -> decltype(std::forward<F>(Func)());

// Help identify processes that have gotten stuck somewhere
class hang_detector {

public:

  hang_detector(comm_view Comm, double Timeout=10.);

  void operator()();

private:

  comm Comm_;
  signal Signal_;
  double Timeout_;

};

}}

#include <ovk/core/CommunicationOps.inl>

#endif
