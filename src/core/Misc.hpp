// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MISC_HPP_INCLUDED
#define OVK_CORE_MISC_HPP_INCLUDED

#include <ovk/core/Global.hpp>

#include <mpi.h>

#include <algorithm>
#include <vector>

namespace ovk {
namespace core {

// Like MPI_Bcast, but can be used when source rank is not known globally
void BroadcastAnySource(void *Data, int Count, MPI_Datatype DataType, bool IsSource, MPI_Comm Comm);

// Wrapper around MPI_Ibarrier (or crappy alternative if MPI_Ibarrier is not supported) 
// representing a global flag that gets set only after all processes call a function
// (StartSignal) to set it
struct signal {
  MPI_Comm Comm_;
#ifdef OVK_HAVE_MPI_IBARRIER
  MPI_Request Request_;
#else
  char SendBuffer_[2], RecvBuffer_[1];
  int CommSize_;
  int CommRank_;
  MPI_Request Requests_[2];
  int NumCompleted_;
#endif
};

void CreateSignal(signal &Signal, MPI_Comm Comm);
void StartSignal(signal &Signal);
void CheckSignal(signal &Signal, bool &Done);
void DestroySignal(signal &Signal);

// Given known list of ranks sending to, get list of ranks receiving from
void DynamicHandshake(MPI_Comm Comm, int NumDestRanks, const std::vector<int> &DestRanks,
  std::vector<int> &SourceRanks);

// Generate permutation corresponding to ascending-order sort
template <typename T> void SortPermutation(long long N, const T *Array, long long *Permutation);

// Sometimes want to avoid bool in templates (types containing std::vector<T>, MPI sends/recvs, etc.)
template <typename T> struct no_bool_ { using type = T; };
template <> struct no_bool_<bool> { using type = unsigned char; };
template <typename T> using no_bool = typename no_bool_<T>::type;

}}

#include <ovk/core/Misc.inl>

#endif
