// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MISC_HPP_INCLUDED
#define OVK_CORE_MISC_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>

#include <mpi.h>

#include <algorithm>

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

// Given known list of ranks on one end of communication, generate list of ranks on other end
array<int> DynamicHandshake(MPI_Comm Comm, array_view<const int> Ranks);

// Generate permutation corresponding to ascending-order sort
template <typename ArrayType, OVK_FUNCDECL_REQUIRES(IsArray<ArrayType>() && ArrayRank<ArrayType>()
  == 1)> void SortPermutation(const ArrayType &Array, array_view<long long> Permutation);

}}

#include <ovk/core/Misc.inl>

#endif
