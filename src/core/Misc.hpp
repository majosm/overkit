// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MISC_HPP_INCLUDED
#define OVK_CORE_MISC_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>

#include <mpi.h>

#include <algorithm>
#include <string>

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
  unsigned char SendBuffer_[2] = {0, 0};
  unsigned char RecvBuffer_[1] = {0};
  MPI_Request Requests_[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
  int NumCompleted_ = 0;
#endif

};

// Given known list of ranks on one end of communication, generate list of ranks on other end
array<int> DynamicHandshake(comm_view Comm, array_view<const int> Ranks);

// Generate permutation corresponding to ascending-order sort
template <typename ArrayType, OVK_FUNCDECL_REQUIRES(IsArray<ArrayType>() && ArrayRank<ArrayType>()
  == 1)> void SortPermutation(const ArrayType &Array, array_view<long long> Permutation);

}}

#include <ovk/core/Misc.inl>

#endif
