// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/RecvMap.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Profiler.hpp"

#include <mpi.h>

#include <utility>

namespace ovk {
namespace core {

recv_map::recv_map(array_view<const int> SourceRanks):
  recv_map(SourceRanks, MakeDefaultRecvOrder_(SourceRanks.Count()))
{}

recv_map::recv_map(array_view<const int> SourceRanks, array<long long> RecvOrder):
  RecvOrder_(std::move(RecvOrder))
{

  long long NumValues = SourceRanks.Count();

  RecvIndices_.Resize({NumValues}, -1);

  map<int,long long> RecvCounts;

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    int Rank = SourceRanks(iValue);
    if (Rank >= 0) {
      ++RecvCounts.Fetch(Rank, 0);
    }
  }

  int NumRecvs = RecvCounts.Count();

  for (auto &Entry : RecvCounts) {
    recv &Recv = Recvs_.Append();
    Recv.Rank = Entry.Key();
    Recv.NumValues = Entry.Value();
  }

  RecvCounts.Clear();

  map<int,int> RankToRecvIndex;
  RankToRecvIndex.Reserve(NumRecvs);

  for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
    RankToRecvIndex.Insert(Recvs_(iRecv).Rank, iRecv);
  }

  RecvIndices_.Resize({NumValues}, -1);

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    int Rank = SourceRanks(iValue);
    if (Rank >= 0) {
      RecvIndices_(iValue) = RankToRecvIndex(Rank);
    }
  }

}

array<long long> recv_map::MakeDefaultRecvOrder_(long long NumValues) {

  array<long long> RecvOrder({NumValues});

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    RecvOrder(iValue) = iValue;
  }

  return RecvOrder;

}

}}
