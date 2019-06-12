// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/RecvMap.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Profiler.hpp"

#include <mpi.h>

#include <map>
#include <utility>

namespace ovk {
namespace core {

recv_map::recv_map(long long NumValues, array<long long> RecvOrder, const array<int> &SourceRanks):
  RecvOrder_(std::move(RecvOrder))
{

  RecvIndices_.Resize({NumValues}, -1);

  std::map<int, long long> RecvCounts;

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    if (RecvOrder_(iValue) >= 0) {
      int Rank = SourceRanks(iValue);
      auto Iter = RecvCounts.lower_bound(Rank);
      if (Iter != RecvCounts.end() && Iter->first == Rank) {
        ++Iter->second;
      } else {
        RecvCounts.emplace_hint(Iter, Rank, 1);
      }
    }
  }

  int NumRecvs = RecvCounts.size();

  for (auto &Pair : RecvCounts) {
    recv &Recv = Recvs_.Append();
    Recv.Rank = Pair.first;
    Recv.NumValues = Pair.second;
  }

  RecvCounts.clear();

  std::map<int, int> RankToRecvIndex;

  for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
    RankToRecvIndex.emplace(Recvs_(iRecv).Rank, iRecv);
  }

  RecvIndices_.Resize({NumValues}, -1);

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    if (RecvOrder_(iValue) >= 0) {
      int Rank = SourceRanks(iValue);
      RecvIndices_(iValue) = RankToRecvIndex[Rank];
    }
  }

}

}}