// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/SendMap.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Profiler.hpp"

#include <mpi.h>

#include <map>
#include <utility>

namespace ovk {
namespace core {

send_map::send_map(long long NumValues, array<long long> SendOrder, const array<int>
  &DestinationRanks):
  SendOrder_(std::move(SendOrder))
{

  SendIndices_.Resize({NumValues}, -1);

  std::map<int, long long> SendCounts;

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    if (SendOrder_(iValue) >= 0) {
      int Rank = DestinationRanks(iValue);
      auto Iter = SendCounts.lower_bound(Rank);
      if (Iter != SendCounts.end() && Iter->first == Rank) {
        ++Iter->second;
      } else {
        SendCounts.emplace_hint(Iter, Rank, 1);
      }
    }
  }

  int NumSends = SendCounts.size();

  for (auto &Pair : SendCounts) {
    send &Send = Sends_.Append();
    Send.Rank = Pair.first;
    Send.NumValues = Pair.second;
  }

  SendCounts.clear();

  std::map<int, int> RankToSendIndex;

  for (int iSend = 0; iSend < NumSends; ++iSend) {
    RankToSendIndex.emplace(Sends_(iSend).Rank, iSend);
  }

  SendIndices_.Resize({NumValues}, -1);

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    if (SendOrder_(iValue) >= 0) {
      int Rank = DestinationRanks(iValue);
      SendIndices_(iValue) = RankToSendIndex[Rank];
    }
  }

}

}}
