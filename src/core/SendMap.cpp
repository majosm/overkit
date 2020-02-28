// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/SendMap.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Map.hpp"
#include "ovk/core/Profiler.hpp"

#include <mpi.h>

#include <utility>

namespace ovk {
namespace core {

send_map::send_map(array_view<const int> DestinationRanks):
  send_map(DestinationRanks, MakeDefaultSendOrder_(DestinationRanks.Count()))
{}

send_map::send_map(array_view<const int> DestinationRanks, array<long long> SendOrder):
  SendOrder_(std::move(SendOrder))
{

  long long NumValues = DestinationRanks.Count();

  SendIndices_.Resize({NumValues}, -1);

  map<int,long long> SendCounts;

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    int Rank = DestinationRanks(iValue);
    if (Rank >= 0) {
      ++SendCounts.Fetch(Rank, 0);
    }
  }

  int NumSends = SendCounts.Count();

  for (auto &Entry : SendCounts) {
    send &Send = Sends_.Append();
    Send.Rank = Entry.Key();
    Send.NumValues = Entry.Value();
  }

  SendCounts.Clear();

  map<int,int> RankToSendIndex;
  RankToSendIndex.Reserve(NumSends);

  for (int iSend = 0; iSend < NumSends; ++iSend) {
    RankToSendIndex.Insert(Sends_(iSend).Rank, iSend);
  }

  SendIndices_.Resize({NumValues}, -1);

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    int Rank = DestinationRanks(iValue);
    if (Rank >= 0) {
      SendIndices_(iValue) = RankToSendIndex(Rank);
    }
  }

}

array<long long> send_map::MakeDefaultSendOrder_(long long NumValues) {

  array<long long> SendOrder({NumValues});

  for (long long iValue = 0; iValue < NumValues; ++iValue) {
    SendOrder(iValue) = iValue;
  }

  return SendOrder;

}

}}
