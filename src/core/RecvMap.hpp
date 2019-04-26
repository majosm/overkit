// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RECV_MAP_HPP_INCLUDED
#define OVK_CORE_RECV_MAP_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Global.hpp>

namespace ovk {
namespace core {

class recv_map {

public:

  struct recv {
    int Rank;
    long long NumValues;
  };

  recv_map() = default;
  recv_map(long long NumValues, array<long long> RecvOrder, const array<int> &SourceRanks);

  long long Count() const { return RecvOrder_.Count(); }

  const array<recv> &Recvs() const { return Recvs_; }

  const array<long long> &RecvOrder() const { return RecvOrder_; }
  const array<int> &RecvIndices() const { return RecvIndices_; }

private:

  array<recv> Recvs_;
  array<long long> RecvOrder_;
  array<int> RecvIndices_;

};

}}

#endif
