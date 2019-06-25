// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RECV_MAP_HPP_INCLUDED
#define OVK_CORE_RECV_MAP_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/FloatingRef.hpp>
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
  recv_map(long long NumValues, array<long long> RecvOrder, array_view<const int> SourceRanks);

  floating_ref<const recv_map> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<recv_map> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  long long Count() const { return RecvOrder_.Count(); }

  const array<recv> &Recvs() const { return Recvs_; }

  const array<long long> &RecvOrder() const { return RecvOrder_; }
  const array<int> &RecvIndices() const { return RecvIndices_; }

private:

  floating_ref_generator FloatingRefGenerator_;

  array<recv> Recvs_;
  array<long long> RecvOrder_;
  array<int> RecvIndices_;

};

}}

#endif
