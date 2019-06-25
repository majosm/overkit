// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_MAP_HPP_INCLUDED
#define OVK_CORE_COLLECT_MAP_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Partition.hpp>

#include <mpi.h>

namespace ovk {
namespace core {

class collect_map {

public:

  struct send {
    int Rank;
    long long NumPoints;
    array<int,2> Points;
  };

  struct recv {
    int Rank;
    long long NumPoints;
  };

  collect_map() = default;
  collect_map(const cart &Cart, const partition &Partition, array<int,3> CellExtents);

  floating_ref<const collect_map> GetFloatingRef() const {
    return FloatingRefGenerator_.Generate(*this);
  }
  floating_ref<collect_map> GetFloatingRef() { return FloatingRefGenerator_.Generate(*this); }

  long long Count() const { return CellExtents_.Size(2); }

  const array<int,3> &CellExtents() const { return CellExtents_; }

  int MaxVertices() const { return MaxVertices_; }

  const array<send> &Sends() const { return Sends_; }
  const array<recv> &Recvs() const { return Recvs_; }

  const array<int> &RemoteVertexCounts() const { return NumRemoteVertices_; }
  const array<long long *> &RemoteVertices() const { return RemoteVertices_; }
  const array<int *> &RemoteVertexRecvs() const { return RemoteVertexRecvs_; }
  const array<long long *> &RemoteVertexRecvBufferIndices() const {
    return RemoteVertexRecvBufferIndices_;
  }

private:

  floating_ref_generator FloatingRefGenerator_;

  array<int,3> CellExtents_;
  int MaxVertices_ = 0;
  array<send> Sends_;
  array<recv> Recvs_;
  array<int> NumRemoteVertices_;
  array<long long *> RemoteVertices_;
  array<long long> RemoteVerticesData_;
  array<int *> RemoteVertexRecvs_;
  array<int> RemoteVertexRecvsData_;
  array<long long *> RemoteVertexRecvBufferIndices_;
  array<long long> RemoteVertexRecvBufferIndicesData_;

  void CreateSendData_(const cart &Cart, const partition &Partition);
  void CreateRecvData_(const cart &Cart, const partition &Partition);

};

}}

#endif
