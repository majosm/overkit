// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Halo.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/Field.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

namespace halo_internal {

halo_map::halo_map(const cart &Cart, const range &LocalRange, const range &ExtendedRange, const
  map<int,decomp_info> &Neighbors) {

  const range &GlobalRange = Cart.Range();
  int NumNeighbors = Neighbors.Count();

  field_indexer ExtendedIndexer(ExtendedRange);

  NeighborRanks_.Resize({NumNeighbors});
  NeighborSendIndices_.Resize({NumNeighbors});
  NeighborRecvIndices_.Resize({NumNeighbors});

  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    NeighborRanks_(iNeighbor) = Neighbors[iNeighbor].Key();
  }

  for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
    const range &NeighborLocalRange = Neighbors[iNeighbor].Value().LocalRange;
    const range &NeighborExtendedRange = Neighbors[iNeighbor].Value().ExtendedRange;
    array<long long> &SendIndices = NeighborSendIndices_(iNeighbor);
    if (Cart.Range().Includes(NeighborExtendedRange)) {
      range SendRange = IntersectRanges(NeighborExtendedRange, LocalRange);
      SendIndices.Reserve(SendRange.Count());
      for (int k = SendRange.Begin(2); k < SendRange.End(2); ++k) {
        for (int j = SendRange.Begin(1); j < SendRange.End(1); ++j) {
          for (int i = SendRange.Begin(0); i < SendRange.End(0); ++i) {
            long long iPoint = ExtendedIndexer.ToIndex(i,j,k);
            SendIndices.Append(iPoint);
          }
        }
      }
    } else {
      field<bool> SendMask(NeighborExtendedRange, false);
      long long NumSendPoints = 0;
      for (int k = NeighborExtendedRange.Begin(2); k < NeighborExtendedRange.End(2); ++k) {
        for (int j = NeighborExtendedRange.Begin(1); j < NeighborExtendedRange.End(1); ++j) {
          for (int i = NeighborExtendedRange.Begin(0); i < NeighborExtendedRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (!NeighborLocalRange.Contains(Point)) {
              if (Cart.Range().Contains(Point)) {
                SendMask(Point) = LocalRange.Contains(Point);
              } else {
                tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
                SendMask(Point) = LocalRange.Contains(AdjustedPoint);
              }
              if (SendMask(Point)) {
                ++NumSendPoints;
              }
            }
          }
        }
      }
      SendIndices.Reserve(NumSendPoints);
      for (int k = NeighborExtendedRange.Begin(2); k < NeighborExtendedRange.End(2); ++k) {
        for (int j = NeighborExtendedRange.Begin(1); j < NeighborExtendedRange.End(1); ++j) {
          for (int i = NeighborExtendedRange.Begin(0); i < NeighborExtendedRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (SendMask(Point)) {
              long long iPoint;
              if (Cart.Range().Contains(Point)) {
                iPoint = ExtendedIndexer.ToIndex(Point);
              } else {
                tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
                iPoint = ExtendedIndexer.ToIndex(AdjustedPoint);
              }
              SendIndices.Append(iPoint);
            }
          }
        }
      }
    }
  }

  if (GlobalRange.Includes(ExtendedRange)) {
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      const range &NeighborLocalRange = Neighbors[iNeighbor].Value().LocalRange;
      array<long long> &RecvIndices = NeighborRecvIndices_(iNeighbor);
      range RecvRange = IntersectRanges(ExtendedRange, NeighborLocalRange);
      RecvIndices.Reserve(RecvRange.Count());
      for (int k = RecvRange.Begin(2); k < RecvRange.End(2); ++k) {
        for (int j = RecvRange.Begin(1); j < RecvRange.End(1); ++j) {
          for (int i = RecvRange.Begin(0); i < RecvRange.End(0); ++i) {
            long long iPoint = ExtendedIndexer.ToIndex(i,j,k);
            RecvIndices.Append(iPoint);
          }
        }
      }
    }
  } else {
    field<bool> RecvMask(ExtendedRange);
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      RecvMask.Fill(false);
      const range &NeighborLocalRange = Neighbors[iNeighbor].Value().LocalRange;
      array<long long> &RecvIndices = NeighborRecvIndices_(iNeighbor);
      long long NumRecvPoints = 0;
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (!LocalRange.Contains(Point)) {
              if (Cart.Range().Contains(Point)) {
                RecvMask(Point) = NeighborLocalRange.Contains(Point);
              } else {
                tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
                RecvMask(Point) = NeighborLocalRange.Contains(AdjustedPoint);
              }
              if (RecvMask(Point)) {
                ++NumRecvPoints;
              }
            }
          }
        }
      }
      RecvIndices.Reserve(NumRecvPoints);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (RecvMask(Point)) {
              long long iPoint = ExtendedIndexer.ToIndex(Point);
              RecvIndices.Append(iPoint);
            }
          }
        }
      }
    }
  }

  // Source and destination may sometimes be on the same rank (e.g. if for a given
  // direction periodic==true and #procs==1)
  if (!GlobalRange.Includes(ExtendedRange)) {
    field<bool> LocalToLocalMask(ExtendedRange, false);
    long long NumLocalToLocal = 0;
    for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
      for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
        for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          if (!LocalRange.Contains(Point)) {
            tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
            if (LocalRange.Contains(AdjustedPoint)) {
              ++NumLocalToLocal;
            }
          }
        }
      }
    }
    LocalToLocalSourceIndices_.Reserve(NumLocalToLocal);
    LocalToLocalDestIndices_.Reserve(NumLocalToLocal);
    for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
      for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
        for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
          tuple<int> Point = {i,j,k};
          if (!LocalRange.Contains(Point)) {
            tuple<int> AdjustedPoint = Cart.PeriodicAdjust(Point);
            if (LocalRange.Contains(AdjustedPoint)) {
              long long iSourcePoint = ExtendedIndexer.ToIndex(AdjustedPoint);
              long long iDestPoint = ExtendedIndexer.ToIndex(Point);
              LocalToLocalSourceIndices_.Append(iSourcePoint);
              LocalToLocalDestIndices_.Append(iDestPoint);
            }
          }
        }
      }
    }
  }

}

}

halo::halo(std::shared_ptr<context> Context, const cart &Cart, comm_view Comm, const range
  &LocalRange, const range &ExtendedRange, const map<int,decomp_info> &Neighbors):
  Context_(std::move(Context)),
  Comm_(Comm)
{

  profiler &Profiler = Context_->core_Profiler();

  Profiler.StartSync(TOTAL_TIME, Comm);
  Profiler.Start(SETUP_TIME);

  HaloMap_ = halo_map(Cart, LocalRange, ExtendedRange, Neighbors);

  Profiler.Stop(SETUP_TIME);
  Profiler.Stop(TOTAL_TIME);

}

}}
