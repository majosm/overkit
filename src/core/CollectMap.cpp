// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/CollectMap.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Cart.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/FloatingRef.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Indexer.hpp"
#include "ovk/core/Partition.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/ScalarOps.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <utility>

namespace ovk {
namespace core {

collect_map::collect_map(const cart &Cart, const partition &Partition, array<int,3> CellExtents):
  CellExtents_(std::move(CellExtents))
{

  MaxVertices_ = 0;
  for (long long iCell = 0; iCell < CellExtents_.Size(2); ++iCell) {
    int NumVertices = 1;
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      NumVertices *= CellExtents_(1,iDim,iCell) - CellExtents_(0,iDim,iCell);
    }
    MaxVertices_ = Max(MaxVertices_, NumVertices);
  }

  CreateSendData_(Cart, Partition);
  CreateRecvData_(Cart, Partition);

}

void collect_map::CreateSendData_(const cart &Cart, const partition &Partition) {

  long long NumCells = CellExtents_.Size(2);

  if (NumCells > 0) {

    int NumDims = Cart.Dimension();
    const range &GlobalRange = Cart.Range();
    const range &LocalRange = Partition.LocalRange();
    const array<core::partition_info> &Neighbors = Partition.Neighbors();
    int NumNeighbors = Neighbors.Count();

    array<range> SendToNeighborRanges({NumNeighbors});
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      SendToNeighborRanges(iNeighbor) = MakeEmptyRange(NumDims);
    }

    for (long long iCell = 0; iCell < NumCells; ++iCell) {
      range CellRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        CellRange.Begin(iDim) = CellExtents_(0,iDim,iCell);
        CellRange.End(iDim) = CellExtents_(1,iDim,iCell);
      }
      bool AwayFromEdge = GlobalRange.Includes(CellRange);
      for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
        if (AwayFromEdge) {
          bool Overlaps = RangesOverlap(Neighbors(iNeighbor).LocalRange, CellRange);
          if (Overlaps) {
            SendToNeighborRanges(iNeighbor) = UnionRanges(SendToNeighborRanges(iNeighbor),
              IntersectRanges(LocalRange, CellRange));
          }
        } else {
          bool Overlaps = false;
          for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
            for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
              for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
                tuple<int> Vertex = Cart.PeriodicAdjust({i,j,k});
                if (Neighbors(iNeighbor).LocalRange.Contains(Vertex)) {
                  Overlaps = true;
                  goto done_checking_for_overlap1;
                }
              }
            }
          }
          done_checking_for_overlap1:;
          if (Overlaps) {
            for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
              for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
                for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
                  tuple<int> Vertex = Cart.PeriodicAdjust({i,j,k});
                  if (LocalRange.Contains(Vertex)) {
                    SendToNeighborRanges(iNeighbor) = ExtendRange(SendToNeighborRanges(iNeighbor),
                      Vertex);
                  }
                }
              }
            }
          }
        }
      }
    }

    using range_indexer = indexer<long long, int, MAX_DIMS, array_layout::COLUMN_MAJOR>;
    array<range_indexer> SendToNeighborIndexers({NumNeighbors});
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      const range &SendToNeighborRange = SendToNeighborRanges(iNeighbor);
      SendToNeighborIndexers(iNeighbor) = range_indexer(SendToNeighborRange);
    }

    array<int> SendIndexToNeighbor;
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!SendToNeighborRanges(iNeighbor).Empty()) {
        SendIndexToNeighbor.Append(iNeighbor);
      }
    }
    int NumSends = SendIndexToNeighbor.Count();

    Sends_.Resize({NumSends});

    for (int iSend = 0; iSend < NumSends; ++iSend) {
      int iNeighbor = SendIndexToNeighbor(iSend);
      Sends_(iSend).Rank = Neighbors(iNeighbor).Rank;
    }

    array<array<bool>> SendMasks({NumSends});
    for (int iSend = 0; iSend < NumSends; ++iSend) {
      int iNeighbor = SendIndexToNeighbor(iSend);
      SendMasks(iSend).Resize({SendToNeighborRanges(iNeighbor).Count()}, false);
    }

    for (long long iCell = 0; iCell < NumCells; ++iCell) {
      range CellRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        CellRange.Begin(iDim) = CellExtents_(0,iDim,iCell);
        CellRange.End(iDim) = CellExtents_(1,iDim,iCell);
      }
      bool AwayFromEdge = GlobalRange.Includes(CellRange);
      for (int iSend = 0; iSend < NumSends; ++iSend) {
        int iNeighbor = SendIndexToNeighbor(iSend);
        const range_indexer &Indexer = SendToNeighborIndexers(iNeighbor);
        if (AwayFromEdge) {
          bool Overlaps = RangesOverlap(Neighbors(iNeighbor).LocalRange, CellRange);
          if (Overlaps) {
            range LocalCellRange = IntersectRanges(LocalRange, CellRange);
            for (int k = LocalCellRange.Begin(2); k < LocalCellRange.End(2); ++k) {
              for (int j = LocalCellRange.Begin(1); j < LocalCellRange.End(1); ++j) {
                for (int i = LocalCellRange.Begin(0); i < LocalCellRange.End(0); ++i) {
                  long long iPoint = Indexer.ToIndex(i,j,k);
                  SendMasks(iSend)(iPoint) = true;
                }
              }
            }
          }
        } else {
          bool Overlaps = false;
          for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
            for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
              for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
                tuple<int> Vertex = Cart.PeriodicAdjust({i,j,k});
                if (Neighbors(iNeighbor).LocalRange.Contains(Vertex)) {
                  Overlaps = true;
                  goto done_checking_for_overlap2;
                }
              }
            }
          }
          done_checking_for_overlap2:;
          if (Overlaps) {
            for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
              for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
                for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
                  tuple<int> Vertex = Cart.PeriodicAdjust({i,j,k});
                  if (LocalRange.Contains(Vertex)) {
                    long long iPoint = Indexer.ToIndex(Vertex);
                    SendMasks(iSend)(iPoint) = true;
                  }
                }
              }
            }
          }
        }
      }
    }

    for (int iSend = 0; iSend < NumSends; ++iSend) {
      send &Send = Sends_(iSend);
      int iNeighbor = SendIndexToNeighbor(iSend);
      Send.NumPoints = 0;
      for (long long iPoint = 0; iPoint < SendToNeighborRanges(iNeighbor).Count(); ++iPoint) {
        if (SendMasks(iSend)(iPoint)) {
          ++Send.NumPoints;
        }
      }
    }

    for (int iSend = 0; iSend < NumSends; ++iSend) {
      send &Send = Sends_(iSend);
      int iNeighbor = SendIndexToNeighbor(iSend);
      Send.Points.Resize({{MAX_DIMS,Send.NumPoints}});
      const range_indexer &Indexer = SendToNeighborIndexers(iNeighbor);
      long long iSendPoint = 0;
      for (long long iPoint = 0; iPoint < SendToNeighborRanges(iNeighbor).Count(); ++iPoint) {
        if (SendMasks(iSend)(iPoint)) {
          tuple<int> Point = Indexer.ToTuple(iPoint);
          for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
            Send.Points(iDim,iSendPoint) = Point(iDim);
          }
          ++iSendPoint;
        }
      }
    }

  }

}

void collect_map::CreateRecvData_(const cart &Cart, const partition &Partition) {

  long long NumCells = CellExtents_.Size(2);

  if (NumCells > 0) {

    int NumDims = Cart.Dimension();
    const range &GlobalRange = Cart.Range();
    const range &LocalRange = Partition.LocalRange();
    const array<core::partition_info> &Neighbors = Partition.Neighbors();
    int NumNeighbors = Neighbors.Count();

    array<range> RecvFromNeighborRanges({NumNeighbors});
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      RecvFromNeighborRanges(iNeighbor) = MakeEmptyRange(NumDims);
    }

    for (long long iCell = 0; iCell < NumCells; ++iCell) {
      range CellRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        CellRange.Begin(iDim) = CellExtents_(0,iDim,iCell);
        CellRange.End(iDim) = CellExtents_(1,iDim,iCell);
      }
      bool AwayFromEdge = GlobalRange.Includes(CellRange);
      for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
        if (AwayFromEdge) {
          RecvFromNeighborRanges(iNeighbor) = UnionRanges(RecvFromNeighborRanges(iNeighbor),
            IntersectRanges(Neighbors(iNeighbor).LocalRange, CellRange));
        } else {
          for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
            for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
              for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
                tuple<int> Vertex = Cart.PeriodicAdjust({i,j,k});
                if (Neighbors(iNeighbor).LocalRange.Contains(Vertex)) {
                  RecvFromNeighborRanges(iNeighbor) = ExtendRange(RecvFromNeighborRanges(iNeighbor),
                    Vertex);
                }
              }
            }
          }
        }
      }
    }

    using range_indexer = indexer<long long, int, MAX_DIMS, array_layout::COLUMN_MAJOR>;
    array<range_indexer> RecvFromNeighborIndexers({NumNeighbors});
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      const range &RecvFromNeighborRange = RecvFromNeighborRanges(iNeighbor);
      RecvFromNeighborIndexers(iNeighbor) = range_indexer(RecvFromNeighborRange);
    }

    array<int> RecvIndexToNeighbor;
    for (int iNeighbor = 0; iNeighbor < NumNeighbors; ++iNeighbor) {
      if (!RecvFromNeighborRanges(iNeighbor).Empty()) {
        RecvIndexToNeighbor.Append(iNeighbor);
      }
    }
    int NumRecvs = RecvIndexToNeighbor.Count();

    Recvs_.Resize({NumRecvs});

    for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
      int iNeighbor = RecvIndexToNeighbor(iRecv);
      Recvs_(iRecv).Rank = Neighbors(iNeighbor).Rank;
    }

    array<array<bool>> RecvMasks({NumRecvs});
    for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
      int iNeighbor = RecvIndexToNeighbor(iRecv);
      RecvMasks(iRecv).Resize({RecvFromNeighborRanges(iNeighbor).Count()}, false);
    }

    for (long long iCell = 0; iCell < NumCells; ++iCell) {
      range CellRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        CellRange.Begin(iDim) = CellExtents_(0,iDim,iCell);
        CellRange.End(iDim) = CellExtents_(1,iDim,iCell);
      }
      bool AwayFromEdge = GlobalRange.Includes(CellRange);
      for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
        int iNeighbor = RecvIndexToNeighbor(iRecv);
        const range_indexer &Indexer = RecvFromNeighborIndexers(iNeighbor);
        if (AwayFromEdge) {
          range RemoteCellRange = IntersectRanges(Neighbors(iNeighbor).LocalRange, CellRange);
          for (int k = RemoteCellRange.Begin(2); k < RemoteCellRange.End(2); ++k) {
            for (int j = RemoteCellRange.Begin(1); j < RemoteCellRange.End(1); ++j) {
              for (int i = RemoteCellRange.Begin(0); i < RemoteCellRange.End(0); ++i) {
                long long iPoint = Indexer.ToIndex(i,j,k);
                RecvMasks(iRecv)(iPoint) = true;
              }
            }
          }
        } else {
          for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
            for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
              for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
                tuple<int> Vertex = Cart.PeriodicAdjust({i,j,k});
                if (Neighbors(iNeighbor).LocalRange.Contains(Vertex)) {
                  long long iPoint = Indexer.ToIndex(Vertex);
                  RecvMasks(iRecv)(iPoint) = true;
                }
              }
            }
          }
        }
      }
    }

    for (int iRecv = 0; iRecv < NumRecvs; ++iRecv) {
      recv &Recv = Recvs_(iRecv);
      int iNeighbor = RecvIndexToNeighbor(iRecv);
      Recv.NumPoints = 0;
      for (long long iPoint = 0; iPoint < RecvFromNeighborRanges(iNeighbor).Count(); ++iPoint) {
        if (RecvMasks(iRecv)(iPoint)) {
          ++Recv.NumPoints;
        }
      }
    }

    NumRemoteVertices_.Resize({NumCells}, 0);
    RemoteVertices_.Resize({NumCells}, nullptr);
    RemoteVertexRecvs_.Resize({NumCells}, nullptr);
    RemoteVertexRecvBufferIndices_.Resize({NumCells}, nullptr);

    long long TotalRemoteVertices = 0;
    for (long long iCell = 0; iCell < NumCells; ++iCell) {
      range CellRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        CellRange.Begin(iDim) = CellExtents_(0,iDim,iCell);
        CellRange.End(iDim) = CellExtents_(1,iDim,iCell);
      }
      bool AwayFromEdge = GlobalRange.Includes(CellRange);
      int NumRemoteVertices;
      if (AwayFromEdge) {
        range LocalCellRange = IntersectRanges(LocalRange, CellRange);
        NumRemoteVertices = CellRange.Count() - LocalCellRange.Count();
      } else {
        NumRemoteVertices = 0;
        for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
          for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
            for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
              tuple<int> Vertex = Cart.PeriodicAdjust({i,j,k});
              if (!LocalRange.Contains(Vertex)) {
                ++NumRemoteVertices;
              }
            }
          }
        }
      }
      NumRemoteVertices_(iCell) = NumRemoteVertices;
      TotalRemoteVertices += NumRemoteVertices;
    }

    RemoteVerticesData_.Resize({TotalRemoteVertices});
    RemoteVertexRecvsData_.Resize({TotalRemoteVertices});
    RemoteVertexRecvBufferIndicesData_.Resize({TotalRemoteVertices});

    long long Offset = 0;
    for (long long iCell = 0; iCell < NumCells; ++iCell) {
      long long NumRemoteVertices = NumRemoteVertices_(iCell);
      RemoteVertices_(iCell) = RemoteVerticesData_.Data(Offset);
      RemoteVertexRecvs_(iCell) = RemoteVertexRecvsData_.Data(Offset);
      RemoteVertexRecvBufferIndices_(iCell) = RemoteVertexRecvBufferIndicesData_.Data(Offset);
      Offset += NumRemoteVertices;
    }

    array<array<long long>> RecvBufferIndices({Recvs_.Count()});
    for (int iRecv = 0; iRecv < Recvs_.Count(); ++iRecv) {
      int iNeighbor = RecvIndexToNeighbor(iRecv);
      long long NumPoints = RecvFromNeighborRanges(iNeighbor).Count();
      RecvBufferIndices(iRecv).Resize({NumPoints}, -1);
      long long iRemotePoint = 0;
      for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
        if (RecvMasks(iRecv)(iPoint)) {
          RecvBufferIndices(iRecv)(iPoint) = iRemotePoint;
          ++iRemotePoint;
        }
      }
    }

    array<int> CellRecvs({MaxVertices_});
    array<long long> CellRecvBufferIndices({MaxVertices_});

    using donor_indexer = indexer<int, int, MAX_DIMS, array_layout::COLUMN_MAJOR>;

    for (long long iCell = 0; iCell < NumCells; ++iCell) {
      for (int iVertex = 0; iVertex < MaxVertices_; ++iVertex) {
        CellRecvs(iVertex) = -1;
        CellRecvBufferIndices(iVertex) = -1;
      }
      range CellRange;
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        CellRange.Begin(iDim) = CellExtents_(0,iDim,iCell);
        CellRange.End(iDim) = CellExtents_(1,iDim,iCell);
      }
      donor_indexer CellIndexer(CellRange);
      bool AwayFromEdge = GlobalRange.Includes(CellRange);
      for (int iRecv = 0; iRecv < Recvs_.Count(); ++iRecv) {
        int iNeighbor = RecvIndexToNeighbor(iRecv);
        const range_indexer &RecvFromNeighborIndexer = RecvFromNeighborIndexers(iNeighbor);
        if (AwayFromEdge) {
          range RemoteCellRange = IntersectRanges(Neighbors(iNeighbor).LocalRange, CellRange);
          for (int k = RemoteCellRange.Begin(2); k < RemoteCellRange.End(2); ++k) {
            for (int j = RemoteCellRange.Begin(1); j < RemoteCellRange.End(1); ++j) {
              for (int i = RemoteCellRange.Begin(0); i < RemoteCellRange.End(0); ++i) {
                int iVertex = CellIndexer.ToIndex(i,j,k);
                long long iPoint = RecvFromNeighborIndexer.ToIndex(i,j,k);
                CellRecvs(iVertex) = iRecv;
                CellRecvBufferIndices(iVertex) = RecvBufferIndices(iRecv)
                  (iPoint);
              }
            }
          }
        } else {
          int iVertex = 0;
          for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
            for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
              for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
                tuple<int> Point = Cart.PeriodicAdjust({i,j,k});
                if (Neighbors(iNeighbor).LocalRange.Contains(Point)) {
                  long long iPoint = RecvFromNeighborIndexer.ToIndex(Point);
                  CellRecvs(iVertex) = iRecv;
                  CellRecvBufferIndices(iVertex) = RecvBufferIndices
                    (iRecv)(iPoint);
                }
                ++iVertex;
              }
            }
          }
        }
      }
      int iRemoteVertex = 0;
      int iVertex = 0;
      for (int k = CellRange.Begin(2); k < CellRange.End(2); ++k) {
        for (int j = CellRange.Begin(1); j < CellRange.End(1); ++j) {
          for (int i = CellRange.Begin(0); i < CellRange.End(0); ++i) {
            if (CellRecvs(iVertex) >= 0) {
              RemoteVertices_(iCell)[iRemoteVertex] = iVertex;
              RemoteVertexRecvs_(iCell)[iRemoteVertex] = CellRecvs(iVertex);
              RemoteVertexRecvBufferIndices_(iCell)[iRemoteVertex] = CellRecvBufferIndices(iVertex);
              ++iRemoteVertex;
            }
            ++iVertex;
          }
        }
      }
    }

  }

}

}}
