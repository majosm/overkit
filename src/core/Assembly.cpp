// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Domain.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/Collect.hpp"
#include "ovk/core/CollectMap.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/ConnectivityD.hpp"
#include "ovk/core/ConnectivityR.hpp"
#include "ovk/core/Disperse.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Misc.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Recv.hpp"
#include "ovk/core/RecvMap.hpp"
#include "ovk/core/Send.hpp"
#include "ovk/core/SendMap.hpp"

#include <mpi.h>

#include <map>
#include <utility>

namespace ovk {
namespace domain_internal {

namespace {

void UpdateSourceDestRanks(connectivity &Connectivity);

}

void UpdateConnectivitySourceDestRanks(domain &Domain) {

  // TODO: Try to find a way to do this that doesn't block for each grid pair (serializes work
  // that could be done in parallel)
  for (auto &MPair : Domain.LocalConnectivities_) {
    for (auto &NPair : MPair.second) {
      connectivity &Connectivity = NPair.second;
      UpdateSourceDestRanks(Connectivity);
    }
  }

}

namespace {

void UpdateSourceDestRanks(connectivity &Connectivity) {

  int NumDims;
  GetConnectivityDimension(Connectivity, NumDims);

  core::comm_view Comm = core::GetConnectivityComm(Connectivity);

  bool DonorGridIsLocal = RankHasConnectivityDonorSide(Connectivity);
  bool ReceiverGridIsLocal = RankHasConnectivityReceiverSide(Connectivity);

  const connectivity_d *Donors;
  const grid *DonorGrid;
  long long NumDonors;
  if (DonorGridIsLocal) {
    GetConnectivityDonorSide(Connectivity, Donors);
    GetConnectivityDonorSideGrid(*Donors, DonorGrid);
    GetConnectivityDonorSideCount(*Donors, NumDonors);
  }

  const connectivity_r *Receivers;
  const grid *ReceiverGrid;
  long long NumReceivers;
  if (ReceiverGridIsLocal) {
    GetConnectivityReceiverSide(Connectivity, Receivers);
    GetConnectivityReceiverSideGrid(*Receivers, ReceiverGrid);
    GetConnectivityReceiverSideCount(*Receivers, NumReceivers);
  }

  const connectivity::edits *Edits;
  core::GetConnectivityEdits(Connectivity, Edits);

  if (Edits->DonorDestinations_) {

    array<bool> DonorCommunicates;
    if (DonorGridIsLocal) {
      const cart &Cart = DonorGrid->Cart();
      const range &LocalRange = DonorGrid->LocalRange();
      DonorCommunicates.Resize({NumDonors});
      for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
        tuple<int> DonorLower = Cart.PeriodicAdjust({
          Donors->Extents_(0,0,iDonor),
          Donors->Extents_(0,1,iDonor),
          Donors->Extents_(0,2,iDonor)
        });
        DonorCommunicates(iDonor) = LocalRange.Contains(DonorLower);
      }
    }

    long long NumUnmapped = 0;
    if (DonorGridIsLocal) {
      for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
        if (DonorCommunicates(iDonor) && Donors->DestinationRanks_(iDonor) < 0) {
          ++NumUnmapped;
        }
      }
    }

    int GenerateDestinations = NumUnmapped > 0;
    MPI_Allreduce(MPI_IN_PLACE, &GenerateDestinations, 1, MPI_INT, MPI_MAX, Comm);

    if (GenerateDestinations) {

      const grid_info *ReceiverGridInfo;
      GetConnectivityReceiverGridInfo(Connectivity, ReceiverGridInfo);

      range ReceiverGridGlobalRange;
      GetGridInfoGlobalRange(*ReceiverGridInfo, ReceiverGridGlobalRange);

      range ReceiverGridLocalRange = MakeEmptyRange(NumDims);
      if (ReceiverGridIsLocal) {
        ReceiverGridLocalRange = ReceiverGrid->LocalRange();
      }

      core::partition_hash DestinationHash(NumDims, Comm, ReceiverGridGlobalRange,
        ReceiverGridLocalRange);

      array<int,2> Destinations;
      array<int> DestinationBinIndices;
      std::map<int, core::partition_hash::bin> Bins;

      if (DonorGridIsLocal) {
        Destinations.Resize({{MAX_DIMS,NumUnmapped}});
        DestinationBinIndices.Resize({NumUnmapped});
        long long iUnmapped = 0;
        for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
          if (DonorCommunicates(iDonor) && Donors->DestinationRanks_(iDonor) < 0) {
            for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
              Destinations(iDim,iUnmapped) = Donors->Destinations_(iDim,iDonor);
            }
            ++iUnmapped;
          }
        }
        DestinationHash.MapToBins(Destinations, DestinationBinIndices);
        for (long long iUnmapped = 0; iUnmapped < NumUnmapped; ++iUnmapped) {
          int BinIndex = DestinationBinIndices(iUnmapped);
          auto Iter = Bins.lower_bound(BinIndex);
          if (Iter == Bins.end() || Iter->first > BinIndex) {
            Bins.emplace_hint(Iter, BinIndex, core::partition_hash::bin());
          }
        }
      }

      DestinationHash.RetrieveBins(Bins);

      if (DonorGridIsLocal) {
        array<int> DestinationRanks({NumUnmapped});
        DestinationHash.FindPartitions(Bins, Destinations, DestinationBinIndices,
          DestinationRanks);
        connectivity_d *DonorsEdit;
        EditConnectivityDonorSideLocal(Connectivity, DonorsEdit);
        int *DestinationRanksEdit;
        EditDonorDestinationRanks(*DonorsEdit, DestinationRanksEdit);
        long long iUnmapped = 0;
        for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
          if (DonorCommunicates(iDonor) && Donors->DestinationRanks_(iDonor) < 0) {
            DestinationRanksEdit[iDonor] = DestinationRanks(iUnmapped);
            ++iUnmapped;
          }
        }
        ReleaseDonorDestinationRanks(*DonorsEdit, DestinationRanksEdit);
        ReleaseConnectivityDonorSideLocal(Connectivity, DonorsEdit);
      } else {
        EditConnectivityDonorSideRemote(Connectivity);
        ReleaseConnectivityDonorSideRemote(Connectivity);
      }

    }

  }

  if (Edits->ReceiverSources_) {

    long long NumUnmapped = 0;
    if (ReceiverGridIsLocal) {
      for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
        if (Receivers->SourceRanks_(iReceiver) < 0) {
          ++NumUnmapped;
        }
      }
    }

    int GenerateSources = NumUnmapped > 0;
    MPI_Allreduce(MPI_IN_PLACE, &GenerateSources, 1, MPI_INT, MPI_MAX, Comm);

    if (GenerateSources) {

      const grid_info *DonorGridInfo;
      GetConnectivityDonorGridInfo(Connectivity, DonorGridInfo);

      range DonorGridGlobalRange;
      GetGridInfoGlobalRange(*DonorGridInfo, DonorGridGlobalRange);

      range DonorGridLocalRange = MakeEmptyRange(NumDims);
      if (DonorGridIsLocal) {
        DonorGridLocalRange = DonorGrid->LocalRange();
      }

      core::partition_hash SourceHash(NumDims, Comm, DonorGridGlobalRange, DonorGridLocalRange);

      array<int,2> Sources;
      array<int> SourceBinIndices;
      std::map<int, core::partition_hash::bin> Bins;

      if (ReceiverGridIsLocal) {
        Sources.Resize({{MAX_DIMS,NumUnmapped}});
        SourceBinIndices.Resize({NumUnmapped});
        long long iUnmapped = 0;
        for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
          if (Receivers->SourceRanks_(iReceiver) < 0) {
            for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
              Sources(iDim,iUnmapped) = Receivers->Sources_(iDim,iReceiver);
            }
            ++iUnmapped;
          }
        }
        SourceHash.MapToBins(Sources, SourceBinIndices);
        for (long long iUnmapped = 0; iUnmapped < NumUnmapped; ++iUnmapped) {
          int BinIndex = SourceBinIndices(iUnmapped);
          auto Iter = Bins.lower_bound(BinIndex);
          if (Iter == Bins.end() || Iter->first > BinIndex) {
            Bins.emplace_hint(Iter, BinIndex, core::partition_hash::bin());
          }
        }
      }

      SourceHash.RetrieveBins(Bins);

      if (ReceiverGridIsLocal) {
        array<int> SourceRanks({NumUnmapped});
        SourceHash.FindPartitions(Bins, Sources, SourceBinIndices,
          SourceRanks);
        connectivity_r *ReceiversEdit;
        EditConnectivityReceiverSideLocal(Connectivity, ReceiversEdit);
        int *SourceRanksEdit;
        EditReceiverSourceRanks(*ReceiversEdit, SourceRanksEdit);
        long long iUnmapped = 0;
        for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
          if (Receivers->SourceRanks_(iReceiver) < 0) {
            SourceRanksEdit[iReceiver] = SourceRanks(iUnmapped);
            ++iUnmapped;
          }
        }
        ReleaseReceiverSourceRanks(*ReceiversEdit, SourceRanksEdit);
        ReleaseConnectivityReceiverSideLocal(Connectivity, ReceiversEdit);
      } else {
        EditConnectivityReceiverSideRemote(Connectivity);
        ReleaseConnectivityReceiverSideRemote(Connectivity);
      }

    }

  }

}

}

void UpdateExchanges(domain &Domain) {

  MPI_Barrier(Domain.Comm_);

//   for (auto &MPair : Domain.LocalExchanges_) {
//     for (auto &NPair : MPair.second) {
//       exchange &Exchange = NPair.second;
//       core::UpdateExchange(Exchange);
//     }
//   }

  int CollectTime = core::GetProfilerTimerID(Domain.Profiler_, "Collect");
  int SendRecvTime = core::GetProfilerTimerID(Domain.Profiler_, "SendRecv");
  int DisperseTime = core::GetProfilerTimerID(Domain.Profiler_, "Disperse");

  for (auto &MPair : Domain.LocalConnectivities_) {
    int DonorGridID = MPair.first;
    for (auto &NPair : MPair.second) {
      int ReceiverGridID = NPair.first;
      const connectivity &Connectivity = NPair.second;
      if (RankHasConnectivityDonorSide(Connectivity)) {
        const connectivity::edits *Edits;
        core::GetConnectivityEdits(Connectivity, Edits);
        core::StartProfile(Domain.Profiler_, CollectTime);
        auto CollectDataRowIter = Domain.CollectData_.find(DonorGridID);
        if (CollectDataRowIter != Domain.CollectData_.end()) {
          std::map<int, domain::collect_data> &CollectDataRow = CollectDataRowIter->second;
          CollectDataRow.erase(ReceiverGridID);
          if (CollectDataRow.empty()) {
            Domain.CollectData_.erase(CollectDataRowIter);
          }
        }
        core::EndProfile(Domain.Profiler_, CollectTime);
        core::StartProfile(Domain.Profiler_, SendRecvTime);
        auto SendDataRowIter = Domain.SendData_.find(DonorGridID);
        if (SendDataRowIter != Domain.SendData_.end()) {
          std::map<int, domain::send_data> &SendDataRow = SendDataRowIter->second;
          SendDataRow.erase(ReceiverGridID);
          if (SendDataRow.empty()) {
            Domain.SendData_.erase(SendDataRowIter);
          }
        }
        auto RecvDataRowIter = Domain.RecvData_.find(DonorGridID);
        if (RecvDataRowIter != Domain.RecvData_.end()) {
          std::map<int, domain::recv_data> &RecvDataRow = RecvDataRowIter->second;
          RecvDataRow.erase(ReceiverGridID);
          if (RecvDataRow.empty()) {
            Domain.RecvData_.erase(RecvDataRowIter);
          }
        }
        core::EndProfile(Domain.Profiler_, SendRecvTime);
        core::StartProfile(Domain.Profiler_, DisperseTime);
        auto DisperseDataRowIter = Domain.DisperseData_.find(DonorGridID);
        if (DisperseDataRowIter != Domain.DisperseData_.end()) {
          std::map<int, domain::disperse_data> &DisperseDataRow = DisperseDataRowIter->second;
          DisperseDataRow.erase(ReceiverGridID);
          if (DisperseDataRow.empty()) {
            Domain.DisperseData_.erase(DisperseDataRowIter);
          }
        }
        core::EndProfile(Domain.Profiler_, DisperseTime);
      }
    }
  }

  MPI_Barrier(Domain.Comm_);

}

}}
