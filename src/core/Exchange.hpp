// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_EXCHANGE_HPP_INCLUDED
#define OVK_CORE_EXCHANGE_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Connectivity.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/ErrorHandler.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Logger.hpp>
#include <ovk/core/PartitionHash.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Request.hpp>

#include <mpi.h>

#include <string>

namespace ovk {

struct exchange {

  struct collect_send {
    int Rank;
    long long NumPoints;
    array<int,2> Points;
  };

  struct collect_recv {
    int Rank;
    long long NumPoints;
  };

  struct send {
    int Rank;
    long long Count;
  };

  struct recv {
    int Rank;
    long long Count;
  };

  const connectivity *Connectivity_;
  mutable core::logger *Logger_;
  mutable core::error_handler *ErrorHandler_;
  mutable core::profiler *Profiler_;
  int NumDims_;
  core::comm_view Comm_;

  array<long long> DonorsSorted_;
  array<long long> ReceiversSorted_;
  array<int> DonorDestRanks_;
  array<int> ReceiverSourceRanks_;
  core::partition_hash SourceHash_;
  core::partition_hash DestinationHash_;

  array<collect_send> CollectSends_;
  array<collect_recv> CollectRecvs_;
  array<int> NumRemoteDonorPoints_;
  array<long long *> RemoteDonorPoints_;
  array<long long> RemoteDonorPointsData_;
  array<int *> RemoteDonorPointCollectRecvs_;
  array<int> RemoteDonorPointCollectRecvsData_;
  array<long long *> RemoteDonorPointCollectRecvBufferIndices_;
  array<long long> RemoteDonorPointCollectRecvBufferIndicesData_;

  array<send> Sends_;
  array<recv> Recvs_;
  array<int> DonorSendIndices_;
  array<int> ReceiverRecvIndices_;

};

struct exchange_info {
  int DonorGridID_;
  int ReceiverGridID_;
  std::string Name_;
  int NumDims_;
  int RootRank_;
  bool IsLocal_;
};

namespace core {

void CreateExchange(exchange &Exchange, const connectivity &Connectivity, logger &Logger,
  error_handler &ErrorHandler, profiler &Profiler);
void DestroyExchange(exchange &Exchange);

void CreateExchangeInfo(exchange_info &Info, const exchange *Exchange, comm_view Comm);
void DestroyExchangeInfo(exchange_info &Info);

}

bool RankHasExchangeDonorSide(const exchange &Exchange);
bool RankHasExchangeReceiverSide(const exchange &Exchange);

namespace core {
void UpdateExchange(exchange &Exchange);
}

void GetExchangeInfoDonorGridID(const exchange_info &Info, int &DonorGridID);
void GetExchangeInfoReceiverGridID(const exchange_info &Info, int &ReceiverGridID);
void GetExchangeInfoName(const exchange_info &Info, std::string &Name);
void GetExchangeInfoDimension(const exchange_info &Info, int &NumDims);
void GetExchangeInfoRootRank(const exchange_info &Info, int &RootRank);
void GetExchangeInfoIsLocal(const exchange_info &Info, bool &IsLocal);

}

#endif
