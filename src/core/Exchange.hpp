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
  int NumDims_;
  core::comm Comm_;

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
  error_handler &ErrorHandler);
void DestroyExchange(exchange &Exchange);

void CreateExchangeInfo(exchange_info &Info, const exchange *Exchange, const comm &Comm);
void DestroyExchangeInfo(exchange_info &Info);

}

bool RankHasExchangeDonorSide(const exchange &Exchange);
bool RankHasExchangeReceiverSide(const exchange &Exchange);

namespace core {

void UpdateExchange(exchange &Exchange);

void Collect(const exchange &Exchange, data_type ValueType, int Count, collect_op CollectOp,
  const range &GridValuesRange, array_layout GridValuesLayout, const void * const *GridValues,
  void **DonorValues);
void Send(const exchange &Exchange, data_type ValueType, int Count, const void * const *DonorValues,
  int Tag, request &Request);
void Receive(const exchange &Exchange, data_type ValueType, int Count, void **ReceiverValues,
  int Tag, request &Request);
void Disperse(const exchange &Exchange, data_type ValueType, int Count, disperse_op DisperseOp,
  const void * const *ReceiverValues, const range &GridValuesRange, array_layout GridValuesLayout,
  void **GridValues);

}

void GetExchangeInfoDonorGridID(const exchange_info &Info, int &DonorGridID);
void GetExchangeInfoReceiverGridID(const exchange_info &Info, int &ReceiverGridID);
void GetExchangeInfoName(const exchange_info &Info, std::string &Name);
void GetExchangeInfoDimension(const exchange_info &Info, int &NumDims);
void GetExchangeInfoRootRank(const exchange_info &Info, int &RootRank);
void GetExchangeInfoIsLocal(const exchange_info &Info, bool &IsLocal);

}

#endif
