// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_EXCHANGE_HPP_INCLUDED
#define OVK_CORE_EXCHANGE_HPP_INCLUDED

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
#include <vector>

namespace ovk {

struct exchange {

  struct collect_send {
    int Rank;
    long long NumPoints;
    int *Points[MAX_DIMS];
    std::vector<int> PointsData;
  };

  struct collect_recv {
    int Rank;
    long long NumPoints;
    int *Points[MAX_DIMS];
    std::vector<int> PointsData;
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
  MPI_Comm Comm_;
  int CommSize_;
  int CommRank_;

  std::vector<long long> DonorsSorted_;
  std::vector<long long> ReceiversSorted_;
  std::vector<int> DonorDestRanks_;
  std::vector<int> ReceiverSourceRanks_;
  core::partition_hash SourceHash_;
  core::partition_hash DestinationHash_;

  std::vector<collect_send> CollectSends_;
  std::vector<collect_recv> CollectRecvs_;
  std::vector<int> NumRemoteDonorPoints_;
  long long **RemoteDonorPoints_;
  std::vector<long long *> RemoteDonorPointsPtrs_;
  std::vector<long long> RemoteDonorPointsData_;
  int **RemoteDonorPointCollectRecvs_;
  std::vector<int *> RemoteDonorPointCollectRecvsPtrs_;
  std::vector<int> RemoteDonorPointCollectRecvsData_;
  long long **RemoteDonorPointCollectRecvBufferIndices_;
  std::vector<long long *> RemoteDonorPointCollectRecvBufferIndicesPtrs_;
  std::vector<long long> RemoteDonorPointCollectRecvBufferIndicesData_;

  std::vector<send> Sends_;
  std::vector<recv> Recvs_;
  std::vector<int> DonorSendIndices_;
  std::vector<int> ReceiverRecvIndices_;

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

void CreateExchangeInfo(exchange_info &Info, const exchange *Exchange, MPI_Comm Comm, int CommRank);
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
