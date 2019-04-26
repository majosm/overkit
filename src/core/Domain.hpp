// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DOMAIN_HPP_INCLUDED
#define OVK_CORE_DOMAIN_HPP_INCLUDED

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/AssemblyOptions.hpp>
#include <ovk/core/Collect.hpp>
#include <ovk/core/CollectMap.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Connectivity.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Disperse.hpp>
#include <ovk/core/ErrorHandler.hpp>
#include <ovk/core/Exchange.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Logger.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Recv.hpp>
#include <ovk/core/RecvMap.hpp>
#include <ovk/core/Send.hpp>
#include <ovk/core/SendMap.hpp>

#include <mpi.h>

#include <list>
#include <map>
#include <string>

namespace ovk {

struct domain_params {
  std::string Name_;
  int NumDims_;
  MPI_Comm Comm_;
};

struct domain {

  struct grid_info : ovk::grid_info {
    int EditRefCount_;
  };

  struct connectivity_info : ovk::connectivity_info {
    int EditRefCount_;
  };

  struct collect_data {
    core::collect_map Map;
    std::map<int, core::collect> Collects;
  };

  struct send_data {
    core::send_map Map;
    std::map<int, core::send> Sends;
  };

  struct recv_data {
    core::recv_map Map;
    std::map<int, core::recv> Recvs;
  };

  struct disperse_data {
    std::map<int, core::disperse> Disperses;
  };

  mutable core::logger *Logger_;
  mutable core::error_handler *ErrorHandler_;
  mutable core::profiler Profiler_;
  std::string Name_;
  int NumDims_;
  core::comm Comm_;
  domain_config Config_;
  int NumGrids_;
  std::map<int, grid_info> GridInfo_;
  int AllGridsEditRefCount_;
  std::map<int, grid> LocalGrids_;
  std::map<int, std::map<int, connectivity_info>> ConnectivityInfo_;
  int AllConnectivitiesEditRefCount_;
  std::map<int, std::map<int, connectivity>> LocalConnectivities_;
  std::map<int, std::map<int, exchange_info>> ExchangeInfo_;
  std::map<int, std::map<int, exchange>> LocalExchanges_;
  std::map<int, std::map<int, collect_data>> CollectData_;
  std::map<int, std::map<int, send_data>> SendData_;
  std::map<int, std::map<int, recv_data>> RecvData_;
  std::map<int, std::map<int, disperse_data>> DisperseData_;

};

namespace core {
void CreateDomain(domain &Domain, const domain_params &Params, logger &Logger, error_handler
  &ErrorHandler);
void DestroyDomain(domain &Domain);
}

void GetDomainName(const domain &Domain, std::string &Name);
void GetDomainDimension(const domain &Domain, int &NumDims);
void GetDomainComm(const domain &Domain, MPI_Comm &Comm);
void GetDomainCommSize(const domain &Domain, int &CommSize);
void GetDomainCommRank(const domain &Domain, int &CommRank);

void ConfigureDomain(domain &Domain, domain_config Config);
void GetDomainConfiguration(const domain &Domain, domain_config &Config);

void GetDomainGridCount(const domain &Domain, int &NumGrids);
void GetDomainGridIDs(const domain &Domain, int *GridIDs);
void GetNextAvailableGridID(const domain &Domain, int &GridID);

namespace core {
comm_view GetDomainComm(const domain &Domain);
logger &GetDomainLogger(const domain &Domain);
error_handler &GetDomainErrorHandler(const domain &Domain);
profiler &GetDomainProfiler(const domain &Domain);
}

void CreateGridLocal(domain &Domain, int GridID, const grid_params &Params);
void CreateGridRemote(domain &Domain, int GridID);
void DestroyGrid(domain &Domain, int GridID);

bool GridExists(const domain &Domain, int GridID);
void GetGridInfo(const domain &Domain, int GridID, const grid_info *&GridInfo);
bool RankHasGrid(const domain &Domain, int GridID);
void GetGrid(const domain &Domain, int GridID, const grid *&Grid);
void EditGridLocal(domain &Domain, int GridID, grid *&Grid);
void EditGridRemote(domain &Domain, int GridID);
void ReleaseGridLocal(domain &Domain, int GridID, grid *&Grid);
void ReleaseGridRemote(domain &Domain, int GridID);

bool ConnectivityExists(const domain &Domain, int DonorGridID, int ReceiverGridID);
void GetConnectivityInfo(const domain &Domain, int DonorGridID, int ReceiverGridID, const
  connectivity_info *&ConnectivityInfo);
bool RankHasConnectivity(const domain &Domain, int DonorGridID, int ReceiverGridID);
void GetConnectivity(const domain &Domain, int DonorGridID, int ReceiverGridID,
  const connectivity *&Connectivity);
void EditConnectivityLocal(domain &Domain, int DonorGridID, int ReceiverGridID,
  connectivity *&Connectivity);
void EditConnectivityRemote(domain &Domain, int DonorGridID, int ReceiverGridID);
void ReleaseConnectivityLocal(domain &Domain, int DonorGridID, int ReceiverGridID,
  connectivity *&Connectivity);
void ReleaseConnectivityRemote(domain &Domain, int DonorGridID, int ReceiverGridID);

bool ExchangeExists(const domain &Domain, int DonorGridID, int ReceiverGridID);
void GetExchangeInfo(const domain &Domain, int DonorGridID, int ReceiverGridID, const exchange_info
  *&ExchangeInfo);
bool RankHasExchange(const domain &Domain, int DonorGridID, int ReceiverGridID);
void GetExchange(const domain &Domain, int DonorGridID, int ReceiverGridID,
  const exchange *&Exchange);

void GetLocalDonorCount(const domain &Domain, int DonorGridID, int ReceiverGridID,
  long long &NumDonors);
void GetLocalReceiverCount(const domain &Domain, int DonorGridID, int ReceiverGridID,
  long long &NumReceivers);

void Assemble(domain &Domain, const assembly_options &Options);

void CreateCollect(domain &Domain, int DonorGridID, int ReceiverGridID, int CollectID, collect_op
  CollectOp, data_type ValueType, int Count, const range &GridValuesRange, array_layout
  GridValuesLayout);
void DestroyCollect(domain &Domain, int DonorGridID, int ReceiverGridID, int CollectID);
void Collect(domain &Domain, int DonorGridID, int ReceiverGridID, int CollectID, const void * const
  *GridValues, void **DonorValues);
bool CollectExists(const domain &Domain, int DonorGridID, int ReceiverGridID, int CollectID);
void GetNextAvailableCollectID(const domain &Domain, int DonorGridID, int ReceiverGridID, int
  &CollectID);

void CreateSend(domain &Domain, int DonorGridID, int ReceiverGridID, int SendID, data_type
  ValueType, int Count, int Tag);
void DestroySend(domain &Domain, int DonorGridID, int ReceiverGridID, int SendID);
request Send(domain &Domain, int DonorGridID, int ReceiverGridID, int SendID, const void * const
  *DonorValues);
bool SendExists(const domain &Domain, int DonorGridID, int ReceiverGridID, int SendID);
void GetNextAvailableSendID(const domain &Domain, int DonorGridID, int ReceiverGridID, int &SendID);

void CreateReceive(domain &Domain, int DonorGridID, int ReceiverGridID, int RecvID, data_type
  ValueType, int Count, int Tag);
void DestroyReceive(domain &Domain, int DonorGridID, int ReceiverGridID, int RecvID);
request Receive(domain &Domain, int DonorGridID, int ReceiverGridID, int RecvID, void
  **ReceiverValues);
bool ReceiveExists(const domain &Domain, int DonorGridID, int ReceiverGridID, int RecvID);
void GetNextAvailableReceiveID(const domain &Domain, int DonorGridID, int ReceiverGridID, int
  &RecvID);

void Wait(const domain &Domain, request &Request);
void WaitAll(const domain &Domain, array_view<request> Requests);
void WaitAny(const domain &Domain, array_view<request> Requests, int &Index);
// Needed for C API
void WaitAll(const domain &Domain, array_view<request *> Requests);
void WaitAny(const domain &Domain, array_view<request *> Requests, int &Index);

void CreateDisperse(domain &Domain, int DonorGridID, int ReceiverGridID, int DisperseID, disperse_op
  DisperseOp, data_type ValueType, int Count, const range &GridValuesRange, array_layout
  GridValuesLayout);
void DestroyDisperse(domain &Domain, int DonorGridID, int ReceiverGridID, int DisperseID);
void Disperse(domain &Domain, int DonorGridID, int ReceiverGridID, int DisperseID, const void *
  const *ReceiverValues, void **GridValues);
bool DisperseExists(const domain &Domain, int DonorGridID, int ReceiverGridID, int DisperseID);
void GetNextAvailableDisperseID(const domain &Domain, int DonorGridID, int ReceiverGridID, int
  &DisperseID);

void CreateDomainParams(domain_params &Params, int NumDims);
void DestroyDomainParams(domain_params &Params);

void GetDomainParamName(const domain_params &Params, std::string &Name);
void SetDomainParamName(domain_params &Params, std::string Name);
void GetDomainParamDimension(const domain_params &Params, int &NumDims);
void GetDomainParamComm(const domain_params &Params, MPI_Comm &Comm);
void SetDomainParamComm(domain_params &Params, MPI_Comm Comm);

}

#endif
