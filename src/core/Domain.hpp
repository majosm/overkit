// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DOMAIN_HPP_INCLUDED
#define OVK_CORE_DOMAIN_HPP_INCLUDED

#include <ovk/core/AssemblyOptions.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Connectivity.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/ErrorHandler.hpp>
#include <ovk/core/Exchange.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Logger.hpp>

#include <mpi.h>

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

  mutable core::logger *Logger_;
  mutable core::error_handler *ErrorHandler_;
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
const comm &GetDomainComm(const domain &Domain);
logger &GetDomainLogger(const domain &Domain);
error_handler &GetDomainErrorHandler(const domain &Domain);
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

void Collect(const domain &Domain, int DonorGridID, int ReceiverGridID, data_type ValueType,
  int Count, collect_op CollectOp, const range &GridValuesRange, array_layout GridValuesLayout,
  const void * const *GridValues, void **DonorValues);

void Send(const domain &Domain, int DonorGridID, int ReceiverGridID, data_type ValueType, int Count,
  const void * const *DonorValues, int Tag, request &Request);

void Receive(const domain &Domain, int DonorGridID, int ReceiverGridID, data_type ValueType,
  int Count, void **ReceiverValues, int Tag, request &Request);

void Disperse(const domain &Domain, int DonorGridID, int ReceiverGridID, data_type ValueType,
  int Count, disperse_op DisperseOp, const void * const *ReceiverValues, const range
  &GridValuesRange, array_layout GridValuesLayout, void **GridValues);

void CreateDomainParams(domain_params &Params, int NumDims);
void DestroyDomainParams(domain_params &Params);

void GetDomainParamName(const domain_params &Params, std::string &Name);
void SetDomainParamName(domain_params &Params, std::string Name);
void GetDomainParamDimension(const domain_params &Params, int &NumDims);
void GetDomainParamComm(const domain_params &Params, MPI_Comm &Comm);
void SetDomainParamComm(domain_params &Params, MPI_Comm Comm);

}

#endif
