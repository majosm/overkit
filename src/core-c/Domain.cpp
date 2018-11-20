// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Domain.h"

#include "ovk/core-c/AssemblyOptions.h"
#include "ovk/core-c/Connectivity.h"
#include "ovk/core-c/Constants.h"
#include "ovk/core-c/Exchange.h"
#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core-c/Request.h"
#include "ovk/core/AssemblyOptions.hpp"
#include "ovk/core/Connectivity.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Request.hpp"

#include <mpi.h>

#include <cstring>
#include <string>

extern "C" {

void ovkGetDomainName(const ovk_domain *Domain, char *Name) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);

  std::string NameCPP;
  ovk::GetDomainName(DomainCPP, NameCPP);

  std::strcpy(Name, NameCPP.c_str());

}

void ovkGetDomainDimension(const ovk_domain *Domain, int *NumDims) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::GetDomainDimension(DomainCPP, *NumDims);

}

void ovkGetDomainComm(const ovk_domain *Domain, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::GetDomainComm(DomainCPP, *Comm);

}

void ovkGetDomainCommSize(const ovk_domain *Domain, int *CommSize) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::GetDomainCommSize(DomainCPP, *CommSize);

}

void ovkGetDomainCommRank(const ovk_domain *Domain, int *CommRank) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::GetDomainCommRank(DomainCPP, *CommRank);

}

void ovkConfigureDomain(ovk_domain *Domain, ovk_domain_config Config) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  ovk::ConfigureDomain(DomainCPP, ovk::domain_config(Config));

}

void ovkGetDomainConfiguration(const ovk_domain *Domain, ovk_domain_config *Config) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Config, "Invalid config pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);

  ovk::domain_config ConfigCPP;
  ovk::GetDomainConfiguration(DomainCPP, ConfigCPP);

  *Config = ovk_domain_config(ConfigCPP);

}

void ovkGetDomainGridCount(const ovk_domain *Domain, int *NumGrids) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(NumGrids, "Invalid num grids pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::GetDomainGridCount(DomainCPP, *NumGrids);

}

void ovkGetDomainGridIDs(const ovk_domain *Domain, int *GridIDs) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridIDs, "Invalid grid IDs pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::GetDomainGridIDs(DomainCPP, GridIDs);

}

void ovkGetNextAvailableGridID(const ovk_domain *Domain, int *GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::GetNextAvailableGridID(DomainCPP, *GridID);

}

void ovkCreateGridLocal(ovk_domain *Domain, int GridID, const ovk_grid_params *Params) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  auto &ParamsCPP = *reinterpret_cast<const ovk::grid_params *>(Params);
  ovk::CreateGridLocal(DomainCPP, GridID, ParamsCPP);

}

void ovkCreateGridRemote(ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  ovk::CreateGridRemote(DomainCPP, GridID);

}

void ovkDestroyGrid(ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  ovk::DestroyGrid(DomainCPP, GridID);

}

bool ovkGridExists(const ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  return ovk::GridExists(DomainCPP, GridID);

}

void ovkGetGridInfo(const ovk_domain *Domain, int GridID, const ovk_grid_info **GridInfo) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridInfo, "Invalid grid info pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  
  const ovk::grid_info *GridInfoCPPPtr;
  ovk::GetGridInfo(DomainCPP, GridID, GridInfoCPPPtr);

  *GridInfo = reinterpret_cast<const ovk_grid_info *>(GridInfoCPPPtr);

}

bool ovkRankHasGrid(const ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  return ovk::RankHasGrid(DomainCPP, GridID);

}

void ovkGetGrid(const ovk_domain *Domain, int GridID, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);

  const ovk::grid *GridCPPPtr;
  ovk::GetGrid(DomainCPP, GridID, GridCPPPtr);

  *Grid = reinterpret_cast<const ovk_grid *>(GridCPPPtr);

}

void ovkEditGridLocal(ovk_domain *Domain, int GridID, ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);

  ovk::grid *GridCPPPtr;
  ovk::EditGridLocal(DomainCPP, GridID, GridCPPPtr);

  *Grid = reinterpret_cast<ovk_grid *>(GridCPPPtr);

}

void ovkEditGridRemote(ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  ovk::EditGridRemote(DomainCPP, GridID);

}

void ovkReleaseGridLocal(ovk_domain *Domain, int GridID, ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid grid pointer.");
  OVK_DEBUG_ASSERT(*Grid, "Invalid grid pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);

  auto GridCPPPtr = reinterpret_cast<ovk::grid *>(*Grid);
  ovk::ReleaseGridLocal(DomainCPP, GridID, GridCPPPtr);

  *Grid = nullptr;

}

void ovkReleaseGridRemote(ovk_domain *Domain, int GridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  ovk::ReleaseGridRemote(DomainCPP, GridID);

}

bool ovkConnectivityExists(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  return ovk::ConnectivityExists(DomainCPP, DonorGridID, ReceiverGridID);

}

void ovkGetConnectivityInfo(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID, const
  ovk_connectivity_info **ConnectivityInfo) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(ConnectivityInfo, "Invalid connectivity info pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  
  const ovk::connectivity_info *ConnectivityInfoCPPPtr;
  ovk::GetConnectivityInfo(DomainCPP, DonorGridID, ReceiverGridID, ConnectivityInfoCPPPtr);

  *ConnectivityInfo = reinterpret_cast<const ovk_connectivity_info *>(ConnectivityInfoCPPPtr);

}

bool ovkRankHasConnectivity(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  return ovk::RankHasConnectivity(DomainCPP, DonorGridID, ReceiverGridID);

}

void ovkGetConnectivity(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  const ovk_connectivity **Connectivity) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);

  const ovk::connectivity *ConnectivityCPPPtr;
  ovk::GetConnectivity(DomainCPP, DonorGridID, ReceiverGridID, ConnectivityCPPPtr);

  *Connectivity = reinterpret_cast<const ovk_connectivity *>(ConnectivityCPPPtr);

}

void ovkEditConnectivityLocal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);

  ovk::connectivity *ConnectivityCPPPtr;
  ovk::EditConnectivityLocal(DomainCPP, DonorGridID, ReceiverGridID, ConnectivityCPPPtr);

  *Connectivity = reinterpret_cast<ovk_connectivity *>(ConnectivityCPPPtr);

}

void ovkEditConnectivityRemote(ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  ovk::EditConnectivityRemote(DomainCPP, DonorGridID, ReceiverGridID);

}

void ovkReleaseConnectivityLocal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Connectivity, "Invalid connectivity pointer.");
  OVK_DEBUG_ASSERT(*Connectivity, "Invalid connectivity pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);

  auto ConnectivityCPPPtr = reinterpret_cast<ovk::connectivity *>(*Connectivity);
  ovk::ReleaseConnectivityLocal(DomainCPP, DonorGridID, ReceiverGridID, ConnectivityCPPPtr);

  *Connectivity = nullptr;

}

void ovkReleaseConnectivityRemote(ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  ovk::ReleaseConnectivityRemote(DomainCPP, DonorGridID, ReceiverGridID);

}

bool ovkExchangeExists(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  return ovk::ExchangeExists(DomainCPP, DonorGridID, ReceiverGridID);

}

void ovkGetExchangeInfo(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID, const
  ovk_exchange_info **ExchangeInfo) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(ExchangeInfo, "Invalid exchange info pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  
  const ovk::exchange_info *ExchangeInfoCPPPtr;
  ovk::GetExchangeInfo(DomainCPP, DonorGridID, ReceiverGridID, ExchangeInfoCPPPtr);

  *ExchangeInfo = reinterpret_cast<const ovk_exchange_info *>(ExchangeInfoCPPPtr);

}

bool ovkRankHasExchange(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  return ovk::RankHasExchange(DomainCPP, DonorGridID, ReceiverGridID);

}

void ovkGetExchange(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  const ovk_exchange **Exchange) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Exchange, "Invalid exchange pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);

  const ovk::exchange *ExchangeCPPPtr;
  ovk::GetExchange(DomainCPP, DonorGridID, ReceiverGridID, ExchangeCPPPtr);

  *Exchange = reinterpret_cast<const ovk_exchange *>(ExchangeCPPPtr);

}

void ovkGetLocalDonorCount(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  long long *NumDonors) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(NumDonors, "Invalid num donors pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::GetLocalDonorCount(DomainCPP, DonorGridID, ReceiverGridID, *NumDonors);

}

void ovkGetLocalReceiverCount(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  long long *NumReceivers) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(NumReceivers, "Invalid num receivers pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::GetLocalReceiverCount(DomainCPP, DonorGridID, ReceiverGridID, *NumReceivers);

}

void ovkAssemble(ovk_domain *Domain, const ovk_assembly_options *Options) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Options, "Invalid domain pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  auto &OptionsCPP = *reinterpret_cast<const ovk::assembly_options *>(Options);
  ovk::Assemble(DomainCPP, OptionsCPP);

}

void ovkCollect(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_data_type ValueType, int Count, ovk_collect_op CollectOp, const ovk_range *GridValuesRange,
  ovk_array_layout GridValuesLayout, const void **GridValues, void **DonorValues) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridValuesRange, "Invalid grid data range pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::Collect(DomainCPP, DonorGridID, ReceiverGridID, ovk::data_type(ValueType), Count,
    ovk::collect_op(CollectOp), *GridValuesRange, ovk::array_layout(GridValuesLayout), GridValues,
    DonorValues);

}

void ovkSend(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID, ovk_data_type ValueType,
  int Count, const void **DonorValues, int Tag, ovk_request **Request) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Request, "Invalid request pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  auto RequestCPPPtr = new ovk::request();

  ovk::Send(DomainCPP, DonorGridID, ReceiverGridID, ovk::data_type(ValueType), Count, DonorValues,
    Tag, *RequestCPPPtr);

  *Request = reinterpret_cast<ovk_request *>(RequestCPPPtr);

}

void ovkReceive(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_data_type ValueType, int Count, void **ReceiverValues, int Tag, ovk_request **Request) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Request, "Invalid request pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  auto RequestCPPPtr = new ovk::request();

  ovk::Receive(DomainCPP, DonorGridID, ReceiverGridID, ovk::data_type(ValueType), Count,
    ReceiverValues, Tag, *RequestCPPPtr);

  *Request = reinterpret_cast<ovk_request *>(RequestCPPPtr);

}

void ovkDisperse(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_data_type ValueType, int Count, ovk_disperse_op DisperseOp, const void **ReceiverValues,
  const ovk_range *GridValuesRange, ovk_array_layout GridValuesLayout, void **GridValues) {

  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(GridValuesRange, "Invalid grid data range pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::Disperse(DomainCPP, DonorGridID, ReceiverGridID, ovk::data_type(ValueType), Count,
    ovk::disperse_op(DisperseOp), ReceiverValues, *GridValuesRange,
    ovk::array_layout(GridValuesLayout), GridValues);

}

void ovkCreateDomainParams(ovk_domain_params **Params, int NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::domain_params();

  ovk::CreateDomainParams(*ParamsCPPPtr, NumDims);

  *Params = reinterpret_cast<ovk_domain_params *>(ParamsCPPPtr);

}

void ovkDestroyDomainParams(ovk_domain_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::domain_params *>(*Params);

  ovk::DestroyDomainParams(*ParamsCPPPtr);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetDomainParamName(const ovk_domain_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::domain_params *>(Params);

  std::string NameCPP;
  ovk::GetDomainParamName(ParamsCPP, NameCPP);

  std::strcpy(Name, NameCPP.c_str());

}

void ovkSetDomainParamName(ovk_domain_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::domain_params *>(Params);
  ovk::SetDomainParamName(ParamsCPP, Name);

}

void ovkGetDomainParamDimension(const ovk_domain_params *Params, int *NumDims) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::domain_params *>(Params);
  ovk::GetDomainParamDimension(ParamsCPP, *NumDims);

}

void ovkGetDomainParamComm(const ovk_domain_params *Params, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::domain_params *>(Params);
  ovk::GetDomainParamComm(ParamsCPP, *Comm);

}

void ovkSetDomainParamComm(ovk_domain_params *Params, MPI_Comm Comm) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::domain_params *>(Params);
  ovk::SetDomainParamComm(ParamsCPP, Comm);

}

}
