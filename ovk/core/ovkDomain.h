// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_DOMAIN_INCLUDED
#define OVK_CORE_PUBLIC_DOMAIN_INCLUDED

#include <ovk/core/ovkConnectivity.h>
#include <ovk/core/ovkGlobal.h>
#include <ovk/core/ovkGrid.h>
#include <ovk/core/ovkExchange.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_domain_params;
typedef struct ovk_domain_params ovk_domain_params;

struct ovk_domain;
typedef struct ovk_domain ovk_domain;

struct ovk_domain_properties;
typedef struct ovk_domain_properties ovk_domain_properties;

void ovkConfigureDomain(ovk_domain *Domain, ovk_domain_config Config);
void ovkGetDomainConfiguration(ovk_domain *Domain, ovk_domain_config *Config);

void ovkGetDomainProperties(const ovk_domain *Domain, const ovk_domain_properties **Properties);
void ovkEditDomainProperties(ovk_domain *Domain, ovk_domain_properties **Properties);
void ovkReleaseDomainProperties(ovk_domain *Domain, ovk_domain_properties **Properties);

void ovkGetNextAvailableGridID(const ovk_domain *Domain, int *GridID);

void ovkCreateGridParams(ovk_domain *Domain, ovk_grid_params **Params);
void ovkDestroyGridParams(ovk_domain *Domain, ovk_grid_params **Params);

void ovkCreateGridLocal(ovk_domain *Domain, int GridID, const ovk_grid_params *Params);
void ovkCreateGridRemote(ovk_domain *Domain, int GridID);
void ovkDestroyGrid(ovk_domain *Domain, int GridID);

bool ovkGridExists(const ovk_domain *Domain, int GridID);
void ovkGetGridInfo(const ovk_domain *Domain, int GridID, const ovk_grid_info **GridInfo);
bool ovkRankHasGrid(const ovk_domain *Domain, int GridID);
void ovkGetGrid(const ovk_domain *Domain, int GridID, const ovk_grid **Grid);
void ovkEditGridLocal(ovk_domain *Domain, int GridID, ovk_grid **Grid);
void ovkEditGridRemote(ovk_domain *Domain, int GridID);
void ovkReleaseGridLocal(ovk_domain *Domain, int GridID, ovk_grid **Grid);
void ovkReleaseGridRemote(ovk_domain *Domain, int GridID);

bool ovkConnectivityExists(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID);
bool ovkRankHasConnectivity(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID);
void ovkGetConnectivity(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  const ovk_connectivity **Connectivity);
void ovkEditConnectivityLocal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity);
void ovkEditConnectivityRemote(ovk_domain *Domain, int DonorGridID, int ReceiverGridID);
void ovkReleaseConnectivityLocal(ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_connectivity **Connectivity);
void ovkReleaseConnectivityRemote(ovk_domain *Domain, int DonorGridID, int ReceiverGridID);

bool ovkExchangeExists(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID);
bool ovkRankHasExchange(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID);
void ovkGetExchange(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  const ovk_exchange **Exchange);

void ovkGetLocalDonorCount(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  size_t *NumDonors);
void ovkGetLocalReceiverCount(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  size_t *NumReceivers);

void ovkAssemble(ovk_domain *Domain);

void ovkCollect(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_data_type DataType, int Count, ovk_collect_op CollectOp, const void **GridData,
  ovk_array_layout GridDataLayout, void **DonorData);

void ovkSend(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID, ovk_data_type DataType,
  int Count, const void **DonorData, int Tag, ovk_request **Request);

void ovkReceive(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_data_type DataType, int Count, void **ReceiverData, int Tag, ovk_request **Request);

void ovkWaitAll(int NumRequests, ovk_request **Requests);

void ovkWaitAny(int NumRequests, ovk_request **Requests, int *Index);

void ovkDisperse(const ovk_domain *Domain, int DonorGridID, int ReceiverGridID,
  ovk_data_type DataType, int Count, ovk_disperse_op DisperseOp, const void **ReceiverData,
  void **GridData, ovk_array_layout GridDataLayout);

void ovkGetDomainParamName(const ovk_domain_params *Params, char *Name);
void ovkSetDomainParamName(ovk_domain_params *Params, const char *Name);
void ovkGetDomainParamDimension(const ovk_domain_params *Params, int *NumDims);
void ovkSetDomainParamDimension(ovk_domain_params *Params, int NumDims);
void ovkGetDomainParamComm(const ovk_domain_params *Params, MPI_Comm *Comm);
void ovkSetDomainParamComm(ovk_domain_params *Params, MPI_Comm Comm);

void ovkGetDomainPropertyName(const ovk_domain_properties *Properties, char *Name);
void ovkGetDomainPropertyDimension(const ovk_domain_properties *Properties, int *NumDims);
void ovkGetDomainPropertyComm(const ovk_domain_properties *Properties, MPI_Comm *Comm);
void ovkGetDomainPropertyCommSize(const ovk_domain_properties *Properties, int *CommSize);
void ovkGetDomainPropertyCommRank(const ovk_domain_properties *Properties, int *CommRank);
void ovkGetDomainPropertyGridCount(const ovk_domain_properties *Properties, int *NumGrids);

#ifdef __cplusplus
}
#endif

#endif
