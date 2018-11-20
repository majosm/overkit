// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_EXCHANGE_H_INCLUDED
#define OVK_CORE_C_EXCHANGE_H_INCLUDED

#include <ovk/core-c/Global.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_exchange;
typedef struct ovk_exchange ovk_exchange;

struct ovk_exchange_info;
typedef struct ovk_exchange_info ovk_exchange_info;

bool ovkRankHasExchangeDonorSide(const ovk_exchange *Exchange);
bool ovkRankHasExchangeReceiverSide(const ovk_exchange *Exchange);

void ovkGetExchangeInfoDonorGridID(const ovk_exchange_info *Info, int *DonorGridID);
void ovkGetExchangeInfoReceiverGridID(const ovk_exchange_info *Info, int *ReceiverGridID);
void ovkGetExchangeInfoName(const ovk_exchange_info *Info, char *Name);
void ovkGetExchangeInfoDimension(const ovk_exchange_info *Info, int *NumDims);
void ovkGetExchangeInfoRootRank(const ovk_exchange_info *Info, int *RootRank);
void ovkGetExchangeInfoIsLocal(const ovk_exchange_info *Info, bool *IsLocal);

#ifdef __cplusplus
}
#endif

#endif
