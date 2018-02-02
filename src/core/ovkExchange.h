// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_EXCHANGE_INCLUDED
#define OVK_CORE_PUBLIC_EXCHANGE_INCLUDED

#include <ovkGlobal.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_exchange;
typedef struct ovk_exchange ovk_exchange;

struct ovk_exchange_info;
typedef struct ovk_exchange_info ovk_exchange_info;

struct ovk_request;
typedef struct ovk_request ovk_request;

bool ovkRankHasExchangeDonorSide(const ovk_exchange *Exchange);
bool ovkRankHasExchangeReceiverSide(const ovk_exchange *Exchange);

#ifdef __cplusplus
}
#endif

#endif
