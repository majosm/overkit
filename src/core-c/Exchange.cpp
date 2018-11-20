// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Exchange.h"

#include "ovk/core-c/Global.h"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"

#include <cstring>
#include <string>

bool ovkRankHasExchangeDonorSide(const ovk_exchange *Exchange) {

  OVK_DEBUG_ASSERT(Exchange, "Invalid exchange pointer.");

  auto &ExchangeCPP = *reinterpret_cast<const ovk::exchange *>(Exchange);
  return ovk::RankHasExchangeDonorSide(ExchangeCPP);

}

bool ovkRankHasExchangeReceiverSide(const ovk_exchange *Exchange) {

  OVK_DEBUG_ASSERT(Exchange, "Invalid exchange pointer.");

  auto &ExchangeCPP = *reinterpret_cast<const ovk::exchange *>(Exchange);
  return ovk::RankHasExchangeReceiverSide(ExchangeCPP);

}

void ovkGetExchangeInfoName(const ovk_exchange_info *Info, char *Name) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::exchange_info *>(Info);

  std::string NameCPP;
  ovk::GetExchangeInfoName(InfoCPP, NameCPP);

  std::strcpy(Name, NameCPP.c_str());

}

void ovkGetExchangeInfoDimension(const ovk_exchange_info *Info, int *NumDims) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::exchange_info *>(Info);
  ovk::GetExchangeInfoDimension(InfoCPP, *NumDims);

}

void ovkGetExchangeInfoRootRank(const ovk_exchange_info *Info, int *RootRank) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(RootRank, "Invalid root rank pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::exchange_info *>(Info);
  ovk::GetExchangeInfoRootRank(InfoCPP, *RootRank);

}

void ovkGetExchangeInfoIsLocal(const ovk_exchange_info *Info, bool *IsLocal) {

  OVK_DEBUG_ASSERT(Info, "Invalid info pointer.");
  OVK_DEBUG_ASSERT(IsLocal, "Invalid is local pointer.");

  auto &InfoCPP = *reinterpret_cast<const ovk::exchange_info *>(Info);
  ovk::GetExchangeInfoIsLocal(InfoCPP, *IsLocal);

}
