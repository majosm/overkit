// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Exchanger.h"

#include "ovk/core-c/Domain.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Domain.hpp"
#include "ovk/core/Exchanger.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/ID.hpp"

#include <mpi.h>

#include <cstring>
#include <string>
#include <utility>

extern "C" {

void ovkCreateExchanger(ovk_exchanger **Exchanger, ovk_shared_context *Context, ovk_exchanger_params
  **Params) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ContextCPP = *reinterpret_cast<std::shared_ptr<ovk::context> *>(Context);

  ovk::exchanger::params *ParamsCPPPtr = nullptr;
  if (Params && *Params) {
    ParamsCPPPtr = reinterpret_cast<ovk::exchanger::params *>(*Params);
  }

  ovk::exchanger *ExchangerCPPPtr;
  if (ParamsCPPPtr) {
    ExchangerCPPPtr = new ovk::exchanger(ovk::CreateExchanger(ContextCPP, std::move(
      *ParamsCPPPtr)));
  } else {
    ExchangerCPPPtr = new ovk::exchanger(ovk::CreateExchanger(ContextCPP));
  }

  *Exchanger = reinterpret_cast<ovk_exchanger *>(ExchangerCPPPtr);

  if (ParamsCPPPtr) {
    delete ParamsCPPPtr;
    *Params = nullptr;
  }

}

void ovkDestroyExchanger(ovk_exchanger **Exchanger) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(*Exchanger, "Invalid exchanger pointer.");

  auto ExchangerCPPPtr = reinterpret_cast<ovk::exchanger *>(*Exchanger);

  delete ExchangerCPPPtr;

  *Exchanger = nullptr;

}

void ovkGetExchangerContextC(const ovk_exchanger *Exchanger, const ovk_context **Context) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  *Context = reinterpret_cast<const ovk_context *>(&ExchangerCPP.Context());

}

void ovkGetExchangerContext(ovk_exchanger *Exchanger, ovk_context **Context) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  *Context = reinterpret_cast<ovk_context *>(&ExchangerCPP.Context());

}

void ovkGetExchangerSharedContext(ovk_exchanger *Exchanger, ovk_shared_context **Context) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(Context, "Invalid context pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  auto &ContextCPP = ExchangerCPP.SharedContext();

  auto ContextCPPPtr = new std::shared_ptr<ovk::context>(ContextCPP);

  *Context = reinterpret_cast<ovk_shared_context *>(ContextCPPPtr);

}

bool ovkExchangerIsBound(const ovk_exchanger *Exchanger) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  return ExchangerCPP.Bound();

}

void ovkBindExchanger(ovk_exchanger *Exchanger, const ovk_domain *Domain, ovk_exchanger_bindings
  **Bindings) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");
  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");
  OVK_DEBUG_ASSERT(*Bindings, "Invalid bindings pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  auto BindingsCPPPtr = reinterpret_cast<ovk::exchanger::bindings *>(*Bindings);

  ExchangerCPP.Bind(DomainCPP, std::move(*BindingsCPPPtr));

  delete BindingsCPPPtr;

  *Bindings = nullptr;

}

void ovkUnbindExchanger(ovk_exchanger *Exchanger) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  ExchangerCPP.Unbind();

}

void ovkGetExchangerDomain(const ovk_exchanger *Exchanger, const ovk_domain **Domain) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(Domain, "Invalid domain pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  *Domain = reinterpret_cast<const ovk_domain *>(&ExchangerCPP.Domain());

}

bool ovkExchangerCollectExists(const ovk_exchanger *Exchanger, int MGridID, int NGridID, int
  CollectID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  return ExchangerCPP.CollectExists(MGridID, NGridID, CollectID);

}

void ovkGetNextAvailableExchangerCollectID(const ovk_exchanger *Exchanger, int MGridID, int NGridID,
  int *CollectID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  *CollectID = ovk::NextAvailableID(ExchangerCPP.CollectIDs(MGridID, NGridID));

}

void ovkCreateExchangerCollect(ovk_exchanger *Exchanger, int MGridID, int NGridID, int CollectID,
  ovk_collect_op CollectOp, ovk_data_type ValueType, int Count, const int *GridValuesBegin, const
  int *GridValuesEnd, ovk_array_layout GridValuesLayout) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(GridValuesBegin, "Invalid grid values begin pointer.");
  OVK_DEBUG_ASSERT(GridValuesEnd, "Invalid grid values end pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);

  int NumDims = ExchangerCPP.Domain().Dimension();

  ovk::range GridValuesRange = ovk::MakeEmptyRange(NumDims);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    GridValuesRange.Begin(iDim) = GridValuesBegin[iDim];
    GridValuesRange.End(iDim) = GridValuesEnd[iDim];
  }

  ExchangerCPP.CreateCollect(MGridID, NGridID, CollectID, ovk::collect_op(CollectOp),
    ovk::data_type(ValueType), Count, GridValuesRange, ovk::array_layout(GridValuesLayout));

}

void ovkDestroyExchangerCollect(ovk_exchanger *Exchanger, int MGridID, int NGridID, int CollectID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  ExchangerCPP.DestroyCollect(MGridID, NGridID, CollectID);

}

void ovkExchangerCollect(ovk_exchanger *Exchanger, int MGridID, int NGridID, int CollectID, const
  void *GridValues, void *DonorValues) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  ExchangerCPP.Collect(MGridID, NGridID, CollectID, GridValues, DonorValues);

}

bool ovkExchangerSendExists(const ovk_exchanger *Exchanger, int MGridID, int NGridID, int SendID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  return ExchangerCPP.SendExists(MGridID, NGridID, SendID);

}

void ovkGetNextAvailableExchangerSendID(const ovk_exchanger *Exchanger, int MGridID, int NGridID,
  int *SendID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  *SendID = ovk::NextAvailableID(ExchangerCPP.SendIDs(MGridID, NGridID));

}

void ovkCreateExchangerSend(ovk_exchanger *Exchanger, int MGridID, int NGridID, int SendID,
  ovk_data_type ValueType, int Count, int Tag) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  ExchangerCPP.CreateSend(MGridID, NGridID, SendID, ovk::data_type(ValueType), Count, Tag);

}

void ovkDestroyExchangerSend(ovk_exchanger *Exchanger, int MGridID, int NGridID, int SendID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  ExchangerCPP.DestroySend(MGridID, NGridID, SendID);

}

void ovkExchangerSend(ovk_exchanger *Exchanger, int MGridID, int NGridID, int SendID, const void
  *DonorValues, ovk_request **Request) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(Request, "Invalid request pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  auto RequestCPPPtr = new ovk::request();

  *RequestCPPPtr = ExchangerCPP.Send(MGridID, NGridID, SendID, DonorValues);

  *Request = reinterpret_cast<ovk_request *>(RequestCPPPtr);

}

bool ovkExchangerReceiveExists(const ovk_exchanger *Exchanger, int MGridID, int NGridID, int RecvID)
  {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  return ExchangerCPP.ReceiveExists(MGridID, NGridID, RecvID);

}

void ovkGetNextAvailableExchangerReceiveID(const ovk_exchanger *Exchanger, int MGridID, int NGridID,
  int *RecvID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  *RecvID = ovk::NextAvailableID(ExchangerCPP.ReceiveIDs(MGridID, NGridID));

}

void ovkCreateExchangerReceive(ovk_exchanger *Exchanger, int MGridID, int NGridID, int RecvID,
  ovk_data_type ValueType, int Count, int Tag) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  ExchangerCPP.CreateReceive(MGridID, NGridID, RecvID, ovk::data_type(ValueType), Count, Tag);

}

void ovkDestroyExchangerReceive(ovk_exchanger *Exchanger, int MGridID, int NGridID, int RecvID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  ExchangerCPP.DestroyReceive(MGridID, NGridID, RecvID);

}

void ovkExchangerReceive(ovk_exchanger *Exchanger, int MGridID, int NGridID, int RecvID, void
  *ReceiverValues, ovk_request **Request) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(Request, "Invalid request pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  auto RequestCPPPtr = new ovk::request();

  *RequestCPPPtr = ExchangerCPP.Receive(MGridID, NGridID, RecvID, ReceiverValues);

  *Request = reinterpret_cast<ovk_request *>(RequestCPPPtr);

}

bool ovkExchangerDisperseExists(const ovk_exchanger *Exchanger, int MGridID, int NGridID, int
  DisperseID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  return ExchangerCPP.DisperseExists(MGridID, NGridID, DisperseID);

}

void ovkGetNextAvailableExchangerDisperseID(const ovk_exchanger *Exchanger, int MGridID, int
  NGridID, int *DisperseID) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<const ovk::exchanger *>(Exchanger);
  *DisperseID = ovk::NextAvailableID(ExchangerCPP.DisperseIDs(MGridID, NGridID));

}

void ovkCreateExchangerDisperse(ovk_exchanger *Exchanger, int MGridID, int NGridID, int DisperseID,
  ovk_disperse_op DisperseOp, ovk_data_type ValueType, int Count, const int *GridValuesBegin, const
  int *GridValuesEnd, ovk_array_layout GridValuesLayout) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");
  OVK_DEBUG_ASSERT(GridValuesBegin, "Invalid grid values begin pointer.");
  OVK_DEBUG_ASSERT(GridValuesEnd, "Invalid grid values end pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);

  int NumDims = ExchangerCPP.Domain().Dimension();

  ovk::range GridValuesRange = ovk::MakeEmptyRange(NumDims);
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    GridValuesRange.Begin(iDim) = GridValuesBegin[iDim];
    GridValuesRange.End(iDim) = GridValuesEnd[iDim];
  }

  ExchangerCPP.CreateDisperse(MGridID, NGridID, DisperseID, ovk::disperse_op(DisperseOp),
    ovk::data_type(ValueType), Count, GridValuesRange, ovk::array_layout(GridValuesLayout));

}

void ovkDestroyExchangerDisperse(ovk_exchanger *Exchanger, int MGridID, int NGridID, int DisperseID)
  {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  ExchangerCPP.DestroyDisperse(MGridID, NGridID, DisperseID);

}

void ovkExchangerDisperse(ovk_exchanger *Exchanger, int MGridID, int NGridID, int DisperseID, const
  void *ReceiverValues, void *GridValues) {

  OVK_DEBUG_ASSERT(Exchanger, "Invalid exchanger pointer.");

  auto &ExchangerCPP = *reinterpret_cast<ovk::exchanger *>(Exchanger);
  ExchangerCPP.Disperse(MGridID, NGridID, DisperseID, ReceiverValues, GridValues);

}

void ovkCreateExchangerParams(ovk_exchanger_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");

  auto ParamsCPPPtr = new ovk::exchanger::params();

  *Params = reinterpret_cast<ovk_exchanger_params *>(ParamsCPPPtr);

}

void ovkDestroyExchangerParams(ovk_exchanger_params **Params) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(*Params, "Invalid params pointer.");

  auto ParamsCPPPtr = reinterpret_cast<ovk::exchanger::params *>(*Params);

  delete ParamsCPPPtr;

  *Params = nullptr;

}

void ovkGetExchangerParamName(const ovk_exchanger_params *Params, char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<const ovk::exchanger::params *>(Params);
  std::strcpy(Name, ParamsCPP.Name().c_str());

}

void ovkSetExchangerParamName(ovk_exchanger_params *Params, const char *Name) {

  OVK_DEBUG_ASSERT(Params, "Invalid params pointer.");
  OVK_DEBUG_ASSERT(Name, "Invalid name pointer.");

  auto &ParamsCPP = *reinterpret_cast<ovk::exchanger::params *>(Params);
  ParamsCPP.SetName(Name);

}

void ovkCreateExchangerBindings(ovk_exchanger_bindings **Bindings) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");

  auto BindingsCPPPtr = new ovk::exchanger::bindings();

  *Bindings = reinterpret_cast<ovk_exchanger_bindings *>(BindingsCPPPtr);

}

void ovkDestroyExchangerBindings(ovk_exchanger_bindings **Bindings) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");
  OVK_DEBUG_ASSERT(*Bindings, "Invalid bindings pointer.");

  auto BindingsCPPPtr = reinterpret_cast<ovk::exchanger::bindings *>(*Bindings);

  delete BindingsCPPPtr;

  *Bindings = nullptr;

}

void ovkGetExchangerBindingsConnectivityComponentID(ovk_exchanger_bindings *Bindings, int
  *ConnectivityComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");
  OVK_DEBUG_ASSERT(ConnectivityComponentID, "Invalid connectivity component ID pointer.");

  auto &BindingsCPP = *reinterpret_cast<const ovk::exchanger::bindings *>(Bindings);
  *ConnectivityComponentID = BindingsCPP.ConnectivityComponentID();

}

void ovkSetExchangerBindingsConnectivityComponentID(ovk_exchanger_bindings *Bindings, int
  ConnectivityComponentID) {

  OVK_DEBUG_ASSERT(Bindings, "Invalid bindings pointer.");

  auto &BindingsCPP = *reinterpret_cast<ovk::exchanger::bindings *>(Bindings);
  BindingsCPP.SetConnectivityComponentID(ConnectivityComponentID);

}

}
