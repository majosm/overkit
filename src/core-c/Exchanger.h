// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_EXCHANGER_H_INCLUDED
#define OVK_CORE_C_EXCHANGER_H_INCLUDED

#include <ovk/core-c/Context.h>
#include <ovk/core-c/DataType.h>
#include <ovk/core-c/Domain.h>
#include <ovk/core-c/Global.h>
#include <ovk/core-c/Request.h>
#include <ovk/core/Exchanger.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_exchanger;
typedef struct ovk_exchanger ovk_exchanger;

struct ovk_exchanger_params;
typedef struct ovk_exchanger_params ovk_exchanger_params;

struct ovk_exchanger_bindings;
typedef struct ovk_exchanger_bindings ovk_exchanger_bindings;

void ovkCreateExchanger(ovk_exchanger **Exchanger, ovk_shared_context *Context, ovk_exchanger_params
  **Params);
void ovkDestroyExchanger(ovk_exchanger **Exchanger);

void ovkGetExchangerContextC(const ovk_exchanger *Exchanger, const ovk_context **Context);
void ovkGetExchangerContext(ovk_exchanger *Exchanger, ovk_context **Context);
void ovkGetExchangerSharedContext(ovk_exchanger *Exchanger, ovk_shared_context **Context);

bool ovkExchangerIsBound(const ovk_exchanger *Exchanger);
void ovkBindExchanger(ovk_exchanger *Exchanger, const ovk_domain *Domain, ovk_exchanger_bindings
  **Bindings);
void ovkUnbindExchanger(ovk_exchanger *Exchanger);

void ovkGetExchangerDomain(const ovk_exchanger *Exchanger, const ovk_domain **Domain);

bool ovkExchangerCollectExists(const ovk_exchanger *Exchanger, int MGridID, int NGridID, int
  CollectID);
void ovkGetNextAvailableExchangerCollectID(const ovk_exchanger *Exchanger, int MGridID, int NGridID,
  int *CollectID);
void ovkCreateExchangerCollect(ovk_exchanger *Exchanger, int MGridID, int NGridID, int CollectID,
  ovk_collect_op CollectOp, ovk_data_type ValueType, int Count, const int *GridValuesBegin, const
  int *GridValuesEnd, ovk_array_layout GridValuesLayout);
void ovkDestroyExchangerCollect(ovk_exchanger *Exchanger, int MGridID, int NGridID, int CollectID);
void ovkExchangerCollect(ovk_exchanger *Exchanger, int MGridID, int NGridID, int CollectID, const
  void **GridValues, void **DonorValues);

bool ovkExchangerSendExists(const ovk_exchanger *Exchanger, int MGridID, int NGridID, int SendID);
void ovkGetNextAvailableExchangerSendID(const ovk_exchanger *Exchanger, int MGridID, int NGridID,
  int *SendID);
void ovkCreateExchangerSend(ovk_exchanger *Exchanger, int MGridID, int NGridID, int SendID,
  ovk_data_type ValueType, int Count, int Tag);
void ovkDestroyExchangerSend(ovk_exchanger *Exchanger, int MGridID, int NGridID, int SendID);
void ovkExchangerSend(ovk_exchanger *Exchanger, int MGridID, int NGridID, int SendID, const void
  **DonorValues, ovk_request **Request);

bool ovkExchangerReceiveExists(const ovk_exchanger *Exchanger, int MGridID, int NGridID, int
  RecvID);
void ovkGetNextAvailableExchangerReceiveID(const ovk_exchanger *Exchanger, int MGridID, int NGridID,
  int *RecvID);
void ovkCreateExchangerReceive(ovk_exchanger *Exchanger, int MGridID, int NGridID, int RecvID,
  ovk_data_type ValueType, int Count, int Tag);
void ovkDestroyExchangerReceive(ovk_exchanger *Exchanger, int MGridID, int NGridID, int RecvID);
void ovkExchangerReceive(ovk_exchanger *Exchanger, int MGridID, int NGridID, int RecvID, void
  **ReceiverValues, ovk_request **Request);

bool ovkExchangerDisperseExists(const ovk_exchanger *Exchanger, int MGridID, int NGridID, int
  DisperseID);
void ovkGetNextAvailableExchangerDisperseID(const ovk_exchanger *Exchanger, int MGridID, int
  NGridID, int *DisperseID);
void ovkCreateExchangerDisperse(ovk_exchanger *Exchanger, int MGridID, int NGridID, int DisperseID,
  ovk_disperse_op DisperseOp, ovk_data_type ValueType, int Count, const int *GridValuesBegin, const
  int *GridValuesEnd, ovk_array_layout GridValuesLayout);
void ovkDestroyExchangerDisperse(ovk_exchanger *Exchanger, int MGridID, int NGridID, int
  DisperseID);
void ovkExchangerDisperse(ovk_exchanger *Exchanger, int MGridID, int NGridID, int DisperseID,
  const void **ReceiverValues, void **GridValues);

void ovkCreateExchangerParams(ovk_exchanger_params **Params);
void ovkDestroyExchangerParams(ovk_exchanger_params **Params);
void ovkGetExchangerParamName(const ovk_exchanger_params *Params, char *Name);
void ovkSetExchangerParamName(ovk_exchanger_params *Params, const char *Name);

void ovkCreateExchangerBindings(ovk_exchanger_bindings **Bindings);
void ovkDestroyExchangerBindings(ovk_exchanger_bindings **Bindings);
void ovkGetExchangerBindingsConnectivityComponentID(ovk_exchanger_bindings *Bindings, int
  *ConnectivityComponentID);
void ovkSetExchangerBindingsConnectivityComponentID(ovk_exchanger_bindings *Bindings, int
  ConnectivityComponentID);

#ifdef __cplusplus
}
#endif

#endif
