// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PARALLEL_UTILS_INCLUDED
#define OVK_CORE_PARALLEL_UTILS_INCLUDED

#include "ovk/core/Global.h"
#include "ovk/core/OrderedMap.h"

#ifdef __cplusplus
extern "C" {
#endif

void PRIVATE(BroadcastAnySource)(void *Data, int Count, MPI_Datatype DataType, bool IsSource,
  MPI_Comm Comm);
#define BroadcastAnySource(...) PRIVATE(BroadcastAnySource)(__VA_ARGS__)

struct t_signal;
typedef struct t_signal t_signal;

void PRIVATE(CreateSignal)(t_signal **Signal, MPI_Comm Comm);
#define CreateSignal(...) PRIVATE(CreateSignal)(__VA_ARGS__)

void PRIVATE(StartSignal)(t_signal *Signal);
#define StartSignal(...) PRIVATE(StartSignal)(__VA_ARGS__)

void PRIVATE(CheckSignal)(t_signal *Signal, bool *Done);
#define CheckSignal(...) PRIVATE(CheckSignal)(__VA_ARGS__)

void PRIVATE(DestroySignal)(t_signal **Signal);
#define DestroySignal(...) PRIVATE(DestroySignal)(__VA_ARGS__)

void PRIVATE(DynamicHandshake)(MPI_Comm Comm, int NumDestRanks, const int *DestRanks,
  t_ordered_map *SourceRanks);
#define DynamicHandshake(...) PRIVATE(DynamicHandshake)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
