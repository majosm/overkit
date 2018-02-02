// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_EXCHANGE_INCLUDED
#define OVK_CORE_EXCHANGE_INCLUDED

#include "ovkExchange.h"

#include "Connectivity.h"
#include "ErrorHandler.h"
#include "Global.h"
#include "Logger.h"
#include "PartitionHash.h"

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_exchange {
  const ovk_connectivity *connectivity;
  int num_dims;
  MPI_Comm comm;
  int comm_size;
  int comm_rank;
  t_logger *logger;
  t_error_handler *error_handler;
  int num_collect_sends;
  int *collect_send_dest_ranks;
  size_t *num_collect_send_points;
  int ***collect_send_points;
  int num_collect_recvs;
  int *collect_recv_source_ranks;
  size_t *num_collect_recv_points;
  int ***collect_recv_points;
  int *num_remote_donor_points;
  size_t **remote_donor_points;
  int **remote_donor_point_collect_recv_indices;
  size_t **remote_donor_point_collect_recv_buffer_offsets;
  size_t *donors_sorted;
  size_t *receivers_sorted;
  int *donor_dest_ranks;
  int *receiver_source_ranks;
  int num_sends;
  int *send_ranks;
  size_t *send_counts;
  int *donor_send_indices;
  int num_recvs;
  int *recv_ranks;
  size_t *recv_counts;
  int *receiver_recv_indices;
  t_partition_hash *source_hash;
  t_partition_hash *destination_hash;
};

struct ovk_exchange_info {
  int donor_grid_id;
  int receiver_grid_id;
  char name[OVK_NAME_LENGTH];
  int num_dims;
  int root_rank;
};

typedef enum {
  SEND_REQUEST,
  RECV_REQUEST
} t_request_type;

struct ovk_request {
  t_request_type type;
  void *data;
};

void PRIVATE(CreateExchange)(ovk_exchange **Exchange, const ovk_connectivity *Connectivity,
  t_logger *Logger, t_error_handler *ErrorHandler);
#define CreateExchange(...) PRIVATE(CreateExchange)(__VA_ARGS__)
void PRIVATE(DestroyExchange)(ovk_exchange **Exchange);
#define DestroyExchange(...) PRIVATE(DestroyExchange)(__VA_ARGS__)

void PRIVATE(CreateExchangeInfo)(ovk_exchange_info **Info, const ovk_exchange *Connectivity,
  MPI_Comm Comm, int CommRank);
#define CreateExchangeInfo(...) PRIVATE(CreateExchangeInfo)(__VA_ARGS__)
void PRIVATE(DestroyExchangeInfo)(ovk_exchange_info **Info);
#define DestroyExchangeInfo(...) PRIVATE(DestroyExchangeInfo)(__VA_ARGS__)

void PRIVATE(UpdateExchange)(ovk_exchange *Exchange);
#define UpdateExchange(...) PRIVATE(UpdateExchange)(__VA_ARGS__)

void PRIVATE(ExchangeCollect)(const ovk_exchange *Exchange, ovk_data_type DataType, int Count,
  ovk_collect_op CollectOp, const void **GridData, ovk_array_layout GridDataLayout,
  void **DonorData);
#define ExchangeCollect(...) PRIVATE(ExchangeCollect)(__VA_ARGS__)
void PRIVATE(ExchangeSend)(const ovk_exchange *Exchange, ovk_data_type DataType, int Count,
  const void **DonorData, int Tag, ovk_request **Request);
#define ExchangeSend(...) PRIVATE(ExchangeSend)(__VA_ARGS__)
void PRIVATE(ExchangeReceive)(const ovk_exchange *Exchange, ovk_data_type DataType, int Count,
  void **ReceiverData, int Tag, ovk_request **Request);
#define ExchangeReceive(...) PRIVATE(ExchangeReceive)(__VA_ARGS__)
void PRIVATE(ExchangeWaitAll)(int NumRequests, ovk_request **Requests);
#define ExchangeWaitAll(...) PRIVATE(ExchangeWaitAll)(__VA_ARGS__)
void PRIVATE(ExchangeWaitAny)(int NumRequests, ovk_request **Requests, int *Index);
#define ExchangeWaitAny(...) PRIVATE(ExchangeWaitAny)(__VA_ARGS__)
void PRIVATE(ExchangeDisperse)(const ovk_exchange *Exchange, ovk_data_type DataType, int Count,
  ovk_disperse_op DisperseOp, const void **ReceiverData, void **GridData,
  ovk_array_layout GridDataLayout);
#define ExchangeDisperse(...) PRIVATE(ExchangeDisperse)(__VA_ARGS__)

#ifdef __cplusplus
}
#endif

#endif
