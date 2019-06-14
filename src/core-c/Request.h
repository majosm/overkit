// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_C_REQUEST_H_INCLUDED
#define OVK_CORE_C_REQUEST_H_INCLUDED

#include <ovk/core-c/Global.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_request;
typedef struct ovk_request ovk_request;

void ovkWait(ovk_request **Request);
void ovkWaitAll(int NumRequests, ovk_request **Requests);
void ovkWaitAny(int NumRequests, ovk_request **Requests, int *Index);

#ifdef __cplusplus
}
#endif

#endif
