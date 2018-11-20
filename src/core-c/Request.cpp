// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Request.h"

#include "ovk/core-c/Global.h"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Request.hpp"

extern "C" {

void ovkWait(ovk_request **Request) {

  OVK_DEBUG_ASSERT(Request, "Invalid request pointer.");

  if (*Request) {
    auto *RequestCPPPtr = reinterpret_cast<ovk::request *>(*Request);
    RequestCPPPtr->Wait();
    delete RequestCPPPtr;
    RequestCPPPtr = nullptr;
  }

}

void ovkWaitAll(int NumRequests, ovk_request **Requests) {

  OVK_DEBUG_ASSERT(NumRequests >= 0, "Invalid request count.");
  OVK_DEBUG_ASSERT(Requests || NumRequests == 0, "Invalid requests pointer.");
  // Note: Not checking Requests[i] here on purpose -- allowed to be null

  auto **RequestsCPPPtr = reinterpret_cast<ovk::request **>(Requests);
  ovk::request::core_WaitAll(NumRequests, RequestsCPPPtr);

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    delete RequestsCPPPtr[iRequest];
    RequestsCPPPtr[iRequest] = nullptr;
  }

}

void ovkWaitAny(int NumRequests, ovk_request **Requests, int *Index) {

  OVK_DEBUG_ASSERT(NumRequests >= 0, "Invalid request count.");
  OVK_DEBUG_ASSERT(Requests || NumRequests == 0, "Invalid requests pointer.");
  // Note: Not checking Requests[i] here on purpose -- allowed to be null
  OVK_DEBUG_ASSERT(Index, "Invalid index pointer.");

  auto **RequestsCPPPtr = reinterpret_cast<ovk::request **>(Requests);
  ovk::request::core_WaitAny(NumRequests, RequestsCPPPtr, *Index);

  if (*Index >= 0) {
    delete RequestsCPPPtr[*Index];
    RequestsCPPPtr[*Index] = nullptr;
  }

}

}
