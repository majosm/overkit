// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/Request.h"

#include "ovk/core-c/Global.h"
#include "ovk/core/Array.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Request.hpp"

extern "C" {

void ovkRequestWait(ovk_request **Request) {

  OVK_DEBUG_ASSERT(Request, "Invalid request pointer.");

  if (*Request) {
    auto *RequestCPPPtr = reinterpret_cast<ovk::request *>(*Request);
    RequestCPPPtr->Wait();
    delete RequestCPPPtr;
    *Request = nullptr;
  }

}

void ovkRequestWaitAll(int NumRequests, ovk_request **Requests) {

  OVK_DEBUG_ASSERT(NumRequests >= 0, "Invalid request count.");
  OVK_DEBUG_ASSERT(Requests || NumRequests == 0, "Invalid requests pointer.");
  // Note: Not checking Requests[i] here on purpose -- allowed to be null

  ovk::array<ovk::request *> RequestCPPPtrs({NumRequests});
  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    RequestCPPPtrs(iRequest) = reinterpret_cast<ovk::request *>(Requests[iRequest]);
  }
  ovk::RequestWaitAll(RequestCPPPtrs);

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    delete RequestCPPPtrs(iRequest);
    Requests[iRequest] = nullptr;
  }

}

void ovkRequestWaitAny(int NumRequests, ovk_request **Requests, int *Index) {

  OVK_DEBUG_ASSERT(NumRequests >= 0, "Invalid request count.");
  OVK_DEBUG_ASSERT(Requests || NumRequests == 0, "Invalid requests pointer.");
  // Note: Not checking Requests[i] here on purpose -- allowed to be null
  OVK_DEBUG_ASSERT(Index, "Invalid index pointer.");

  ovk::array<ovk::request *> RequestCPPPtrs({NumRequests});
  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    RequestCPPPtrs(iRequest) = reinterpret_cast<ovk::request *>(Requests[iRequest]);
  }
  ovk::RequestWaitAny(RequestCPPPtrs, *Index);

  if (*Index >= 0) {
    delete RequestCPPPtrs(*Index);
    Requests[*Index] = nullptr;
  }

}

}
