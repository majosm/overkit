// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Request.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <memory>
#include <tuple>
#include <utility>

namespace ovk {

void RequestWaitAll(array_view<request> Requests) {

  request::internal_WaitAll(Requests);

}

void RequestWaitAny(array_view<request> Requests, int &Index) {

  request::internal_WaitAny(Requests, Index);

}

void RequestWaitAll(array_view<request *> Requests) {

  request::internal_WaitAll(Requests);

}

void RequestWaitAny(array_view<request *> Requests, int &Index) {

  request::internal_WaitAny(Requests, Index);

}

void request::internal_WaitAll(array_view<request> Requests) {

  OVK_DEBUG_ASSERT(Requests || Requests.Count() == 0, "Invalid requests array.");

  int NumRequests = Requests.Count();

  array<request *> RequestPtrs({NumRequests});

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    RequestPtrs[iRequest] = Requests.Data(iRequest);
  }

  request::internal_WaitAll(RequestPtrs);

}

void request::internal_WaitAny(array_view<request> Requests, int &Index) {

  OVK_DEBUG_ASSERT(Requests || Requests.Count() == 0, "Invalid requests array.");

  int NumRequests = Requests.Count();

  array<request *> RequestPtrs({NumRequests});

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    RequestPtrs[iRequest] = Requests.Data(iRequest);
  }

  request::internal_WaitAny(RequestPtrs, Index);

}

void request::internal_WaitAll(array_view<request *> Requests) {

  OVK_DEBUG_ASSERT(Requests || Requests.Count() == 0, "Invalid requests array.");

  for (auto Request : Requests) {
    if (Request && *Request) {
      Request->StartWaitTime_();
    }
  }

  auto StartMPITime = [&]() {
    for (auto Request : Requests) {
      if (Request && *Request) {
        Request->StartMPITime_();
      }
    }
  };

  auto StopMPITime = [&]() {
    for (auto Request : Requests) {
      if (Request && *Request) {
        Request->StopMPITime_();
      }
    }
  };

  int NumRequests = Requests.Count();

  array<int> NumRemainingMPIRequests({NumRequests});

  int TotalMPIRequests = 0;

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests(iRequest)) {
      request &Request = *Requests(iRequest);
      if (Request) {
        int NumMPIRequests = Request.MPIRequests_().Count();
        NumRemainingMPIRequests(iRequest) = NumMPIRequests;
        TotalMPIRequests += NumMPIRequests;
      }
    }
  }

  array<MPI_Request> AllMPIRequests;
  array<std::tuple<int,int>> AllMPIRequestToRequest;

  AllMPIRequests.Reserve(TotalMPIRequests);
  AllMPIRequestToRequest.Reserve(TotalMPIRequests);

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests(iRequest)) {
      request &Request = *Requests(iRequest);
      if (Request) {
        array_view<MPI_Request> MPIRequests = Request.MPIRequests_();
        for (int iMPIRequest = 0; iMPIRequest < MPIRequests.Count(); ++iMPIRequest) {
          AllMPIRequests.Append(MPIRequests(iMPIRequest));
          AllMPIRequestToRequest.Append(iRequest, iMPIRequest);
        }
      }
    }
  }

  while (true) {
    int iCompleted;
    StartMPITime();
    MPI_Waitany(TotalMPIRequests, AllMPIRequests.Data(), &iCompleted, MPI_STATUSES_IGNORE);
    StopMPITime();
    if (iCompleted == MPI_UNDEFINED) {
      break;
    }
    int iRequest, iMPIRequest;
    std::tie(iRequest, iMPIRequest) = AllMPIRequestToRequest(iCompleted);
    Requests(iRequest)->MPIRequests_()(iMPIRequest) = MPI_REQUEST_NULL;
    Requests(iRequest)->Finish_(iMPIRequest);
    --NumRemainingMPIRequests(iRequest);
    if (NumRemainingMPIRequests(iRequest) == 0) {
      Requests(iRequest)->Wait();
    }
  }

  for (auto Request : Requests) {
    if (Request && *Request) {
      Request->StopWaitTime_();
    }
  }

}

void request::internal_WaitAny(array_view<request *> Requests, int &Index) {

  OVK_DEBUG_ASSERT(Requests || Requests.Count() == 0, "Invalid requests array.");

  for (auto Request : Requests) {
    if (Request && *Request) {
      Request->StartWaitTime_();
    }
  }

  auto StartMPITime = [&]() {
    for (auto Request : Requests) {
      if (Request && *Request) {
        Request->StartMPITime_();
      }
    }
  };

  auto StopMPITime = [&]() {
    for (auto Request : Requests) {
      if (Request && *Request) {
        Request->StopMPITime_();
      }
    }
  };

  int NumRequests = Requests.Count();

  array<int> NumRemainingMPIRequests({NumRequests});

  int TotalMPIRequests = 0;

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests(iRequest)) {
      request &Request = *Requests(iRequest);
      if (Request) {
        int NumMPIRequests = Request.MPIRequests_().Count();
        NumRemainingMPIRequests(iRequest) = NumMPIRequests;
        TotalMPIRequests += NumMPIRequests;
      }
    }
  }

  array<MPI_Request> AllMPIRequests;
  array<std::tuple<int,int>> AllMPIRequestToRequest;

  AllMPIRequests.Reserve(TotalMPIRequests);
  AllMPIRequestToRequest.Reserve(TotalMPIRequests);

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests(iRequest)) {
      request &Request = *Requests(iRequest);
      if (Request) {
        array_view<MPI_Request> MPIRequests = Request.MPIRequests_();
        for (int iMPIRequest = 0; iMPIRequest < MPIRequests.Count(); ++iMPIRequest) {
          AllMPIRequests.Append(MPIRequests(iMPIRequest));
          AllMPIRequestToRequest.Append(iRequest, iMPIRequest);
        }
      }
    }
  }

  while (true) {
    int iCompleted;
    StartMPITime();
    MPI_Waitany(TotalMPIRequests, AllMPIRequests.Data(), &iCompleted, MPI_STATUSES_IGNORE);
    StopMPITime();
    if (iCompleted == MPI_UNDEFINED) {
      Index = -1;
      break;
    }
    int iRequest, iMPIRequest;
    std::tie(iRequest, iMPIRequest) = AllMPIRequestToRequest(iCompleted);
    Requests(iRequest)->MPIRequests_()(iMPIRequest) = MPI_REQUEST_NULL;
    Requests(iRequest)->Finish_(iMPIRequest);
    --NumRemainingMPIRequests(iRequest);
    if (NumRemainingMPIRequests(iRequest) == 0) {
      Requests(iRequest)->Wait();
      Index = iRequest;
      break;
    }
  }

  for (auto Request : Requests) {
    if (Request && *Request) {
      Request->StopWaitTime_();
    }
  }

}

}
