// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Request.hpp"

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

void request::Wait() {

  if (*this) {

    StartWaitTime_();

    array_view<MPI_Request> MPIRequests = MPIRequests_();

    while (true) {
      int iMPIRequest;
      StartMPITime_();
      MPI_Waitany(MPIRequests.Count(), MPIRequests.Data(), &iMPIRequest, MPI_STATUSES_IGNORE);
      StopMPITime_();
      if (iMPIRequest == MPI_UNDEFINED) {
        break;
      }
      OnMPIRequestComplete_(iMPIRequest);
    }

    OnComplete_();

    StopWaitTime_();

    Reset_();

  }

}

void WaitAll(array_view<request> Requests) {

  request::internal_WaitAll(Requests);

}

void WaitAny(array_view<request> Requests, int &Index) {

  request::internal_WaitAny(Requests, Index);

}

void WaitAll(array_view<request *> Requests) {

  request::internal_WaitAll(Requests);

}

void WaitAny(array_view<request *> Requests, int &Index) {

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
        int NumMPIRequests = 0;
        for (auto &MPIRequest : Request.MPIRequests_()) {
          if (MPIRequest != MPI_REQUEST_NULL) {
            ++NumMPIRequests;
          }
        }
        NumRemainingMPIRequests(iRequest) = NumMPIRequests;
        TotalMPIRequests += NumMPIRequests;
      }
    }
  }

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests(iRequest)) {
      request &Request = *Requests(iRequest);
      if (Request && NumRemainingMPIRequests(iRequest) == 0) {
        Request.OnComplete_();
        Request.StopWaitTime_();
        Request.Reset_();
      }
    }
  }

  array<MPI_Request> AllMPIRequests;
  array<elem<int,2>> AllMPIRequestToRequest;

  AllMPIRequests.Reserve(TotalMPIRequests);
  AllMPIRequestToRequest.Reserve(TotalMPIRequests);

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests(iRequest)) {
      request &Request = *Requests(iRequest);
      if (Request) {
        array_view<MPI_Request> MPIRequests = Request.MPIRequests_();
        for (int iMPIRequest = 0; iMPIRequest < MPIRequests.Count(); ++iMPIRequest) {
          MPI_Request MPIRequest = MPIRequests(iMPIRequest);
          if (MPIRequest != MPI_REQUEST_NULL) {
            AllMPIRequests.Append(MPIRequests(iMPIRequest));
            AllMPIRequestToRequest.Append({iRequest,iMPIRequest});
          }
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
    int iRequest = AllMPIRequestToRequest(iCompleted)(0);
    int iMPIRequest = AllMPIRequestToRequest(iCompleted)(1);
    request &Request = *Requests(iRequest);
    Request.MPIRequests_()(iMPIRequest) = MPI_REQUEST_NULL;
    Request.OnMPIRequestComplete_(iMPIRequest);
    --NumRemainingMPIRequests(iRequest);
    if (NumRemainingMPIRequests(iRequest) == 0) {
      Request.OnComplete_();
      Request.StopWaitTime_();
      Request.Reset_();
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
        int NumMPIRequests = 0;
        for (auto &MPIRequest : Request.MPIRequests_()) {
          if (MPIRequest != MPI_REQUEST_NULL) {
            ++NumMPIRequests;
          }
        }
        NumRemainingMPIRequests(iRequest) = NumMPIRequests;
        TotalMPIRequests += NumMPIRequests;
      }
    }
  }

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests(iRequest)) {
      request &Request = *Requests(iRequest);
      if (Request && NumRemainingMPIRequests(iRequest) == 0) {
        Request.OnComplete_();
        Request.StopWaitTime_();
        Request.Reset_();
        Index = iRequest;
        return;
      }
    }
  }

  array<MPI_Request> AllMPIRequests;
  array<elem<int,2>> AllMPIRequestToRequest;

  AllMPIRequests.Reserve(TotalMPIRequests);
  AllMPIRequestToRequest.Reserve(TotalMPIRequests);

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests(iRequest)) {
      request &Request = *Requests(iRequest);
      if (Request) {
        array_view<MPI_Request> MPIRequests = Request.MPIRequests_();
        for (int iMPIRequest = 0; iMPIRequest < MPIRequests.Count(); ++iMPIRequest) {
          MPI_Request MPIRequest = MPIRequests(iMPIRequest);
          if (MPIRequest != MPI_REQUEST_NULL) {
            AllMPIRequests.Append(MPIRequests(iMPIRequest));
            AllMPIRequestToRequest.Append({iRequest,iMPIRequest});
          }
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
    int iRequest = AllMPIRequestToRequest(iCompleted)(0);
    int iMPIRequest = AllMPIRequestToRequest(iCompleted)(1);
    request &Request = *Requests(iRequest);
    Request.MPIRequests_()(iMPIRequest) = MPI_REQUEST_NULL;
    Request.OnMPIRequestComplete_(iMPIRequest);
    --NumRemainingMPIRequests(iRequest);
    if (NumRemainingMPIRequests(iRequest) == 0) {
      Request.OnComplete_();
      Request.StopWaitTime_();
      Request.Reset_();
      Index = iRequest;
      return;
    }
  }

}

}
