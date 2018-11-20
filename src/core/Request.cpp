// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Request.hpp"

#include "ovk/core/Global.hpp"

#include <mpi.h>

#include <memory>
#include <utility>
#include <vector>

namespace ovk {

void WaitAll(int NumRequests, request *Requests) {

  std::vector<request *> RequestPtrs(NumRequests);

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    RequestPtrs[iRequest] = &Requests[iRequest];
  }

  request::core_WaitAll(NumRequests, RequestPtrs.data());

}

void WaitAny(int NumRequests, request *Requests, int &Index) {

  std::vector<request *> RequestPtrs(NumRequests);

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    RequestPtrs[iRequest] = &Requests[iRequest];
  }

  request::core_WaitAny(NumRequests, RequestPtrs.data(), Index);

}

void request::core_WaitAll(int NumRequests, request **Requests) {

  std::vector<int> NumRemainingMPIRequests(NumRequests);
  int TotalMPIRequests = 0;

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests[iRequest]) {
      request &Request = *Requests[iRequest];
      if (Request) {
        int NumMPIRequests = Request.NumMPIRequests();
        NumRemainingMPIRequests[iRequest] = NumMPIRequests;
        TotalMPIRequests += NumMPIRequests;
      }
    }
  }

  std::vector<MPI_Request> AllMPIRequests(TotalMPIRequests);
  std::vector<int> MPIRequestToRequest(TotalMPIRequests);

  int iNextMPIRequest = 0;
  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests[iRequest]) {
      request &Request = *Requests[iRequest];
      if (Request) {
        int NumMPIRequests = Request.NumMPIRequests();
        MPI_Request *MPIRequests = Request.MPIRequests();
        for (int iMPIRequest = 0; iMPIRequest < NumMPIRequests; ++iMPIRequest) {
          AllMPIRequests[iNextMPIRequest] = MPIRequests[iMPIRequest];
          MPIRequestToRequest[iNextMPIRequest] = iRequest;
          ++iNextMPIRequest;
        }
      }
    }
  }

  while (true) {
    int iMPIRequest;
    MPI_Waitany(TotalMPIRequests, AllMPIRequests.data(), &iMPIRequest, MPI_STATUSES_IGNORE);
    if (iMPIRequest == MPI_UNDEFINED) {
      break;
    }
    int iRequest = MPIRequestToRequest[iMPIRequest];
    --NumRemainingMPIRequests[iRequest];
    if (NumRemainingMPIRequests[iRequest] == 0) {
      request &Request = *Requests[iRequest];
      int NumMPIRequests = Request.NumMPIRequests();
      MPI_Request *MPIRequests = Request.MPIRequests();
      for (int iMPIRequest = 0; iMPIRequest < NumMPIRequests; ++iMPIRequest) {
        MPIRequests[iMPIRequest] = MPI_REQUEST_NULL;
      }
      Request.Wait();
      Request = request();
    }
  }

}

void request::core_WaitAny(int NumRequests, request **Requests, int &Index) {

  std::vector<int> NumRemainingMPIRequests(NumRequests);
  int TotalMPIRequests = 0;

  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests[iRequest]) {
      request &Request = *Requests[iRequest];
      if (Request) {
        int NumMPIRequests = Request.NumMPIRequests();
        NumRemainingMPIRequests[iRequest] = NumMPIRequests;
        TotalMPIRequests += NumMPIRequests;
      }
    }
  }

  std::vector<MPI_Request> AllMPIRequests(TotalMPIRequests);
  std::vector<int> MPIRequestToRequest(TotalMPIRequests);

  int iNextMPIRequest = 0;
  for (int iRequest = 0; iRequest < NumRequests; ++iRequest) {
    if (Requests[iRequest]) {
      request &Request = *Requests[iRequest];
      if (Request) {
        int NumMPIRequests = Request.NumMPIRequests();
        MPI_Request *MPIRequests = Request.MPIRequests();
        for (int iMPIRequest = 0; iMPIRequest < NumMPIRequests; ++iMPIRequest) {
          AllMPIRequests[iNextMPIRequest] = MPIRequests[iMPIRequest];
          MPIRequestToRequest[iNextMPIRequest] = iRequest;
          ++iNextMPIRequest;
        }
      }
    }
  }

  while (true) {
    int iMPIRequest;
    MPI_Waitany(TotalMPIRequests, AllMPIRequests.data(), &iMPIRequest, MPI_STATUSES_IGNORE);
    if (iMPIRequest == MPI_UNDEFINED) {
      Index = -1;
      break;
    }
    int iRequest = MPIRequestToRequest[iMPIRequest];
    --NumRemainingMPIRequests[iRequest];
    if (NumRemainingMPIRequests[iRequest] == 0) {
      request &Request = *Requests[iRequest];
      int NumMPIRequests = Request.NumMPIRequests();
      MPI_Request *MPIRequests = Request.MPIRequests();
      for (int iMPIRequest = 0; iMPIRequest < NumMPIRequests; ++iMPIRequest) {
        MPIRequests[iMPIRequest] = MPI_REQUEST_NULL;
      }
      Request.Wait();
      Request = request();
      Index = iRequest;
      break;
    }
  }

}

}
