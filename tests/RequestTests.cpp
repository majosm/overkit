// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Request.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Array.hpp>

#include <mpi.h>

#include <functional>
#include <utility>

class RequestTests : public tests::mpi_test {};

namespace {

class mock_request {

public:

  template <typename F1, typename F2> mock_request(ovk::array<MPI_Request> MPIRequests, F1
    OnMPIRequestComplete, F2 OnComplete):
    MPIRequests_(std::move(MPIRequests)),
    OnMPIRequestComplete_(std::move(OnMPIRequestComplete)),
    OnComplete_(std::move(OnComplete))
  {}

  ovk::array<MPI_Request> &MPIRequests() { return MPIRequests_; }

  void OnMPIRequestComplete(int iMPIRequest) {
    OnMPIRequestComplete_(iMPIRequest);
  }

  void OnComplete() {
    OnComplete_();
  }

  void StartWaitTime() const {}
  void StopWaitTime() const {}
  void StartMPITime() const {}
  void StopMPITime() const {}

private:

  ovk::array<MPI_Request> MPIRequests_;
  std::function<void(int)> OnMPIRequestComplete_;
  std::function<void()> OnComplete_;

};

int MathMod(int Value, int Period) {
  int CMod = Value % Period;
  return CMod + Period * (CMod < 0);
}

}

TEST_F(RequestTests, Wait) {

  ovk::comm_view Comm = TestComm();

  // Empty MPI requests array
  {
    ovk::request Request = mock_request({}, [](int) {}, [] {});
    Request.Wait();
    EXPECT_FALSE(Request);
  }

  // Non-empty MPI requests array
  {

    int LeftRank = MathMod(Comm.Rank()-1, Comm.Size());
    int RightRank = MathMod(Comm.Rank()+1, Comm.Size());

    ovk::array<MPI_Request> MPIRequests;

    int ReceivedValue = 0;
    MPI_Irecv(&ReceivedValue, 1, MPI_INT, LeftRank, 0, Comm, &MPIRequests.Append());

    int SentValue = Comm.Rank();
    MPI_Isend(&SentValue, 1, MPI_INT, RightRank, 0, Comm, &MPIRequests.Append());

    auto OnMPIRequestComplete = [&](int iMPIRequest) {
      if (iMPIRequest == 0) {
        ReceivedValue *= 2;
      } else {
        SentValue += 1;
      }
    };

    auto OnComplete = [&] {
      ReceivedValue *= 2;
      SentValue += 1;
    };

    ovk::request Request = mock_request(MPIRequests, OnMPIRequestComplete, OnComplete);

    Request.Wait();

    EXPECT_FALSE(Request);
    EXPECT_EQ(ReceivedValue, 4*(MathMod(Comm.Rank()-1, Comm.Size())));
    EXPECT_EQ(SentValue, Comm.Rank()+2);

  }

}

TEST_F(RequestTests, WaitAll) {

  ovk::comm_view Comm = TestComm();

  // Empty MPI requests arrays
  {
    ovk::array<ovk::request> Requests;
    Requests.Append(mock_request({}, [](int) {}, [] {}));
    Requests.Append(mock_request({}, [](int) {}, [] {}));
    ovk::WaitAll(Requests);
    EXPECT_FALSE(Requests(0));
    EXPECT_FALSE(Requests(1));
  }

  // Non-empty MPI requests arrays
  {

    int LeftRank = MathMod(Comm.Rank()-1, Comm.Size());
    int RightRank = MathMod(Comm.Rank()+1, Comm.Size());

    int ReceivedValueLeftToRight = 0;
    int SentValueLeftToRight = Comm.Rank();

    int ReceivedValueRightToLeft = 0;
    int SentValueRightToLeft = Comm.Rank();

    ovk::array<MPI_Request> MPIRequestsLeftToRight;
    ovk::array<MPI_Request> MPIRequestsRightToLeft;

    MPI_Irecv(&ReceivedValueLeftToRight, 1, MPI_INT, LeftRank, 0, Comm, &MPIRequestsLeftToRight.
      Append());
    MPI_Isend(&SentValueLeftToRight, 1, MPI_INT, RightRank, 0, Comm, &MPIRequestsLeftToRight.
      Append());

    MPI_Irecv(&ReceivedValueRightToLeft, 1, MPI_INT, RightRank, 0, Comm, &MPIRequestsRightToLeft.
      Append());
    MPI_Isend(&SentValueRightToLeft, 1, MPI_INT, LeftRank, 0, Comm, &MPIRequestsRightToLeft.
      Append());

    auto OnMPIRequestCompleteLeftToRight = [&](int iMPIRequest) {
      if (iMPIRequest == 0) {
        ReceivedValueLeftToRight *= 2;
      } else {
        SentValueLeftToRight += 1;
      }
    };

    auto OnCompleteLeftToRight = [&] {
      ReceivedValueLeftToRight *= 2;
      SentValueLeftToRight += 1;
    };

    auto OnMPIRequestCompleteRightToLeft = [&](int iMPIRequest) {
      if (iMPIRequest == 0) {
        ReceivedValueRightToLeft *= 2;
      } else {
        SentValueRightToLeft += 1;
      }
    };

    auto OnCompleteRightToLeft = [&] {
      ReceivedValueRightToLeft *= 2;
      SentValueRightToLeft += 1;
    };

    ovk::array<ovk::request> Requests;
    Requests.Append(mock_request(MPIRequestsLeftToRight, OnMPIRequestCompleteLeftToRight,
      OnCompleteLeftToRight));
    Requests.Append(mock_request(MPIRequestsRightToLeft, OnMPIRequestCompleteRightToLeft,
      OnCompleteRightToLeft));

    ovk::WaitAll(Requests);

    EXPECT_FALSE(Requests(0));
    EXPECT_FALSE(Requests(1));
    EXPECT_EQ(ReceivedValueLeftToRight, 4*(MathMod(Comm.Rank()-1, Comm.Size())));
    EXPECT_EQ(ReceivedValueRightToLeft, 4*(MathMod(Comm.Rank()+1, Comm.Size())));
    EXPECT_EQ(SentValueLeftToRight, Comm.Rank()+2);
    EXPECT_EQ(SentValueRightToLeft, Comm.Rank()+2);

  }

}

TEST_F(RequestTests, WaitAny) {

  ovk::comm_view Comm = TestComm();

  // Empty MPI requests arrays
  {

    ovk::array<ovk::request> Requests;
    Requests.Append(mock_request({}, [](int) {}, [] {}));
    Requests.Append(mock_request({}, [](int) {}, [] {}));

    int iRequest1;
    ovk::WaitAny(Requests, iRequest1);

    EXPECT_GE(iRequest1, 0);
    EXPECT_LT(iRequest1, 2);
    EXPECT_FALSE(Requests(iRequest1));

    int iRequest2;
    ovk::WaitAny(Requests, iRequest2);

    EXPECT_GE(iRequest2, 0);
    EXPECT_LT(iRequest2, 2);
    EXPECT_NE(iRequest2, iRequest1);
    EXPECT_FALSE(Requests(iRequest2));

    int iRequest3;
    ovk::WaitAny(Requests, iRequest3);

    EXPECT_EQ(iRequest3, -1);

  }

  // Non-empty MPI requests arrays
  {

    int LeftRank = MathMod(Comm.Rank()-1, Comm.Size());
    int RightRank = MathMod(Comm.Rank()+1, Comm.Size());

    int ReceivedValueLeftToRight = 0;
    int SentValueLeftToRight = Comm.Rank();

    int ReceivedValueRightToLeft = 0;
    int SentValueRightToLeft = Comm.Rank();

    ovk::array<MPI_Request> MPIRequestsLeftToRight;
    ovk::array<MPI_Request> MPIRequestsRightToLeft;

    MPI_Irecv(&ReceivedValueLeftToRight, 1, MPI_INT, LeftRank, 0, Comm, &MPIRequestsLeftToRight.
      Append());
    MPI_Isend(&SentValueLeftToRight, 1, MPI_INT, RightRank, 0, Comm, &MPIRequestsLeftToRight.
      Append());

    MPI_Irecv(&ReceivedValueRightToLeft, 1, MPI_INT, RightRank, 0, Comm, &MPIRequestsRightToLeft.
      Append());
    MPI_Isend(&SentValueRightToLeft, 1, MPI_INT, LeftRank, 0, Comm, &MPIRequestsRightToLeft.
      Append());

    auto OnMPIRequestCompleteLeftToRight = [&](int iMPIRequest) {
      if (iMPIRequest == 0) {
        ReceivedValueLeftToRight *= 2;
      } else {
        SentValueLeftToRight += 1;
      }
    };

    auto OnCompleteLeftToRight = [&] {
      ReceivedValueLeftToRight *= 2;
      SentValueLeftToRight += 1;
    };

    auto OnMPIRequestCompleteRightToLeft = [&](int iMPIRequest) {
      if (iMPIRequest == 0) {
        ReceivedValueRightToLeft *= 2;
      } else {
        SentValueRightToLeft += 1;
      }
    };

    auto OnCompleteRightToLeft = [&] {
      ReceivedValueRightToLeft *= 2;
      SentValueRightToLeft += 1;
    };

    ovk::array<ovk::request> Requests;
    Requests.Append(mock_request(MPIRequestsLeftToRight, OnMPIRequestCompleteLeftToRight,
      OnCompleteLeftToRight));
    Requests.Append(mock_request(MPIRequestsRightToLeft, OnMPIRequestCompleteRightToLeft,
      OnCompleteRightToLeft));

    int iRequest1;
    ovk::WaitAny(Requests, iRequest1);

    EXPECT_GE(iRequest1, 0);
    EXPECT_LT(iRequest1, 2);
    EXPECT_FALSE(Requests(iRequest1));
    if (iRequest1 == 0) {
      EXPECT_EQ(ReceivedValueLeftToRight, 4*(MathMod(Comm.Rank()-1, Comm.Size())));
      EXPECT_EQ(SentValueLeftToRight, Comm.Rank()+2);
    } else {
      EXPECT_EQ(ReceivedValueRightToLeft, 4*(MathMod(Comm.Rank()+1, Comm.Size())));
      EXPECT_EQ(SentValueRightToLeft, Comm.Rank()+2);
    }

    int iRequest2;
    ovk::WaitAny(Requests, iRequest2);

    EXPECT_GE(iRequest2, 0);
    EXPECT_LT(iRequest2, 2);
    EXPECT_NE(iRequest2, iRequest1);
    EXPECT_FALSE(Requests(iRequest2));
    if (iRequest2 == 0) {
      EXPECT_EQ(ReceivedValueLeftToRight, 4*(MathMod(Comm.Rank()-1, Comm.Size())));
      EXPECT_EQ(SentValueLeftToRight, Comm.Rank()+2);
    } else {
      EXPECT_EQ(ReceivedValueRightToLeft, 4*(MathMod(Comm.Rank()+1, Comm.Size())));
      EXPECT_EQ(SentValueRightToLeft, Comm.Rank()+2);
    }

    int iRequest3;
    ovk::WaitAny(Requests, iRequest3);

    EXPECT_EQ(iRequest3, -1);

  }

}
