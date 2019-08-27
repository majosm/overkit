// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Exchanger.hpp>

#include "tests/MPITest.hpp"
#include "tests/fixtures/Interface.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/ConnectivityComponent.hpp>
#include <ovk/core/ConnectivityM.hpp>
#include <ovk/core/ConnectivityN.hpp>
#include <ovk/core/Domain.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <utility>

using testing::ElementsAreArray;
using testing::Not;

class ExchangerTests : public tests::mpi_test {};

using tests::Interface2DManualConnectivity;
using tests::Interface3DManualConnectivity;

TEST_F(ExchangerTests, Exchange2D) {

  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < 16);

  if (Comm) {

    ovk::tuple<int> Size = {32,32,1};

    ovk::domain Domain = Interface2DManualConnectivity(Comm, {{-1.,-1.,0.}, {1.,1.,0.}}, Size,
      {false, false, false}, ovk::periodic_storage::UNIQUE);

    bool LowerIsLocal = Domain.GridIsLocal(1);
    bool UpperIsLocal = Domain.GridIsLocal(2);

    ovk::tuple<int> LowerSize = Domain.GridInfo(1).GlobalRange().Size();

    ovk::exchanger Exchanger = ovk::CreateExchanger(Domain.SharedContext());

    Exchanger.Bind(Domain, ovk::exchanger::bindings()
      .SetConnectivityComponentID(1)
    );

    ovk::field<double> LowerFieldValues;
    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      LowerFieldValues.Resize(LocalRange, 0.);
      for (int j = LocalRange.Begin(1); j < ovk::Min(LocalRange.End(1),LowerSize(1)-1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(j);
          LowerFieldValues(i,j,0) = U*V;
        }
      }
    }

    ovk::field<double> UpperFieldValues;
    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      UpperFieldValues.Resize(LocalRange, 0.);
      for (int j = ovk::Max(LocalRange.Begin(1), 1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(LowerSize(1)-2+j);
          UpperFieldValues(i,j,0) = U*V;
        }
      }
    }

    ovk::field<double> ExpectedLowerFieldValues;
    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      ExpectedLowerFieldValues.Resize(LocalRange);
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(j);
          ExpectedLowerFieldValues(i,j,0) = U*V;
        }
      }
      // Sanity check
      if (LocalRange.End(1) == LowerSize(1)) {
        EXPECT_THAT(LowerFieldValues, Not(ElementsAreArray(ExpectedLowerFieldValues)));
      }
    }

    ovk::field<double> ExpectedUpperFieldValues;
    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      ExpectedUpperFieldValues.Resize(LocalRange);
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(LowerSize(1)-2+j);
          ExpectedUpperFieldValues(i,j,0) = U*V;
        }
      }
      // Sanity check
      if (LocalRange.Begin(1) == 0) {
        EXPECT_THAT(UpperFieldValues, Not(ElementsAreArray(ExpectedUpperFieldValues)));
      }
    }

    ovk::array<double> LowerDonorValues, LowerReceiverValues;
    ovk::array<double> ExpectedLowerDonorValues, ExpectedLowerReceiverValues;
    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      if (LocalRange.End(1) == LowerSize(1)) {
        LowerDonorValues.Resize({LocalRange.Size(0)}, 0.);
        LowerReceiverValues.Resize({LocalRange.Size(0)}, 0.);
        ExpectedLowerDonorValues.Resize({LocalRange.Size(0)});
        long long iDonor = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(LowerSize(1)-2);
          ExpectedLowerDonorValues(iDonor) = U*V;
          ++iDonor;
        }
        ExpectedLowerReceiverValues.Resize({LocalRange.Size(0)});
        long long iReceiver = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(LowerSize(1)-1);
          ExpectedLowerReceiverValues(iReceiver) = U*V;
          ++iReceiver;
        }
      }
    }

    ovk::array<double> UpperDonorValues, UpperReceiverValues;
    ovk::array<double> ExpectedUpperDonorValues, ExpectedUpperReceiverValues;
    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      if (LocalRange.Begin(1) == 0) {
        UpperDonorValues.Resize({LocalRange.Size(0)}, 0.);
        UpperReceiverValues.Resize({LocalRange.Size(0)}, 0.);
        ExpectedUpperDonorValues.Resize({LocalRange.Size(0)});
        long long iDonor = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(LowerSize(1)-1);
          ExpectedUpperDonorValues(iDonor) = U*V;
          ++iDonor;
        }
        ExpectedUpperReceiverValues.Resize({LocalRange.Size(0)});
        long long iReceiver = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(LowerSize(1)-2);
          ExpectedUpperReceiverValues(iReceiver) = U*V;
          ++iReceiver;
        }
      }
    }

    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      Exchanger.CreateCollect({1,2}, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
      Exchanger.CreateSend({1,2}, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateReceive({2,1}, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateDisperse({2,1}, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
    }

    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      Exchanger.CreateCollect({2,1}, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
      Exchanger.CreateSend({2,1}, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateReceive({1,2}, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateDisperse({1,2}, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
    }

    if (LowerIsLocal) {
      const double *FieldValues = LowerFieldValues.Data();
      double *DonorValues = LowerDonorValues.Data();
      Exchanger.Collect({1,2}, 1, &FieldValues, &DonorValues);
      EXPECT_THAT(LowerDonorValues, ElementsAreArray(ExpectedLowerDonorValues));
    }

    if (UpperIsLocal) {
      const double *FieldValues = UpperFieldValues.Data();
      double *DonorValues = UpperDonorValues.Data();
      Exchanger.Collect({2,1}, 1, &FieldValues, &DonorValues);
      EXPECT_THAT(UpperDonorValues, ElementsAreArray(ExpectedUpperDonorValues));
    }

    ovk::array<ovk::request> Requests;

    if (LowerIsLocal) {
      double *ReceiverValues = LowerReceiverValues.Data();
      ovk::request Request = Exchanger.Receive({2,1}, 1, &ReceiverValues);
      Requests.Append(std::move(Request));
    }

    if (UpperIsLocal) {
      double *ReceiverValues = UpperReceiverValues.Data();
      ovk::request Request = Exchanger.Receive({1,2}, 1, &ReceiverValues);
      Requests.Append(std::move(Request));
    }

    if (LowerIsLocal) {
      const double *DonorValues = LowerDonorValues.Data();
      ovk::request Request = Exchanger.Send({1,2}, 1, &DonorValues);
      Requests.Append(std::move(Request));
    }

    if (UpperIsLocal) {
      const double *DonorValues = UpperDonorValues.Data();
      ovk::request Request = Exchanger.Send({2,1}, 1, &DonorValues);
      Requests.Append(std::move(Request));
    }

    ovk::WaitAll(Requests);

    if (LowerIsLocal) {
      EXPECT_THAT(LowerReceiverValues, ElementsAreArray(ExpectedLowerReceiverValues));
      const double *ReceiverValues = LowerReceiverValues.Data();
      double *FieldValues = LowerFieldValues.Data();
      Exchanger.Disperse({2,1}, 1, &ReceiverValues, &FieldValues);
      EXPECT_THAT(LowerFieldValues, ElementsAreArray(ExpectedLowerFieldValues));
    }

    if (UpperIsLocal) {
      EXPECT_THAT(UpperReceiverValues, ElementsAreArray(ExpectedUpperReceiverValues));
      const double *ReceiverValues = UpperReceiverValues.Data();
      double *FieldValues = UpperFieldValues.Data();
      Exchanger.Disperse({1,2}, 1, &ReceiverValues, &FieldValues);
      EXPECT_THAT(UpperFieldValues, ElementsAreArray(ExpectedUpperFieldValues));
    }

  }

}

TEST_F(ExchangerTests, Exchange3D) {

  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < 64);

  if (Comm) {

    ovk::tuple<int> Size = {32,32,32};

    ovk::domain Domain = Interface3DManualConnectivity(Comm, {{-1.,-1.,-1}, {1.,1.,1.}}, Size,
      {false, false, false}, ovk::periodic_storage::UNIQUE);

    bool LowerIsLocal = Domain.GridIsLocal(1);
    bool UpperIsLocal = Domain.GridIsLocal(2);

    ovk::tuple<int> LowerSize = Domain.GridInfo(1).GlobalRange().Size();

    ovk::exchanger Exchanger = ovk::CreateExchanger(Domain.SharedContext());

    Exchanger.Bind(Domain, ovk::exchanger::bindings()
      .SetConnectivityComponentID(1)
    );

    ovk::field<double> LowerFieldValues;
    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      LowerFieldValues.Resize(LocalRange, 0.);
      for (int k = LocalRange.Begin(2); k < ovk::Min(LocalRange.End(2),LowerSize(2)-1); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(k);
            LowerFieldValues(i,j,k) = U*V*W;
          }
        }
      }
    }

    ovk::field<double> UpperFieldValues;
    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      UpperFieldValues.Resize(LocalRange, 0.);
      for (int k = ovk::Max(LocalRange.Begin(2), 1); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(LowerSize(2)-2+k);
            UpperFieldValues(i,j,k) = U*V*W;
          }
        }
      }
    }

    ovk::field<double> ExpectedLowerFieldValues;
    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      ExpectedLowerFieldValues.Resize(LocalRange);
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(k);
            ExpectedLowerFieldValues(i,j,k) = U*V*W;
          }
        }
      }
      // Sanity check
      if (LocalRange.End(2) == LowerSize(2)) {
        EXPECT_THAT(LowerFieldValues, Not(ElementsAreArray(ExpectedLowerFieldValues)));
      }
    }

    ovk::field<double> ExpectedUpperFieldValues;
    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      ExpectedUpperFieldValues.Resize(LocalRange);
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(LowerSize(2)-2+k);
            ExpectedUpperFieldValues(i,j,k) = U*V*W;
          }
        }
      }
      // Sanity check
      if (LocalRange.Begin(2) == 0) {
        EXPECT_THAT(UpperFieldValues, Not(ElementsAreArray(ExpectedUpperFieldValues)));
      }
    }

    ovk::array<double> LowerDonorValues, LowerReceiverValues;
    ovk::array<double> ExpectedLowerDonorValues, ExpectedLowerReceiverValues;
    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      if (LocalRange.End(2) == LowerSize(2)) {
        LowerDonorValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)}, 0.);
        LowerReceiverValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)}, 0.);
        ExpectedLowerDonorValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)});
        long long iDonor = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(LowerSize(2)-2);
            ExpectedLowerDonorValues(iDonor) = U*V*W;
            ++iDonor;
          }
        }
        ExpectedLowerReceiverValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)});
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(LowerSize(2)-1);
            ExpectedLowerReceiverValues(iReceiver) = U*V*W;
            ++iReceiver;
          }
        }
      }
    }

    ovk::array<double> UpperDonorValues, UpperReceiverValues;
    ovk::array<double> ExpectedUpperDonorValues, ExpectedUpperReceiverValues;
    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      if (LocalRange.Begin(2) == 0) {
        UpperDonorValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)}, 0.);
        UpperReceiverValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)}, 0.);
        ExpectedUpperDonorValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)});
        long long iDonor = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(LowerSize(2)-1);
            ExpectedUpperDonorValues(iDonor) = U*V*W;
            ++iDonor;
          }
        }
        ExpectedUpperReceiverValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)});
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(LowerSize(2)-2);
            ExpectedUpperReceiverValues(iReceiver) = U*V*W;
            ++iReceiver;
          }
        }
      }
    }

    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      Exchanger.CreateCollect({1,2}, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
      Exchanger.CreateSend({1,2}, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateReceive({2,1}, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateDisperse({2,1}, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
    }

    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      Exchanger.CreateCollect({2,1}, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
      Exchanger.CreateSend({2,1}, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateReceive({1,2}, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateDisperse({1,2}, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
    }

    if (LowerIsLocal) {
      const double *FieldValues = LowerFieldValues.Data();
      double *DonorValues = LowerDonorValues.Data();
      Exchanger.Collect({1,2}, 1, &FieldValues, &DonorValues);
      EXPECT_THAT(LowerDonorValues, ElementsAreArray(ExpectedLowerDonorValues));
    }

    if (UpperIsLocal) {
      const double *FieldValues = UpperFieldValues.Data();
      double *DonorValues = UpperDonorValues.Data();
      Exchanger.Collect({2,1}, 1, &FieldValues, &DonorValues);
      EXPECT_THAT(UpperDonorValues, ElementsAreArray(ExpectedUpperDonorValues));
    }

    ovk::array<ovk::request> Requests;

    if (LowerIsLocal) {
      double *ReceiverValues = LowerReceiverValues.Data();
      ovk::request Request = Exchanger.Receive({2,1}, 1, &ReceiverValues);
      Requests.Append(std::move(Request));
    }

    if (UpperIsLocal) {
      double *ReceiverValues = UpperReceiverValues.Data();
      ovk::request Request = Exchanger.Receive({1,2}, 1, &ReceiverValues);
      Requests.Append(std::move(Request));
    }

    if (LowerIsLocal) {
      const double *DonorValues = LowerDonorValues.Data();
      ovk::request Request = Exchanger.Send({1,2}, 1, &DonorValues);
      Requests.Append(std::move(Request));
    }

    if (UpperIsLocal) {
      const double *DonorValues = UpperDonorValues.Data();
      ovk::request Request = Exchanger.Send({2,1}, 1, &DonorValues);
      Requests.Append(std::move(Request));
    }

    ovk::WaitAll(Requests);

    if (LowerIsLocal) {
      EXPECT_THAT(LowerReceiverValues, ElementsAreArray(ExpectedLowerReceiverValues));
      const double *ReceiverValues = LowerReceiverValues.Data();
      double *FieldValues = LowerFieldValues.Data();
      Exchanger.Disperse({2,1}, 1, &ReceiverValues, &FieldValues);
      EXPECT_THAT(LowerFieldValues, ElementsAreArray(ExpectedLowerFieldValues));
    }

    if (UpperIsLocal) {
      EXPECT_THAT(UpperReceiverValues, ElementsAreArray(ExpectedUpperReceiverValues));
      const double *ReceiverValues = UpperReceiverValues.Data();
      double *FieldValues = UpperFieldValues.Data();
      Exchanger.Disperse({1,2}, 1, &ReceiverValues, &FieldValues);
      EXPECT_THAT(UpperFieldValues, ElementsAreArray(ExpectedUpperFieldValues));
    }

  }

}
