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

    bool Grid1IsLocal = Domain.GridIsLocal(1);
    bool Grid2IsLocal = Domain.GridIsLocal(2);

    ovk::tuple<int> Grid1Size = Domain.GridInfo(1).Size();

    ovk::exchanger Exchanger = ovk::CreateExchanger(Domain.SharedContext());

    Exchanger.Bind(Domain, ovk::exchanger::bindings()
      .SetConnectivityComponentID(1)
    );

    using field_values = ovk::array<double,ovk::MAX_DIMS,ovk::array_layout::COLUMN_MAJOR>;

    field_values Grid1FieldValues;
    if (Grid1IsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      Grid1FieldValues.Resize(LocalRange, 0.);
      for (int j = LocalRange.Begin(1); j < ovk::Min(LocalRange.End(1),Grid1Size(1)-1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(j);
          Grid1FieldValues(i,j,0) = U*V;
        }
      }
    }

    field_values Grid2FieldValues;
    if (Grid2IsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      Grid2FieldValues.Resize(LocalRange, 0.);
      for (int j = ovk::Max(LocalRange.Begin(1), 1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(Grid1Size(1)-2+j);
          Grid2FieldValues(i,j,0) = U*V;
        }
      }
    }

    field_values ExpectedGrid1FieldValues;
    if (Grid1IsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      ExpectedGrid1FieldValues.Resize(LocalRange);
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(j);
          ExpectedGrid1FieldValues(i,j,0) = U*V;
        }
      }
      // Sanity check
      if (LocalRange.End(1) == Grid1Size(1)) {
        EXPECT_THAT(Grid1FieldValues, Not(ElementsAreArray(ExpectedGrid1FieldValues)));
      }
    }

    field_values ExpectedGrid2FieldValues;
    if (Grid2IsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      ExpectedGrid2FieldValues.Resize(LocalRange);
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(Grid1Size(1)-2+j);
          ExpectedGrid2FieldValues(i,j,0) = U*V;
        }
      }
      // Sanity check
      if (LocalRange.Begin(1) == 0) {
        EXPECT_THAT(Grid2FieldValues, Not(ElementsAreArray(ExpectedGrid2FieldValues)));
      }
    }

    ovk::array<double> Grid1DonorValues, Grid1ReceiverValues;
    ovk::array<double> ExpectedGrid1DonorValues, ExpectedGrid1ReceiverValues;
    if (Grid1IsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      if (LocalRange.End(1) == Grid1Size(1)) {
        Grid1DonorValues.Resize({LocalRange.Size(0)}, 0.);
        Grid1ReceiverValues.Resize({LocalRange.Size(0)}, 0.);
        ExpectedGrid1DonorValues.Resize({LocalRange.Size(0)});
        long long iDonor = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(Grid1Size(1)-2);
          ExpectedGrid1DonorValues(iDonor) = U*V;
          ++iDonor;
        }
        ExpectedGrid1ReceiverValues.Resize({LocalRange.Size(0)});
        long long iReceiver = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(Grid1Size(1)-1);
          ExpectedGrid1ReceiverValues(iReceiver) = U*V;
          ++iReceiver;
        }
      }
    }

    ovk::array<double> Grid2DonorValues, Grid2ReceiverValues;
    ovk::array<double> ExpectedGrid2DonorValues, ExpectedGrid2ReceiverValues;
    if (Grid2IsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      if (LocalRange.Begin(1) == 0) {
        Grid2DonorValues.Resize({LocalRange.Size(0)}, 0.);
        Grid2ReceiverValues.Resize({LocalRange.Size(0)}, 0.);
        ExpectedGrid2DonorValues.Resize({LocalRange.Size(0)});
        long long iDonor = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(Grid1Size(1)-1);
          ExpectedGrid2DonorValues(iDonor) = U*V;
          ++iDonor;
        }
        ExpectedGrid2ReceiverValues.Resize({LocalRange.Size(0)});
        long long iReceiver = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i);
          double V = double(Grid1Size(1)-2);
          ExpectedGrid2ReceiverValues(iReceiver) = U*V;
          ++iReceiver;
        }
      }
    }

    if (Grid1IsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      Exchanger.CreateCollect(1, 2, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
      Exchanger.CreateSend(1, 2, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateReceive(2, 1, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateDisperse(2, 1, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
    }

    if (Grid2IsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      Exchanger.CreateCollect(2, 1, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
      Exchanger.CreateSend(2, 1, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateReceive(1, 2, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateDisperse(1, 2, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
    }

    if (Grid1IsLocal) {
      const double *FieldValues = Grid1FieldValues.Data();
      double *DonorValues = Grid1DonorValues.Data();
      Exchanger.Collect(1, 2, 1, &FieldValues, &DonorValues);
      EXPECT_THAT(Grid1DonorValues, ElementsAreArray(ExpectedGrid1DonorValues));
    }

    if (Grid2IsLocal) {
      const double *FieldValues = Grid2FieldValues.Data();
      double *DonorValues = Grid2DonorValues.Data();
      Exchanger.Collect(2, 1, 1, &FieldValues, &DonorValues);
      EXPECT_THAT(Grid2DonorValues, ElementsAreArray(ExpectedGrid2DonorValues));
    }

    ovk::array<ovk::request> Requests;

    if (Grid1IsLocal) {
      double *ReceiverValues = Grid1ReceiverValues.Data();
      ovk::request Request = Exchanger.Receive(2, 1, 1, &ReceiverValues);
      Requests.Append(std::move(Request));
    }

    if (Grid2IsLocal) {
      double *ReceiverValues = Grid2ReceiverValues.Data();
      ovk::request Request = Exchanger.Receive(1, 2, 1, &ReceiverValues);
      Requests.Append(std::move(Request));
    }

    if (Grid1IsLocal) {
      const double *DonorValues = Grid1DonorValues.Data();
      ovk::request Request = Exchanger.Send(1, 2, 1, &DonorValues);
      Requests.Append(std::move(Request));
    }

    if (Grid2IsLocal) {
      const double *DonorValues = Grid2DonorValues.Data();
      ovk::request Request = Exchanger.Send(2, 1, 1, &DonorValues);
      Requests.Append(std::move(Request));
    }

    ovk::WaitAll(Requests);

    if (Grid1IsLocal) {
      EXPECT_THAT(Grid1ReceiverValues, ElementsAreArray(ExpectedGrid1ReceiverValues));
      const double *ReceiverValues = Grid1ReceiverValues.Data();
      double *FieldValues = Grid1FieldValues.Data();
      Exchanger.Disperse(2, 1, 1, &ReceiverValues, &FieldValues);
      EXPECT_THAT(Grid1FieldValues, ElementsAreArray(ExpectedGrid1FieldValues));
    }

    if (Grid2IsLocal) {
      EXPECT_THAT(Grid2ReceiverValues, ElementsAreArray(ExpectedGrid2ReceiverValues));
      const double *ReceiverValues = Grid2ReceiverValues.Data();
      double *FieldValues = Grid2FieldValues.Data();
      Exchanger.Disperse(1, 2, 1, &ReceiverValues, &FieldValues);
      EXPECT_THAT(Grid2FieldValues, ElementsAreArray(ExpectedGrid2FieldValues));
    }

  }

}

TEST_F(ExchangerTests, Exchange3D) {

  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < 64);

  if (Comm) {

    ovk::tuple<int> Size = {32,32,32};

    ovk::domain Domain = Interface3DManualConnectivity(Comm, {{-1.,-1.,-1}, {1.,1.,1.}}, Size,
      {false, false, false}, ovk::periodic_storage::UNIQUE);

    bool Grid1IsLocal = Domain.GridIsLocal(1);
    bool Grid2IsLocal = Domain.GridIsLocal(2);

    ovk::tuple<int> Grid1Size = Domain.GridInfo(1).Size();

    ovk::exchanger Exchanger = ovk::CreateExchanger(Domain.SharedContext());

    Exchanger.Bind(Domain, ovk::exchanger::bindings()
      .SetConnectivityComponentID(1)
    );

    using field_values = ovk::array<double,ovk::MAX_DIMS,ovk::array_layout::COLUMN_MAJOR>;

    field_values Grid1FieldValues;
    if (Grid1IsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      Grid1FieldValues.Resize(LocalRange, 0.);
      for (int k = LocalRange.Begin(2); k < ovk::Min(LocalRange.End(2),Grid1Size(2)-1); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(k);
            Grid1FieldValues(i,j,k) = U*V*W;
          }
        }
      }
    }

    field_values Grid2FieldValues;
    if (Grid2IsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      Grid2FieldValues.Resize(LocalRange, 0.);
      for (int k = ovk::Max(LocalRange.Begin(2), 1); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(Grid1Size(2)-2+k);
            Grid2FieldValues(i,j,k) = U*V*W;
          }
        }
      }
    }

    field_values ExpectedGrid1FieldValues;
    if (Grid1IsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      ExpectedGrid1FieldValues.Resize(LocalRange);
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(k);
            ExpectedGrid1FieldValues(i,j,k) = U*V*W;
          }
        }
      }
      // Sanity check
      if (LocalRange.End(2) == Grid1Size(2)) {
        EXPECT_THAT(Grid1FieldValues, Not(ElementsAreArray(ExpectedGrid1FieldValues)));
      }
    }

    field_values ExpectedGrid2FieldValues;
    if (Grid2IsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      ExpectedGrid2FieldValues.Resize(LocalRange);
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(Grid1Size(2)-2+k);
            ExpectedGrid2FieldValues(i,j,k) = U*V*W;
          }
        }
      }
      // Sanity check
      if (LocalRange.Begin(2) == 0) {
        EXPECT_THAT(Grid2FieldValues, Not(ElementsAreArray(ExpectedGrid2FieldValues)));
      }
    }

    ovk::array<double> Grid1DonorValues, Grid1ReceiverValues;
    ovk::array<double> ExpectedGrid1DonorValues, ExpectedGrid1ReceiverValues;
    if (Grid1IsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      if (LocalRange.End(2) == Grid1Size(2)) {
        Grid1DonorValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)}, 0.);
        Grid1ReceiverValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)}, 0.);
        ExpectedGrid1DonorValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)});
        long long iDonor = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(Grid1Size(2)-2);
            ExpectedGrid1DonorValues(iDonor) = U*V*W;
            ++iDonor;
          }
        }
        ExpectedGrid1ReceiverValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)});
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(Grid1Size(2)-1);
            ExpectedGrid1ReceiverValues(iReceiver) = U*V*W;
            ++iReceiver;
          }
        }
      }
    }

    ovk::array<double> Grid2DonorValues, Grid2ReceiverValues;
    ovk::array<double> ExpectedGrid2DonorValues, ExpectedGrid2ReceiverValues;
    if (Grid2IsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      if (LocalRange.Begin(2) == 0) {
        Grid2DonorValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)}, 0.);
        Grid2ReceiverValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)}, 0.);
        ExpectedGrid2DonorValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)});
        long long iDonor = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(Grid1Size(2)-1);
            ExpectedGrid2DonorValues(iDonor) = U*V*W;
            ++iDonor;
          }
        }
        ExpectedGrid2ReceiverValues.Resize({LocalRange.Size(0)*LocalRange.Size(1)});
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i);
            double V = double(j);
            double W = double(Grid1Size(2)-2);
            ExpectedGrid2ReceiverValues(iReceiver) = U*V*W;
            ++iReceiver;
          }
        }
      }
    }

    if (Grid1IsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      Exchanger.CreateCollect(1, 2, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
      Exchanger.CreateSend(1, 2, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateReceive(2, 1, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateDisperse(2, 1, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
    }

    if (Grid2IsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      Exchanger.CreateCollect(2, 1, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
      Exchanger.CreateSend(2, 1, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateReceive(1, 2, 1, ovk::data_type::DOUBLE, 1, 1);
      Exchanger.CreateDisperse(1, 2, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
        LocalRange, ovk::array_layout::COLUMN_MAJOR);
    }

    if (Grid1IsLocal) {
      const double *FieldValues = Grid1FieldValues.Data();
      double *DonorValues = Grid1DonorValues.Data();
      Exchanger.Collect(1, 2, 1, &FieldValues, &DonorValues);
      EXPECT_THAT(Grid1DonorValues, ElementsAreArray(ExpectedGrid1DonorValues));
    }

    if (Grid2IsLocal) {
      const double *FieldValues = Grid2FieldValues.Data();
      double *DonorValues = Grid2DonorValues.Data();
      Exchanger.Collect(2, 1, 1, &FieldValues, &DonorValues);
      EXPECT_THAT(Grid2DonorValues, ElementsAreArray(ExpectedGrid2DonorValues));
    }

    ovk::array<ovk::request> Requests;

    if (Grid1IsLocal) {
      double *ReceiverValues = Grid1ReceiverValues.Data();
      ovk::request Request = Exchanger.Receive(2, 1, 1, &ReceiverValues);
      Requests.Append(std::move(Request));
    }

    if (Grid2IsLocal) {
      double *ReceiverValues = Grid2ReceiverValues.Data();
      ovk::request Request = Exchanger.Receive(1, 2, 1, &ReceiverValues);
      Requests.Append(std::move(Request));
    }

    if (Grid1IsLocal) {
      const double *DonorValues = Grid1DonorValues.Data();
      ovk::request Request = Exchanger.Send(1, 2, 1, &DonorValues);
      Requests.Append(std::move(Request));
    }

    if (Grid2IsLocal) {
      const double *DonorValues = Grid2DonorValues.Data();
      ovk::request Request = Exchanger.Send(2, 1, 1, &DonorValues);
      Requests.Append(std::move(Request));
    }

    ovk::WaitAll(Requests);

    if (Grid1IsLocal) {
      EXPECT_THAT(Grid1ReceiverValues, ElementsAreArray(ExpectedGrid1ReceiverValues));
      const double *ReceiverValues = Grid1ReceiverValues.Data();
      double *FieldValues = Grid1FieldValues.Data();
      Exchanger.Disperse(2, 1, 1, &ReceiverValues, &FieldValues);
      EXPECT_THAT(Grid1FieldValues, ElementsAreArray(ExpectedGrid1FieldValues));
    }

    if (Grid2IsLocal) {
      EXPECT_THAT(Grid2ReceiverValues, ElementsAreArray(ExpectedGrid2ReceiverValues));
      const double *ReceiverValues = Grid2ReceiverValues.Data();
      double *FieldValues = Grid2FieldValues.Data();
      Exchanger.Disperse(1, 2, 1, &ReceiverValues, &FieldValues);
      EXPECT_THAT(Grid2FieldValues, ElementsAreArray(ExpectedGrid2FieldValues));
    }

  }

}
