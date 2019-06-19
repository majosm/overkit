// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/extras/XINTOUT.hpp>

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
#include <ovk/core/Exchanger.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

using testing::ElementsAre;
using testing::ElementsAreArray;

class XINTOUTTests : public tests::mpi_test {};

using tests::Interface2D;
using tests::Interface3D;

TEST_F(XINTOUTTests, ImportStandard2D) {

  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < 16);

  if (Comm) {

    ovk::domain Domain = Interface2D(Comm, {{-1.,-1.,0.}, {1.,1.,0.}}, {32,32,1}, {false, false,
      false}, ovk::periodic_storage::UNIQUE);

    Domain.CreateComponent<ovk::connectivity_component>(1);

    ovk::ImportXINTOUT(Domain, 1, "data/XINTOUTTests/XINTOUT_Standard2D.HO.2D",
      "data/XINTOUTTests/XINTOUT_Standard2D.X.2D", 0, MPI_INFO_NULL);

    bool Grid1IsLocal = Domain.GridIsLocal(1);
    bool Grid2IsLocal = Domain.GridIsLocal(2);

    ovk::exchanger Exchanger = ovk::CreateExchanger(Domain.SharedContext());

    Exchanger.Bind(Domain, ovk::exchanger::bindings()
      .SetConnectivityComponentID(1)
    );

    auto &ConnectivityComponent = Domain.Component<ovk::connectivity_component>(1);

    if (Grid1IsLocal) {

      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.End(1) == GlobalRange.End(1);

      const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(1,2);
      const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(2,1);

      if (HasInterface) {

        long long NumDonors = ConnectivityM.Count();
        EXPECT_EQ(NumDonors, LocalRange.Size(0));
        const ovk::array<int,3> &Extents = ConnectivityM.Extents();
        const ovk::array<double,2> &Coords = ConnectivityM.Coords();
        const ovk::array<double,3> &InterpCoefs = ConnectivityM.InterpCoefs();
        const ovk::array<int,2> &Destinations = ConnectivityM.Destinations();
        ovk::array<int,3> ExpectedExtents({{2,ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,2> ExpectedCoords({{ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,3> ExpectedInterpCoefs({{ovk::MAX_DIMS,1,NumDonors}});
        ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumDonors}});
        long long iDonor = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedExtents(0,0,iDonor) = i;
          ExpectedExtents(0,1,iDonor) = GlobalRange.End(1)-2;
          ExpectedExtents(0,2,iDonor) = 0;
          ExpectedExtents(1,0,iDonor) = ExpectedExtents(0,0,iDonor)+1;
          ExpectedExtents(1,1,iDonor) = ExpectedExtents(0,1,iDonor)+1;
          ExpectedExtents(1,2,iDonor) = 1;
          ExpectedCoords(0,iDonor) = 0.;
          ExpectedCoords(1,iDonor) = 0.;
          ExpectedCoords(2,iDonor) = 0.;
          ExpectedInterpCoefs(0,0,iDonor) = 1.;
          ExpectedInterpCoefs(1,0,iDonor) = 1.;
          ExpectedInterpCoefs(2,0,iDonor) = 1.;
          ExpectedDestinations(0,iDonor) = i;
          ExpectedDestinations(1,iDonor) = 0;
          ExpectedDestinations(2,iDonor) = 0;
          ++iDonor;
        }
        EXPECT_THAT(Extents, ElementsAreArray(ExpectedExtents));
        EXPECT_THAT(Coords, ElementsAreArray(ExpectedCoords));
        EXPECT_THAT(InterpCoefs, ElementsAreArray(ExpectedInterpCoefs));
        EXPECT_THAT(Destinations, ElementsAreArray(ExpectedDestinations));

        long long NumReceivers = ConnectivityN.Count();
        EXPECT_EQ(NumReceivers, LocalRange.Size(0));
        const ovk::array<int,2> &Points = ConnectivityN.Points();
        const ovk::array<int,2> &Sources = ConnectivityN.Sources();
        ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumReceivers}});
        ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumReceivers}});
        long long iReceiver = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedPoints(0,iReceiver) = i;
          ExpectedPoints(1,iReceiver) = GlobalRange.End(1)-1;
          ExpectedPoints(2,iReceiver) = 0;
          ExpectedSources(0,iReceiver) = i;
          ExpectedSources(1,iReceiver) = 1;
          ExpectedSources(2,iReceiver) = 0;
          ++iReceiver;
        }
        EXPECT_THAT(Points, ElementsAreArray(Points));
        EXPECT_THAT(Sources, ElementsAreArray(Sources));

      } else {

        EXPECT_EQ(ConnectivityM.Count(), 0);
        EXPECT_EQ(ConnectivityN.Count(), 0);

      }

    }

    if (Grid2IsLocal) {

      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.Begin(1) == GlobalRange.Begin(1);

      const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(2,1);
      const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(1,2);

      if (HasInterface) {

        long long NumDonors = ConnectivityM.Count();
        EXPECT_EQ(NumDonors, LocalRange.Size(0));
        const ovk::array<int,3> &Extents = ConnectivityM.Extents();
        const ovk::array<double,2> &Coords = ConnectivityM.Coords();
        const ovk::array<double,3> &InterpCoefs = ConnectivityM.InterpCoefs();
        const ovk::array<int,2> &Destinations = ConnectivityM.Destinations();
        ovk::array<int,3> ExpectedExtents({{2,ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,2> ExpectedCoords({{ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,3> ExpectedInterpCoefs({{ovk::MAX_DIMS,1,NumDonors}});
        ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumDonors}});
        long long iDonor = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedExtents(0,0,iDonor) = i;
          ExpectedExtents(0,1,iDonor) = 1;
          ExpectedExtents(0,2,iDonor) = 0;
          ExpectedExtents(1,0,iDonor) = ExpectedExtents(0,0,iDonor)+1;
          ExpectedExtents(1,1,iDonor) = ExpectedExtents(0,1,iDonor)+1;
          ExpectedExtents(1,2,iDonor) = 1;
          ExpectedCoords(0,iDonor) = 0.;
          ExpectedCoords(1,iDonor) = 0.;
          ExpectedCoords(2,iDonor) = 0.;
          ExpectedInterpCoefs(0,0,iDonor) = 1.;
          ExpectedInterpCoefs(1,0,iDonor) = 1.;
          ExpectedInterpCoefs(2,0,iDonor) = 1.;
          ExpectedDestinations(0,iDonor) = i;
          ExpectedDestinations(1,iDonor) = Domain.GridInfo(1).GlobalRange().End(1)-1;
          ExpectedDestinations(2,iDonor) = 0;
          ++iDonor;
        }
        EXPECT_THAT(Extents, ElementsAreArray(ExpectedExtents));
        EXPECT_THAT(Coords, ElementsAreArray(ExpectedCoords));
        EXPECT_THAT(InterpCoefs, ElementsAreArray(ExpectedInterpCoefs));
        EXPECT_THAT(Destinations, ElementsAreArray(ExpectedDestinations));

        long long NumReceivers = ConnectivityN.Count();
        EXPECT_EQ(NumReceivers, LocalRange.Size(0));
        const ovk::array<int,2> &Points = ConnectivityN.Points();
        const ovk::array<int,2> &Sources = ConnectivityN.Sources();
        ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumReceivers}});
        ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumReceivers}});
        long long iReceiver = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedPoints(0,iReceiver) = i;
          ExpectedPoints(1,iReceiver) = 0;
          ExpectedPoints(2,iReceiver) = 0;
          ExpectedSources(0,iReceiver) = i;
          ExpectedSources(1,iReceiver) = Domain.GridInfo(1).GlobalRange().End(1)-2;
          ExpectedSources(2,iReceiver) = 0;
          ++iReceiver;
        }
        EXPECT_THAT(Points, ElementsAreArray(Points));
        EXPECT_THAT(Sources, ElementsAreArray(Sources));

      } else {

        EXPECT_EQ(ConnectivityM.Count(), 0);
        EXPECT_EQ(ConnectivityN.Count(), 0);

      }

    }

  }

}

TEST_F(XINTOUTTests, ImportStandard3D) {

  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < 64);

  if (Comm) {

    ovk::domain Domain = Interface3D(Comm, {{-1.,-1.,-1}, {1.,1.,1.}}, {32,32,32}, {false, false,
      false}, ovk::periodic_storage::UNIQUE);

    Domain.CreateComponent<ovk::connectivity_component>(1);

    ovk::ImportXINTOUT(Domain, 1, "data/XINTOUTTests/XINTOUT_Standard3D.HO.2D",
      "data/XINTOUTTests/XINTOUT_Standard3D.X.2D", 0, MPI_INFO_NULL);

    bool Grid1IsLocal = Domain.GridIsLocal(1);
    bool Grid2IsLocal = Domain.GridIsLocal(2);

    ovk::exchanger Exchanger = ovk::CreateExchanger(Domain.SharedContext());

    Exchanger.Bind(Domain, ovk::exchanger::bindings()
      .SetConnectivityComponentID(1)
    );

    auto &ConnectivityComponent = Domain.Component<ovk::connectivity_component>(1);

    if (Grid1IsLocal) {

      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.End(2) == GlobalRange.End(2);

      const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(1,2);
      const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(2,1);

      if (HasInterface) {

        long long NumDonors = ConnectivityM.Count();
        EXPECT_EQ(NumDonors, LocalRange.Size(0)*LocalRange.Size(1));
        const ovk::array<int,3> &Extents = ConnectivityM.Extents();
        const ovk::array<double,2> &Coords = ConnectivityM.Coords();
        const ovk::array<double,3> &InterpCoefs = ConnectivityM.InterpCoefs();
        const ovk::array<int,2> &Destinations = ConnectivityM.Destinations();
        ovk::array<int,3> ExpectedExtents({{2,ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,2> ExpectedCoords({{ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,3> ExpectedInterpCoefs({{ovk::MAX_DIMS,1,NumDonors}});
        ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumDonors}});
        long long iDonor = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ExpectedExtents(0,0,iDonor) = i;
            ExpectedExtents(0,1,iDonor) = j;
            ExpectedExtents(0,2,iDonor) = GlobalRange.End(2)-2;
            ExpectedExtents(1,0,iDonor) = ExpectedExtents(0,0,iDonor)+1;
            ExpectedExtents(1,1,iDonor) = ExpectedExtents(0,1,iDonor)+1;
            ExpectedExtents(1,2,iDonor) = ExpectedExtents(0,2,iDonor)+1;
            ExpectedCoords(0,iDonor) = 0.;
            ExpectedCoords(1,iDonor) = 0.;
            ExpectedCoords(2,iDonor) = 0.;
            ExpectedInterpCoefs(0,0,iDonor) = 1.;
            ExpectedInterpCoefs(1,0,iDonor) = 1.;
            ExpectedInterpCoefs(2,0,iDonor) = 1.;
            ExpectedDestinations(0,iDonor) = i;
            ExpectedDestinations(1,iDonor) = j;
            ExpectedDestinations(2,iDonor) = 0;
            ++iDonor;
          }
        }
        EXPECT_THAT(Extents, ElementsAreArray(ExpectedExtents));
        EXPECT_THAT(Coords, ElementsAreArray(ExpectedCoords));
        EXPECT_THAT(InterpCoefs, ElementsAreArray(ExpectedInterpCoefs));
        EXPECT_THAT(Destinations, ElementsAreArray(ExpectedDestinations));

        long long NumReceivers = ConnectivityN.Count();
        EXPECT_EQ(NumReceivers, LocalRange.Size(0)*LocalRange.Size(1));
        const ovk::array<int,2> &Points = ConnectivityN.Points();
        const ovk::array<int,2> &Sources = ConnectivityN.Sources();
        ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumReceivers}});
        ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumReceivers}});
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ExpectedPoints(0,iReceiver) = i;
            ExpectedPoints(1,iReceiver) = j;
            ExpectedPoints(2,iReceiver) = GlobalRange.End(2)-1;
            ExpectedSources(0,iReceiver) = i;
            ExpectedSources(1,iReceiver) = j;
            ExpectedSources(2,iReceiver) = 1;
            ++iReceiver;
          }
        }
        EXPECT_THAT(Points, ElementsAreArray(Points));
        EXPECT_THAT(Sources, ElementsAreArray(Sources));

      } else {

        EXPECT_EQ(ConnectivityM.Count(), 0);
        EXPECT_EQ(ConnectivityN.Count(), 0);

      }
    }

    if (Grid2IsLocal) {

      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.Begin(2) == GlobalRange.Begin(2);

      const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(2,1);
      const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(1,2);

      if (HasInterface) {

        long long NumDonors = ConnectivityM.Count();
        EXPECT_EQ(NumDonors, LocalRange.Size(0)*LocalRange.Size(1));
        const ovk::array<int,3> &Extents = ConnectivityM.Extents();
        const ovk::array<double,2> &Coords = ConnectivityM.Coords();
        const ovk::array<double,3> &InterpCoefs = ConnectivityM.InterpCoefs();
        const ovk::array<int,2> &Destinations = ConnectivityM.Destinations();
        ovk::array<int,3> ExpectedExtents({{2,ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,2> ExpectedCoords({{ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,3> ExpectedInterpCoefs({{ovk::MAX_DIMS,1,NumDonors}});
        ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumDonors}});
        long long iDonor = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ExpectedExtents(0,0,iDonor) = i;
            ExpectedExtents(0,1,iDonor) = j;
            ExpectedExtents(0,2,iDonor) = 1;
            ExpectedExtents(1,0,iDonor) = ExpectedExtents(0,0,iDonor)+1;
            ExpectedExtents(1,1,iDonor) = ExpectedExtents(0,1,iDonor)+1;
            ExpectedExtents(1,2,iDonor) = ExpectedExtents(0,2,iDonor)+1;
            ExpectedCoords(0,iDonor) = 0.;
            ExpectedCoords(1,iDonor) = 0.;
            ExpectedCoords(2,iDonor) = 0.;
            ExpectedInterpCoefs(0,0,iDonor) = 1.;
            ExpectedInterpCoefs(1,0,iDonor) = 1.;
            ExpectedInterpCoefs(2,0,iDonor) = 1.;
            ExpectedDestinations(0,iDonor) = i;
            ExpectedDestinations(1,iDonor) = j;
            ExpectedDestinations(2,iDonor) = Domain.GridInfo(1).GlobalRange().End(2)-1;
            ++iDonor;
          }
        }
        EXPECT_THAT(Extents, ElementsAreArray(ExpectedExtents));
        EXPECT_THAT(Coords, ElementsAreArray(ExpectedCoords));
        EXPECT_THAT(InterpCoefs, ElementsAreArray(ExpectedInterpCoefs));
        EXPECT_THAT(Destinations, ElementsAreArray(ExpectedDestinations));

        long long NumReceivers = ConnectivityN.Count();
        EXPECT_EQ(NumReceivers, LocalRange.Size(0)*LocalRange.Size(1));
        const ovk::array<int,2> &Points = ConnectivityN.Points();
        const ovk::array<int,2> &Sources = ConnectivityN.Sources();
        ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumReceivers}});
        ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumReceivers}});
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ExpectedPoints(0,iReceiver) = i;
            ExpectedPoints(1,iReceiver) = j;
            ExpectedPoints(2,iReceiver) = 0;
            ExpectedSources(0,iReceiver) = i;
            ExpectedSources(1,iReceiver) = j;
            ExpectedSources(2,iReceiver) = Domain.GridInfo(1).GlobalRange().End(2)-2;
            ++iReceiver;
          }
        }
        EXPECT_THAT(Points, ElementsAreArray(Points));
        EXPECT_THAT(Sources, ElementsAreArray(Sources));

      } else {

        EXPECT_EQ(ConnectivityM.Count(), 0);
        EXPECT_EQ(ConnectivityN.Count(), 0);

      }

    }

  }

}

TEST_F(XINTOUTTests, ImportExtended2D) {

  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < 16);

  if (Comm) {

    ovk::domain Domain = Interface2D(Comm, {{-1.,-1.,0.}, {1.,1.,0.}}, {32,32,1}, {false, false,
      false}, ovk::periodic_storage::UNIQUE);

    Domain.CreateComponent<ovk::connectivity_component>(1);

    ovk::ImportXINTOUT(Domain, 1, "data/XINTOUTTests/XINTOUT_Extended2D.HO.2D",
      "data/XINTOUTTests/XINTOUT_Extended2D.X.2D", 0, MPI_INFO_NULL);

    bool Grid1IsLocal = Domain.GridIsLocal(1);
    bool Grid2IsLocal = Domain.GridIsLocal(2);

    ovk::exchanger Exchanger = ovk::CreateExchanger(Domain.SharedContext());

    Exchanger.Bind(Domain, ovk::exchanger::bindings()
      .SetConnectivityComponentID(1)
    );

    auto &ConnectivityComponent = Domain.Component<ovk::connectivity_component>(1);

    if (Grid1IsLocal) {

      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.End(1) == GlobalRange.End(1);

      const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(1,2);
      const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(2,1);

      if (HasInterface) {

        long long NumDonors = ConnectivityM.Count();
        EXPECT_EQ(NumDonors, LocalRange.Size(0));
        const ovk::array<int,3> &Extents = ConnectivityM.Extents();
        const ovk::array<double,2> &Coords = ConnectivityM.Coords();
        const ovk::array<double,3> &InterpCoefs = ConnectivityM.InterpCoefs();
        const ovk::array<int,2> &Destinations = ConnectivityM.Destinations();
        ovk::array<int,3> ExpectedExtents({{2,ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,2> ExpectedCoords({{ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,3> ExpectedInterpCoefs({{ovk::MAX_DIMS,1,NumDonors}});
        ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumDonors}});
        long long iDonor = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedExtents(0,0,iDonor) = i;
          ExpectedExtents(0,1,iDonor) = GlobalRange.End(1)-2;
          ExpectedExtents(0,2,iDonor) = 0;
          ExpectedExtents(1,0,iDonor) = ExpectedExtents(0,0,iDonor)+1;
          ExpectedExtents(1,1,iDonor) = ExpectedExtents(0,1,iDonor)+1;
          ExpectedExtents(1,2,iDonor) = 1;
          ExpectedCoords(0,iDonor) = 0.;
          ExpectedCoords(1,iDonor) = 0.;
          ExpectedCoords(2,iDonor) = 0.;
          ExpectedInterpCoefs(0,0,iDonor) = 1.;
          ExpectedInterpCoefs(1,0,iDonor) = 1.;
          ExpectedInterpCoefs(2,0,iDonor) = 1.;
          ExpectedDestinations(0,iDonor) = i;
          ExpectedDestinations(1,iDonor) = 0;
          ExpectedDestinations(2,iDonor) = 0;
          ++iDonor;
        }
        EXPECT_THAT(Extents, ElementsAreArray(ExpectedExtents));
        EXPECT_THAT(Coords, ElementsAreArray(ExpectedCoords));
        EXPECT_THAT(InterpCoefs, ElementsAreArray(ExpectedInterpCoefs));
        EXPECT_THAT(Destinations, ElementsAreArray(ExpectedDestinations));

        long long NumReceivers = ConnectivityN.Count();
        EXPECT_EQ(NumReceivers, LocalRange.Size(0));
        const ovk::array<int,2> &Points = ConnectivityN.Points();
        const ovk::array<int,2> &Sources = ConnectivityN.Sources();
        ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumReceivers}});
        ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumReceivers}});
        long long iReceiver = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedPoints(0,iReceiver) = i;
          ExpectedPoints(1,iReceiver) = GlobalRange.End(1)-1;
          ExpectedPoints(2,iReceiver) = 0;
          ExpectedSources(0,iReceiver) = i;
          ExpectedSources(1,iReceiver) = 1;
          ExpectedSources(2,iReceiver) = 0;
          ++iReceiver;
        }
        EXPECT_THAT(Points, ElementsAreArray(Points));
        EXPECT_THAT(Sources, ElementsAreArray(Sources));

      } else {

        EXPECT_EQ(ConnectivityM.Count(), 0);
        EXPECT_EQ(ConnectivityN.Count(), 0);

      }

    }

    if (Grid2IsLocal) {

      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.Begin(1) == GlobalRange.Begin(1);

      const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(2,1);
      const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(1,2);

      if (HasInterface) {

        long long NumDonors = ConnectivityM.Count();
        EXPECT_EQ(NumDonors, LocalRange.Size(0));
        const ovk::array<int,3> &Extents = ConnectivityM.Extents();
        const ovk::array<double,2> &Coords = ConnectivityM.Coords();
        const ovk::array<double,3> &InterpCoefs = ConnectivityM.InterpCoefs();
        const ovk::array<int,2> &Destinations = ConnectivityM.Destinations();
        ovk::array<int,3> ExpectedExtents({{2,ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,2> ExpectedCoords({{ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,3> ExpectedInterpCoefs({{ovk::MAX_DIMS,1,NumDonors}});
        ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumDonors}});
        long long iDonor = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedExtents(0,0,iDonor) = i;
          ExpectedExtents(0,1,iDonor) = 1;
          ExpectedExtents(0,2,iDonor) = 0;
          ExpectedExtents(1,0,iDonor) = ExpectedExtents(0,0,iDonor)+1;
          ExpectedExtents(1,1,iDonor) = ExpectedExtents(0,1,iDonor)+1;
          ExpectedExtents(1,2,iDonor) = 1;
          ExpectedCoords(0,iDonor) = 0.;
          ExpectedCoords(1,iDonor) = 0.;
          ExpectedCoords(2,iDonor) = 0.;
          ExpectedInterpCoefs(0,0,iDonor) = 1.;
          ExpectedInterpCoefs(1,0,iDonor) = 1.;
          ExpectedInterpCoefs(2,0,iDonor) = 1.;
          ExpectedDestinations(0,iDonor) = i;
          ExpectedDestinations(1,iDonor) = Domain.GridInfo(1).GlobalRange().End(1)-1;
          ExpectedDestinations(2,iDonor) = 0;
          ++iDonor;
        }
        EXPECT_THAT(Extents, ElementsAreArray(ExpectedExtents));
        EXPECT_THAT(Coords, ElementsAreArray(ExpectedCoords));
        EXPECT_THAT(InterpCoefs, ElementsAreArray(ExpectedInterpCoefs));
        EXPECT_THAT(Destinations, ElementsAreArray(ExpectedDestinations));

        long long NumReceivers = ConnectivityN.Count();
        EXPECT_EQ(NumReceivers, LocalRange.Size(0));
        const ovk::array<int,2> &Points = ConnectivityN.Points();
        const ovk::array<int,2> &Sources = ConnectivityN.Sources();
        ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumReceivers}});
        ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumReceivers}});
        long long iReceiver = 0;
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedPoints(0,iReceiver) = i;
          ExpectedPoints(1,iReceiver) = 0;
          ExpectedPoints(2,iReceiver) = 0;
          ExpectedSources(0,iReceiver) = i;
          ExpectedSources(1,iReceiver) = Domain.GridInfo(1).GlobalRange().End(1)-2;
          ExpectedSources(2,iReceiver) = 0;
          ++iReceiver;
        }
        EXPECT_THAT(Points, ElementsAreArray(Points));
        EXPECT_THAT(Sources, ElementsAreArray(Sources));

      } else {

        EXPECT_EQ(ConnectivityM.Count(), 0);
        EXPECT_EQ(ConnectivityN.Count(), 0);

      }

    }

  }

}

TEST_F(XINTOUTTests, ImportExtended3D) {

  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < 64);

  if (Comm) {

    ovk::domain Domain = Interface3D(Comm, {{-1.,-1.,-1}, {1.,1.,1.}}, {32,32,32}, {false, false,
      false}, ovk::periodic_storage::UNIQUE);

    Domain.CreateComponent<ovk::connectivity_component>(1);

    ovk::ImportXINTOUT(Domain, 1, "data/XINTOUTTests/XINTOUT_Extended3D.HO.2D",
      "data/XINTOUTTests/XINTOUT_Extended3D.X.2D", 0, MPI_INFO_NULL);

    bool Grid1IsLocal = Domain.GridIsLocal(1);
    bool Grid2IsLocal = Domain.GridIsLocal(2);

    ovk::exchanger Exchanger = ovk::CreateExchanger(Domain.SharedContext());

    Exchanger.Bind(Domain, ovk::exchanger::bindings()
      .SetConnectivityComponentID(1)
    );

    auto &ConnectivityComponent = Domain.Component<ovk::connectivity_component>(1);

    if (Grid1IsLocal) {

      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.End(2) == GlobalRange.End(2);

      const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(1,2);
      const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(2,1);

      if (HasInterface) {

        long long NumDonors = ConnectivityM.Count();
        EXPECT_EQ(NumDonors, LocalRange.Size(0)*LocalRange.Size(1));
        const ovk::array<int,3> &Extents = ConnectivityM.Extents();
        const ovk::array<double,2> &Coords = ConnectivityM.Coords();
        const ovk::array<double,3> &InterpCoefs = ConnectivityM.InterpCoefs();
        const ovk::array<int,2> &Destinations = ConnectivityM.Destinations();
        ovk::array<int,3> ExpectedExtents({{2,ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,2> ExpectedCoords({{ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,3> ExpectedInterpCoefs({{ovk::MAX_DIMS,1,NumDonors}});
        ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumDonors}});
        long long iDonor = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ExpectedExtents(0,0,iDonor) = i;
            ExpectedExtents(0,1,iDonor) = j;
            ExpectedExtents(0,2,iDonor) = GlobalRange.End(2)-2;
            ExpectedExtents(1,0,iDonor) = ExpectedExtents(0,0,iDonor)+1;
            ExpectedExtents(1,1,iDonor) = ExpectedExtents(0,1,iDonor)+1;
            ExpectedExtents(1,2,iDonor) = ExpectedExtents(0,2,iDonor)+1;
            ExpectedCoords(0,iDonor) = 0.;
            ExpectedCoords(1,iDonor) = 0.;
            ExpectedCoords(2,iDonor) = 0.;
            ExpectedInterpCoefs(0,0,iDonor) = 1.;
            ExpectedInterpCoefs(1,0,iDonor) = 1.;
            ExpectedInterpCoefs(2,0,iDonor) = 1.;
            ExpectedDestinations(0,iDonor) = i;
            ExpectedDestinations(1,iDonor) = j;
            ExpectedDestinations(2,iDonor) = 0;
            ++iDonor;
          }
        }
        EXPECT_THAT(Extents, ElementsAreArray(ExpectedExtents));
        EXPECT_THAT(Coords, ElementsAreArray(ExpectedCoords));
        EXPECT_THAT(InterpCoefs, ElementsAreArray(ExpectedInterpCoefs));
        EXPECT_THAT(Destinations, ElementsAreArray(ExpectedDestinations));

        long long NumReceivers = ConnectivityN.Count();
        EXPECT_EQ(NumReceivers, LocalRange.Size(0)*LocalRange.Size(1));
        const ovk::array<int,2> &Points = ConnectivityN.Points();
        const ovk::array<int,2> &Sources = ConnectivityN.Sources();
        ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumReceivers}});
        ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumReceivers}});
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ExpectedPoints(0,iReceiver) = i;
            ExpectedPoints(1,iReceiver) = j;
            ExpectedPoints(2,iReceiver) = GlobalRange.End(2)-1;
            ExpectedSources(0,iReceiver) = i;
            ExpectedSources(1,iReceiver) = j;
            ExpectedSources(2,iReceiver) = 1;
            ++iReceiver;
          }
        }
        EXPECT_THAT(Points, ElementsAreArray(Points));
        EXPECT_THAT(Sources, ElementsAreArray(Sources));

      } else {

        EXPECT_EQ(ConnectivityM.Count(), 0);
        EXPECT_EQ(ConnectivityN.Count(), 0);

      }
    }

    if (Grid2IsLocal) {

      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.Begin(2) == GlobalRange.Begin(2);

      const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(2,1);
      const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(1,2);

      if (HasInterface) {

        long long NumDonors = ConnectivityM.Count();
        EXPECT_EQ(NumDonors, LocalRange.Size(0)*LocalRange.Size(1));
        const ovk::array<int,3> &Extents = ConnectivityM.Extents();
        const ovk::array<double,2> &Coords = ConnectivityM.Coords();
        const ovk::array<double,3> &InterpCoefs = ConnectivityM.InterpCoefs();
        const ovk::array<int,2> &Destinations = ConnectivityM.Destinations();
        ovk::array<int,3> ExpectedExtents({{2,ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,2> ExpectedCoords({{ovk::MAX_DIMS,NumDonors}});
        ovk::array<double,3> ExpectedInterpCoefs({{ovk::MAX_DIMS,1,NumDonors}});
        ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumDonors}});
        long long iDonor = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ExpectedExtents(0,0,iDonor) = i;
            ExpectedExtents(0,1,iDonor) = j;
            ExpectedExtents(0,2,iDonor) = 1;
            ExpectedExtents(1,0,iDonor) = ExpectedExtents(0,0,iDonor)+1;
            ExpectedExtents(1,1,iDonor) = ExpectedExtents(0,1,iDonor)+1;
            ExpectedExtents(1,2,iDonor) = ExpectedExtents(0,2,iDonor)+1;
            ExpectedCoords(0,iDonor) = 0.;
            ExpectedCoords(1,iDonor) = 0.;
            ExpectedCoords(2,iDonor) = 0.;
            ExpectedInterpCoefs(0,0,iDonor) = 1.;
            ExpectedInterpCoefs(1,0,iDonor) = 1.;
            ExpectedInterpCoefs(2,0,iDonor) = 1.;
            ExpectedDestinations(0,iDonor) = i;
            ExpectedDestinations(1,iDonor) = j;
            ExpectedDestinations(2,iDonor) = Domain.GridInfo(1).GlobalRange().End(2)-1;
            ++iDonor;
          }
        }
        EXPECT_THAT(Extents, ElementsAreArray(ExpectedExtents));
        EXPECT_THAT(Coords, ElementsAreArray(ExpectedCoords));
        EXPECT_THAT(InterpCoefs, ElementsAreArray(ExpectedInterpCoefs));
        EXPECT_THAT(Destinations, ElementsAreArray(ExpectedDestinations));

        long long NumReceivers = ConnectivityN.Count();
        EXPECT_EQ(NumReceivers, LocalRange.Size(0)*LocalRange.Size(1));
        const ovk::array<int,2> &Points = ConnectivityN.Points();
        const ovk::array<int,2> &Sources = ConnectivityN.Sources();
        ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumReceivers}});
        ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumReceivers}});
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ExpectedPoints(0,iReceiver) = i;
            ExpectedPoints(1,iReceiver) = j;
            ExpectedPoints(2,iReceiver) = 0;
            ExpectedSources(0,iReceiver) = i;
            ExpectedSources(1,iReceiver) = j;
            ExpectedSources(2,iReceiver) = Domain.GridInfo(1).GlobalRange().End(2)-2;
            ++iReceiver;
          }
        }
        EXPECT_THAT(Points, ElementsAreArray(Points));
        EXPECT_THAT(Sources, ElementsAreArray(Sources));

      } else {

        EXPECT_EQ(ConnectivityM.Count(), 0);
        EXPECT_EQ(ConnectivityN.Count(), 0);

      }

    }

  }

}
