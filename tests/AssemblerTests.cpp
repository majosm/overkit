// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Assembler.hpp>

#include "tests/MPITest.hpp"
#include "tests/fixtures/CylinderInCylinder.hpp"
#include "tests/fixtures/WavyInWavy.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/CommunicationOps.hpp>
#include <ovk/core/ConnectivityComponent.hpp>
#include <ovk/core/ConnectivityM.hpp>
#include <ovk/core/ConnectivityN.hpp>
#include <ovk/core/Domain.hpp>
#include <ovk/core/FieldOps.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <utility>

using testing::ElementsAreArray;
using testing::Matcher;
using testing::DoubleNear;

class AssemblerTests : public tests::mpi_test {};

using tests::CylinderInCylinder;
using tests::WavyInWavy;

namespace {

ovk::array<ovk::array<int,2>> GetCartesianDecompIntervals(const ovk::domain &Domain, int GridID) {

  bool IsLocal = Domain.GridIsLocal(GridID);

  bool IsGridRoot = false;
  if (IsLocal) {
    const ovk::grid &Grid = Domain.Grid(GridID);
    IsGridRoot = Grid.Comm().Rank() == 0;
  }

  ovk::tuple<int> CartDims;
  if (IsGridRoot) {
    const ovk::grid &Grid = Domain.Grid(GridID);
    CartDims = ovk::GetCartCommDims(Grid.Comm());
  }
  ovk::core::BroadcastAnySource(CartDims.Data(), ovk::MAX_DIMS, MPI_INT, IsGridRoot, Domain.Comm());

  ovk::array<ovk::array<int,2>> DecompIntervals({Domain.Dimension()});

  for (int iDim = 0; iDim < Domain.Dimension(); ++iDim) {
    DecompIntervals(iDim).Resize({{CartDims(iDim),2}});
    if (IsLocal) {
      const ovk::grid &Grid = Domain.Grid(GridID);
      ovk::tuple<int> CartCoords = ovk::GetCartCommCoords(Grid.Comm());
      ovk::comm LineComm = ovk::CreateSubsetComm(Grid.Comm(), CartCoords((iDim+1) % ovk::MAX_DIMS)
        == 0 && CartCoords((iDim+2) % ovk::MAX_DIMS) == 0);
      if (LineComm) {
        int Interval[] = {Grid.LocalRange().Begin(iDim), Grid.LocalRange().End(iDim)};
        MPI_Gather(Interval, 2, MPI_INT, DecompIntervals(iDim).Data(), 2, MPI_INT, 0, LineComm);
      }
    }
    // Grid root is also the root of each line comm, so no further communication is necessary to get
    // the data to it
    ovk::core::BroadcastAnySource(DecompIntervals(iDim).Data(), DecompIntervals(iDim).Count(),
      MPI_INT, IsGridRoot, Domain.Comm());
  }

  return DecompIntervals;

};

int FindRank(const ovk::array<ovk::array<int,2>> &DecompIntervals, const ovk::tuple<int> &Point) {

  ovk::tuple<int> CartDims = {1,1,1};
  for (int iDim = 0; iDim < DecompIntervals.Count(); ++iDim) {
    CartDims(iDim) = int(DecompIntervals(iDim).Size(0));
  }

  ovk::tuple<int> CartCoords = {0,0,0};

  for (int iDim = 0; iDim < DecompIntervals.Count(); ++iDim) {
    while (Point(iDim) >= DecompIntervals(iDim)(CartCoords(iDim),1)) {
      ++CartCoords(iDim);
    }
  }

  return (CartCoords(0)*CartDims(1) + CartCoords(1))*CartDims(2) + CartCoords(2);

};

ovk::range CreateCellCoverRange(const ovk::grid &Grid) {

  const ovk::range &CellGlobalRange = Grid.CellGlobalRange();
  const ovk::range &CellLocalRange = Grid.CellLocalRange();
  const ovk::tuple<bool> &Periodic = Grid.Periodic();

  ovk::range CellCoverRange = ovk::MakeEmptyRange(Grid.Dimension());

  for (int iDim = 0; iDim < Grid.Dimension(); ++iDim) {
    if (CellLocalRange.Begin(iDim) > 0 || (Periodic(iDim) && CellLocalRange.End(iDim) !=
      CellGlobalRange.End(iDim))) {
      CellCoverRange.Begin(iDim) = CellLocalRange.Begin(iDim)-1;
    } else {
      CellCoverRange.Begin(iDim) = CellLocalRange.Begin(iDim);
    }
    CellCoverRange.End(iDim) = CellLocalRange.End(iDim);
  }

  return CellCoverRange;

}

}

TEST_F(AssemblerTests, Overlap2D) {

  int NumProc = TestComm().Size();
  // Avoid sizes that make decomposition too small
  int AllowedSubsetSizes[] = {1, 2, 4, 6, 8, 12, 16, 18};
  int SubsetSize = 1;
  for (int Size : AllowedSubsetSizes) {
    if (Size > NumProc) break;
    SubsetSize = Size;
  }
  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < SubsetSize);

  if (Comm) {

    ovk::domain Domain = WavyInWavy(2, Comm, 10, true);

    Domain.CreateComponent<ovk::overlap_component>(3);
    Domain.CreateComponent<ovk::connectivity_component>(4);

    bool BackgroundIsLocal = Domain.GridIsLocal(1);
    bool ForegroundIsLocal = Domain.GridIsLocal(2);

    ovk::assembler Assembler = ovk::CreateAssembler(Domain.SharedContext());

    Assembler.Bind(Domain, ovk::assembler::bindings()
      .SetGeometryComponentID(1)
      .SetStateComponentID(2)
      .SetOverlapComponentID(3)
      .SetConnectivityComponentID(4)
    );

    {
      auto OptionsEditHandle = Assembler.EditOptions();
      ovk::assembler::options &Options = *OptionsEditHandle;
      Options.SetOverlappable({2,1}, true);
      Options.SetOverlappable({1,2}, true);
    }

    Assembler.Assemble();

    auto BackgroundDecompIntervals = GetCartesianDecompIntervals(Domain, 1);
    auto ForegroundDecompIntervals = GetCartesianDecompIntervals(Domain, 2);

    auto FindBackgroundRank = [&](const ovk::tuple<int> &Point) -> int {
      return FindRank(BackgroundDecompIntervals, Point);
    };

    auto FindForegroundRank = [&](const ovk::tuple<int> &Point) -> int {
      return Comm.Size()/2 + FindRank(ForegroundDecompIntervals, Point);
    };

    // Sanity check
    int LowerCornerRank = FindBackgroundRank({0,0,0});
    int UpperCornerRank = FindBackgroundRank(Domain.GridInfo(1).GlobalRange().Size() -
      ovk::tuple<int>(1,1,0));
    EXPECT_EQ(LowerCornerRank, 0);
    EXPECT_EQ(UpperCornerRank, ovk::Max(Comm.Size()/2,1)-1);
    LowerCornerRank = FindForegroundRank({0,0,0});
    UpperCornerRank = FindForegroundRank(Domain.GridInfo(2).GlobalRange().Size() -
      ovk::tuple<int>(1,1,0));
    EXPECT_EQ(LowerCornerRank, Comm.Size()/2);
    EXPECT_EQ(UpperCornerRank, Comm.Size()-1);

    ovk::range HoleRange = {{4,4,0}, {6,6,1}};
    ovk::range CellHoleRange = {{3,3,0}, {6,6,1}};

    auto &OverlapComponent = Domain.Component<ovk::overlap_component>(3);

    if (BackgroundIsLocal) {

      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      const ovk::range &CellLocalRange = Grid.CellLocalRange();

      const ovk::overlap_m &OverlapM = OverlapComponent.OverlapM({1,2});

      ovk::range ForegroundRange = ovk::MakeEmptyRange(2);
      for (int iDim = 0; iDim < 2; ++iDim) {
        ForegroundRange.Begin(iDim) = 2*(CellLocalRange.Begin(iDim)-2)-2;
        ForegroundRange.End(iDim) = 2*(CellLocalRange.End(iDim)-2);
      }
      ForegroundRange = ovk::IntersectRanges(ForegroundRange, Domain.GridInfo(2).GlobalRange());

      long long NumOverlapping = 0;
      for (int j = ForegroundRange.Begin(1); j < ForegroundRange.End(1); ++j) {
        for (int i = ForegroundRange.Begin(0); i < ForegroundRange.End(0); ++i) {
          ovk::tuple<int> Cell = {2+i/2,2+j/2,0};
          if (!CellHoleRange.Contains(Cell)) {
            ++NumOverlapping;
          }
        }
      }

      EXPECT_EQ(OverlapM.Size(), NumOverlapping);

      ovk::array<int,2> ExpectedCells({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<Matcher<double>,2> ExpectedCoords({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedDestinationRanks({NumOverlapping}, -1);
      long long iOverlapping = 0;
      for (int j = ForegroundRange.Begin(1); j < ForegroundRange.End(1); ++j) {
        for (int i = ForegroundRange.Begin(0); i < ForegroundRange.End(0); ++i) {
          ovk::tuple<int> Cell = {2+i/2,2+j/2,0};
          if (!CellHoleRange.Contains(Cell)) {
            ExpectedCells(0,iOverlapping) = Cell(0);
            ExpectedCells(1,iOverlapping) = Cell(1);
            ExpectedCells(2,iOverlapping) = 0;
            ExpectedCoords(0,iOverlapping) = DoubleNear(i % 2 == 0 ? 0.25 : 0.75, 1.e-2);
            ExpectedCoords(1,iOverlapping) = DoubleNear(j % 2 == 0 ? 0.25 : 0.75, 1.e-2);
            ExpectedCoords(2,iOverlapping) = 0.;
            ExpectedDestinations(0,iOverlapping) = i;
            ExpectedDestinations(1,iOverlapping) = j;
            ExpectedDestinations(2,iOverlapping) = 0;
            if (CellLocalRange.Contains(Cell)) {
              ExpectedDestinationRanks(iOverlapping) = FindForegroundRank({i,j,0});
            }
            ++iOverlapping;
          }
        }
      }

      EXPECT_THAT(OverlapM.Cells(), ElementsAreArray(ExpectedCells));
      EXPECT_THAT(OverlapM.Coords(), ElementsAreArray(ExpectedCoords));
      EXPECT_THAT(OverlapM.Destinations(), ElementsAreArray(ExpectedDestinations));
      EXPECT_THAT(OverlapM.DestinationRanks(), ElementsAreArray(ExpectedDestinationRanks));

      const ovk::overlap_n &OverlapN = OverlapComponent.OverlapN({2,1});

      ovk::field<bool> ExpectedMask(Grid.LocalRange());
      NumOverlapping = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedMask(i,j,0) = i >= 3 && i < 7 && j >= 3 && j < 7 && !HoleRange.Contains({i,j,0});
          if (ExpectedMask(i,j,0)) {
            ++NumOverlapping;
          }
        }
      }

      EXPECT_EQ(OverlapN.Size(), NumOverlapping);
      EXPECT_THAT(OverlapN.Mask(), ElementsAreArray(ExpectedMask));

      ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedSourceRanks({NumOverlapping});
      iOverlapping = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          if (ExpectedMask(i,j,0)) {
            ExpectedPoints(0,iOverlapping) = i;
            ExpectedPoints(1,iOverlapping) = j;
            ExpectedPoints(2,iOverlapping) = 0;
            ovk::tuple<int> Source = {2*(i-3)+1, 2*(j-3)+1, 0};
            ExpectedSources(0,iOverlapping) = Source(0);
            ExpectedSources(1,iOverlapping) = Source(1);
            ExpectedSources(2,iOverlapping) = 0;
            ExpectedSourceRanks(iOverlapping) = FindForegroundRank(Source);
            ++iOverlapping;
          }
        }
      }

      EXPECT_THAT(OverlapN.Points(), ElementsAreArray(ExpectedPoints));
      EXPECT_THAT(OverlapN.Sources(), ElementsAreArray(ExpectedSources));
      EXPECT_THAT(OverlapN.SourceRanks(), ElementsAreArray(ExpectedSourceRanks));

    }

    if (ForegroundIsLocal) {

      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      const ovk::range &CellLocalRange = Grid.CellLocalRange();

      const ovk::overlap_m &OverlapM = OverlapComponent.OverlapM({2,1});

      ovk::range BackgroundRange = ovk::MakeEmptyRange(2);
      for (int iDim = 0; iDim < 2; ++iDim) {
        BackgroundRange.Begin(iDim) = ovk::Max(3+(CellLocalRange.Begin(iDim)-1)/2, 3);
        BackgroundRange.End(iDim) = ovk::Min(3+(CellLocalRange.End(iDim))/2, 7);
      }

      long long NumOverlapping = 0;
      for (int j = BackgroundRange.Begin(1); j < BackgroundRange.End(1); ++j) {
        for (int i = BackgroundRange.Begin(0); i < BackgroundRange.End(0); ++i) {
          if (!HoleRange.Contains({i,j,0})) {
            ++NumOverlapping;
          }
        }
      }

      EXPECT_EQ(OverlapM.Size(), NumOverlapping);

      ovk::array<int,2> ExpectedCells({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<Matcher<double>,2> ExpectedCoords({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedDestinationRanks({NumOverlapping}, -1);
      long long iOverlapping = 0;
      for (int j = BackgroundRange.Begin(1); j < BackgroundRange.End(1); ++j) {
        for (int i = BackgroundRange.Begin(0); i < BackgroundRange.End(0); ++i) {
          if (!HoleRange.Contains({i,j,0})) {
            ovk::tuple<int> Cell = {2*(i-3)+1,2*(j-3)+1,0};
            ExpectedCells(0,iOverlapping) = Cell(0);
            ExpectedCells(1,iOverlapping) = Cell(1);
            ExpectedCells(2,iOverlapping) = 0;
            ExpectedCoords(0,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(1,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(2,iOverlapping) = 0.;
            ExpectedDestinations(0,iOverlapping) = i;
            ExpectedDestinations(1,iOverlapping) = j;
            ExpectedDestinations(2,iOverlapping) = 0;
            if (CellLocalRange.Contains(Cell)) {
              ExpectedDestinationRanks(iOverlapping) = FindBackgroundRank({i,j,0});
            }
            ++iOverlapping;
          }
        }
      }

      EXPECT_THAT(OverlapM.Cells(), ElementsAreArray(ExpectedCells));
      EXPECT_THAT(OverlapM.Coords(), ElementsAreArray(ExpectedCoords));
      EXPECT_THAT(OverlapM.Destinations(), ElementsAreArray(ExpectedDestinations));
      EXPECT_THAT(OverlapM.DestinationRanks(), ElementsAreArray(ExpectedDestinationRanks));

      const ovk::overlap_n &OverlapN = OverlapComponent.OverlapN({1,2});

      ovk::field<bool> ExpectedMask(Grid.LocalRange());
      NumOverlapping = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ovk::tuple<int> Source = {2+i/2, 2+j/2, 0};
          ExpectedMask(i,j,0) = !CellHoleRange.Contains(Source);
          if (ExpectedMask(i,j,0)) {
            ++NumOverlapping;
          }
        }
      }

      EXPECT_EQ(OverlapN.Size(), NumOverlapping);
      EXPECT_THAT(OverlapN.Mask(), ElementsAreArray(ExpectedMask));

      ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedSourceRanks({NumOverlapping});
      iOverlapping = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          if (ExpectedMask(i,j,0)) {
            ExpectedPoints(0,iOverlapping) = i;
            ExpectedPoints(1,iOverlapping) = j;
            ExpectedPoints(2,iOverlapping) = 0;
            ovk::tuple<int> Source = {2+i/2, 2+j/2, 0};
            ExpectedSources(0,iOverlapping) = Source(0);
            ExpectedSources(1,iOverlapping) = Source(1);
            ExpectedSources(2,iOverlapping) = 0;
            ExpectedSourceRanks(iOverlapping) = FindBackgroundRank(Source);
            ++iOverlapping;
          }
        }
      }

      EXPECT_THAT(OverlapN.Points(), ElementsAreArray(ExpectedPoints));
      EXPECT_THAT(OverlapN.Sources(), ElementsAreArray(ExpectedSources));
      EXPECT_THAT(OverlapN.SourceRanks(), ElementsAreArray(ExpectedSourceRanks));

    }

  }

}

TEST_F(AssemblerTests, Overlap3D) {

  int NumProc = TestComm().Size();
  // Avoid sizes that make decomposition too small
  int AllowedSubsetSizes[] = {1, 2, 4, 6, 8, 12, 16, 18, 32, 36, 48, 54};
  int SubsetSize = 1;
  for (int Size : AllowedSubsetSizes) {
    if (Size > NumProc) break;
    SubsetSize = Size;
  }
  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < SubsetSize);

  if (Comm) {

    ovk::domain Domain = WavyInWavy(3, Comm, 10, true);

    Domain.CreateComponent<ovk::overlap_component>(3);
    Domain.CreateComponent<ovk::connectivity_component>(4);

    bool BackgroundIsLocal = Domain.GridIsLocal(1);
    bool ForegroundIsLocal = Domain.GridIsLocal(2);

    ovk::assembler Assembler = ovk::CreateAssembler(Domain.SharedContext());

    Assembler.Bind(Domain, ovk::assembler::bindings()
      .SetGeometryComponentID(1)
      .SetStateComponentID(2)
      .SetOverlapComponentID(3)
      .SetConnectivityComponentID(4)
    );

    {
      auto OptionsEditHandle = Assembler.EditOptions();
      ovk::assembler::options &Options = *OptionsEditHandle;
      Options.SetOverlappable({2,1}, true);
      Options.SetOverlappable({1,2}, true);
    }

    Assembler.Assemble();

    auto BackgroundDecompIntervals = GetCartesianDecompIntervals(Domain, 1);
    auto ForegroundDecompIntervals = GetCartesianDecompIntervals(Domain, 2);

    auto FindBackgroundRank = [&](const ovk::tuple<int> &Point) -> int {
      return FindRank(BackgroundDecompIntervals, Point);
    };

    auto FindForegroundRank = [&](const ovk::tuple<int> &Point) -> int {
      return Comm.Size()/2 + FindRank(ForegroundDecompIntervals, Point);
    };

    // Sanity check
    int LowerCornerRank = FindBackgroundRank({0,0,0});
    int UpperCornerRank = FindBackgroundRank(Domain.GridInfo(1).GlobalRange().Size() -
      ovk::tuple<int>(1,1,1));
    EXPECT_EQ(LowerCornerRank, 0);
    EXPECT_EQ(UpperCornerRank, ovk::Max(Comm.Size()/2,1)-1);
    LowerCornerRank = FindForegroundRank({0,0,0});
    UpperCornerRank = FindForegroundRank(Domain.GridInfo(2).GlobalRange().Size() -
      ovk::tuple<int>(1,1,1));
    EXPECT_EQ(LowerCornerRank, Comm.Size()/2);
    EXPECT_EQ(UpperCornerRank, Comm.Size()-1);

    ovk::range HoleRange = {{4,4,4}, {6,6,6}};
    ovk::range CellHoleRange = {{3,3,3}, {6,6,6}};

    auto &OverlapComponent = Domain.Component<ovk::overlap_component>(3);

    if (BackgroundIsLocal) {

      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      const ovk::range &CellLocalRange = Grid.CellLocalRange();

      const ovk::overlap_m &OverlapM = OverlapComponent.OverlapM({1,2});

      ovk::range ForegroundRange;
      for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
        ForegroundRange.Begin(iDim) = 2*(CellLocalRange.Begin(iDim)-2)-2;
        ForegroundRange.End(iDim) = 2*(CellLocalRange.End(iDim)-2);
      }
      ForegroundRange = ovk::IntersectRanges(ForegroundRange, Domain.GridInfo(2).GlobalRange());

      long long NumOverlapping = 0;
      for (int k = ForegroundRange.Begin(2); k < ForegroundRange.End(2); ++k) {
        for (int j = ForegroundRange.Begin(1); j < ForegroundRange.End(1); ++j) {
          for (int i = ForegroundRange.Begin(0); i < ForegroundRange.End(0); ++i) {
            ovk::tuple<int> Cell = {2+i/2,2+j/2,2+k/2};
            if (!CellHoleRange.Contains(Cell)) {
              ++NumOverlapping;
            }
          }
        }
      }

      EXPECT_EQ(OverlapM.Size(), NumOverlapping);

      ovk::array<int,2> ExpectedCells({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<Matcher<double>,2> ExpectedCoords({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedDestinationRanks({NumOverlapping}, -1);
      long long iOverlapping = 0;
      for (int k = ForegroundRange.Begin(2); k < ForegroundRange.End(2); ++k) {
        for (int j = ForegroundRange.Begin(1); j < ForegroundRange.End(1); ++j) {
          for (int i = ForegroundRange.Begin(0); i < ForegroundRange.End(0); ++i) {
            ovk::tuple<int> Cell = {2+i/2,2+j/2,2+k/2};
            if (!CellHoleRange.Contains(Cell)) {
              ExpectedCells(0,iOverlapping) = Cell(0);
              ExpectedCells(1,iOverlapping) = Cell(1);
              ExpectedCells(2,iOverlapping) = Cell(2);
              ExpectedCoords(0,iOverlapping) = DoubleNear(i % 2 == 0 ? 0.25 : 0.75, 1.e-2);
              ExpectedCoords(1,iOverlapping) = DoubleNear(j % 2 == 0 ? 0.25 : 0.75, 1.e-2);
              ExpectedCoords(2,iOverlapping) = DoubleNear(k % 2 == 0 ? 0.25 : 0.75, 1.e-2);
              ExpectedDestinations(0,iOverlapping) = i;
              ExpectedDestinations(1,iOverlapping) = j;
              ExpectedDestinations(2,iOverlapping) = k;
              if (CellLocalRange.Contains(Cell)) {
                ExpectedDestinationRanks(iOverlapping) = FindForegroundRank({i,j,k});
              }
              ++iOverlapping;
            }
          }
        }
      }

      EXPECT_THAT(OverlapM.Cells(), ElementsAreArray(ExpectedCells));
      EXPECT_THAT(OverlapM.Coords(), ElementsAreArray(ExpectedCoords));
      EXPECT_THAT(OverlapM.Destinations(), ElementsAreArray(ExpectedDestinations));
      EXPECT_THAT(OverlapM.DestinationRanks(), ElementsAreArray(ExpectedDestinationRanks));

      const ovk::overlap_n &OverlapN = OverlapComponent.OverlapN({2,1});

      ovk::field<bool> ExpectedMask(Grid.LocalRange());
      NumOverlapping = 0;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ExpectedMask(i,j,k) = i >= 3 && i < 7 && j >= 3 && j < 7 && k >= 3 && k < 7 &&
              !HoleRange.Contains({i,j,k});
            if (ExpectedMask(i,j,k)) {
              ++NumOverlapping;
            }
          }
        }
      }

      EXPECT_EQ(OverlapN.Size(), NumOverlapping);
      EXPECT_THAT(OverlapN.Mask(), ElementsAreArray(ExpectedMask));

      ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedSourceRanks({NumOverlapping});
      iOverlapping = 0;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            if (ExpectedMask(i,j,k)) {
              ExpectedPoints(0,iOverlapping) = i;
              ExpectedPoints(1,iOverlapping) = j;
              ExpectedPoints(2,iOverlapping) = k;
              ovk::tuple<int> Source = {2*(i-3)+1, 2*(j-3)+1, 2*(k-3)+1};
              ExpectedSources(0,iOverlapping) = Source(0);
              ExpectedSources(1,iOverlapping) = Source(1);
              ExpectedSources(2,iOverlapping) = Source(2);
              ExpectedSourceRanks(iOverlapping) = FindForegroundRank(Source);
              ++iOverlapping;
            }
          }
        }
      }

      EXPECT_THAT(OverlapN.Points(), ElementsAreArray(ExpectedPoints));
      EXPECT_THAT(OverlapN.Sources(), ElementsAreArray(ExpectedSources));
      EXPECT_THAT(OverlapN.SourceRanks(), ElementsAreArray(ExpectedSourceRanks));

    }

    if (ForegroundIsLocal) {

      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      const ovk::range &CellLocalRange = Grid.CellLocalRange();

      const ovk::overlap_m &OverlapM = OverlapComponent.OverlapM({2,1});

      ovk::range BackgroundRange;
      for (int iDim = 0; iDim < ovk::MAX_DIMS; ++iDim) {
        BackgroundRange.Begin(iDim) = ovk::Max(3+(CellLocalRange.Begin(iDim)-1)/2, 3);
        BackgroundRange.End(iDim) = ovk::Min(3+(CellLocalRange.End(iDim))/2, 7);
      }

      long long NumOverlapping = 0;
      for (int k = BackgroundRange.Begin(2); k < BackgroundRange.End(2); ++k) {
        for (int j = BackgroundRange.Begin(1); j < BackgroundRange.End(1); ++j) {
          for (int i = BackgroundRange.Begin(0); i < BackgroundRange.End(0); ++i) {
            if (!HoleRange.Contains({i,j,k})) {
              ++NumOverlapping;
            }
          }
        }
      }

      EXPECT_EQ(OverlapM.Size(), NumOverlapping);

      ovk::array<int,2> ExpectedCells({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<Matcher<double>,2> ExpectedCoords({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedDestinationRanks({NumOverlapping}, -1);
      long long iOverlapping = 0;
      for (int k = BackgroundRange.Begin(2); k < BackgroundRange.End(2); ++k) {
        for (int j = BackgroundRange.Begin(1); j < BackgroundRange.End(1); ++j) {
          for (int i = BackgroundRange.Begin(0); i < BackgroundRange.End(0); ++i) {
            if (!HoleRange.Contains({i,j,k})) {
              ovk::tuple<int> Cell = {2*(i-3)+1,2*(j-3)+1,2*(k-3)+1};
              ExpectedCells(0,iOverlapping) = Cell(0);
              ExpectedCells(1,iOverlapping) = Cell(1);
              ExpectedCells(2,iOverlapping) = Cell(2);
              ExpectedCoords(0,iOverlapping) = DoubleNear(0.5, 1.e-2);
              ExpectedCoords(1,iOverlapping) = DoubleNear(0.5, 1.e-2);
              ExpectedCoords(2,iOverlapping) = DoubleNear(0.5, 1.e-2);
              ExpectedDestinations(0,iOverlapping) = i;
              ExpectedDestinations(1,iOverlapping) = j;
              ExpectedDestinations(2,iOverlapping) = k;
              if (CellLocalRange.Contains(Cell)) {
                ExpectedDestinationRanks(iOverlapping) = FindBackgroundRank({i,j,k});
              }
              ++iOverlapping;
            }
          }
        }
      }

      EXPECT_THAT(OverlapM.Cells(), ElementsAreArray(ExpectedCells));
      EXPECT_THAT(OverlapM.Coords(), ElementsAreArray(ExpectedCoords));
      EXPECT_THAT(OverlapM.Destinations(), ElementsAreArray(ExpectedDestinations));
      EXPECT_THAT(OverlapM.DestinationRanks(), ElementsAreArray(ExpectedDestinationRanks));

      const ovk::overlap_n &OverlapN = OverlapComponent.OverlapN({1,2});

      ovk::field<bool> ExpectedMask(Grid.LocalRange());
      NumOverlapping = 0;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ovk::tuple<int> Source = {2+i/2, 2+j/2, 2+k/2};
            ExpectedMask(i,j,k) = !CellHoleRange.Contains(Source);
            if (ExpectedMask(i,j,k)) {
              ++NumOverlapping;
            }
          }
        }
      }

      EXPECT_EQ(OverlapN.Size(), NumOverlapping);
      EXPECT_THAT(OverlapN.Mask(), ElementsAreArray(ExpectedMask));

      ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedSourceRanks({NumOverlapping});
      iOverlapping = 0;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            if (ExpectedMask(i,j,k)) {
              ExpectedPoints(0,iOverlapping) = i;
              ExpectedPoints(1,iOverlapping) = j;
              ExpectedPoints(2,iOverlapping) = k;
              ovk::tuple<int> Source = {2+i/2, 2+j/2, 2+k/2};
              ExpectedSources(0,iOverlapping) = Source(0);
              ExpectedSources(1,iOverlapping) = Source(1);
              ExpectedSources(2,iOverlapping) = Source(2);
              ExpectedSourceRanks(iOverlapping) = FindBackgroundRank(Source);
              ++iOverlapping;
            }
          }
        }
      }

      EXPECT_THAT(OverlapN.Points(), ElementsAreArray(ExpectedPoints));
      EXPECT_THAT(OverlapN.Sources(), ElementsAreArray(ExpectedSources));
      EXPECT_THAT(OverlapN.SourceRanks(), ElementsAreArray(ExpectedSourceRanks));

    }

  }

}

TEST_F(AssemblerTests, OverlapPeriodic) {

  int NumProc = TestComm().Size();
  int SubsetSize = ovk::Min(16, NumProc % 2 == 0 ? NumProc : ovk::Max(NumProc-1,1));
  ovk::comm Comm = CreateSubsetComm(TestComm(), TestComm().Rank() < SubsetSize);

  // Unique
  if (Comm) {

    ovk::domain Domain = CylinderInCylinder(Comm, 8, 4, true, {false,true,false},
      ovk::periodic_storage::UNIQUE);

    Domain.CreateComponent<ovk::overlap_component>(3);
    Domain.CreateComponent<ovk::connectivity_component>(4);

    bool InnerIsLocal = Domain.GridIsLocal(1);
    bool OuterIsLocal = Domain.GridIsLocal(2);

    ovk::assembler Assembler = ovk::CreateAssembler(Domain.SharedContext());

    Assembler.Bind(Domain, ovk::assembler::bindings()
      .SetGeometryComponentID(1)
      .SetStateComponentID(2)
      .SetOverlapComponentID(3)
      .SetConnectivityComponentID(4)
    );

    {
      auto OptionsEditHandle = Assembler.EditOptions();
      ovk::assembler::options &Options = *OptionsEditHandle;
      Options.SetOverlappable({2,1}, true);
      Options.SetOverlappable({1,2}, true);
    }

    Assembler.Assemble();

    auto &OverlapComponent = Domain.Component<ovk::overlap_component>(3);

    if (InnerIsLocal) {

      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      const ovk::range &CellGlobalRange = Grid.CellGlobalRange();
      const ovk::range &CellLocalRange = Grid.CellLocalRange();

      const ovk::overlap_m &OverlapM = OverlapComponent.OverlapM({1,2});

      ovk::range CellCoverRange = CreateCellCoverRange(Grid);
      ovk::range OverlapMRange = ovk::UnionRanges(CellGlobalRange, CellCoverRange);
      OverlapMRange.Begin(0) = 2;
      OverlapMRange = ovk::IntersectRanges(OverlapMRange, CellCoverRange);

      long long NumOverlapping = OverlapMRange.Count();

      EXPECT_EQ(OverlapM.Size(), NumOverlapping);

      ovk::array<int,2> ExpectedCells({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<Matcher<double>,2> ExpectedCoords({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedDestinationRanks({NumOverlapping}, -1);
      long long iOverlapping = 0;
      for (int j = OverlapMRange.Begin(1); j < OverlapMRange.End(1); ++j) {
        for (int i = OverlapMRange.Begin(0); i < OverlapMRange.End(0); ++i) {
          if (j >= 0) {
            ExpectedCells(0,iOverlapping) = i;
            ExpectedCells(1,iOverlapping) = j;
            ExpectedCells(2,iOverlapping) = 0;
            ExpectedCoords(0,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(1,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(2,iOverlapping) = 0.;
            ExpectedDestinations(0,iOverlapping) = i-2;
            ExpectedDestinations(1,iOverlapping) = j;
            ExpectedDestinations(2,iOverlapping) = 0;
            if (CellLocalRange.Contains({i,j,0})) {
              ExpectedDestinationRanks(iOverlapping) = Comm.Rank() + Comm.Size()/2;
            }
            ++iOverlapping;
          }
        }
      }
      // Must do these after the other cells because overlap M data is ordered by destination
      // point index
      for (int j = OverlapMRange.Begin(1); j < OverlapMRange.End(1); ++j) {
        for (int i = OverlapMRange.Begin(0); i < OverlapMRange.End(0); ++i) {
          if (j < 0) {
            ExpectedCells(0,iOverlapping) = i;
            ExpectedCells(1,iOverlapping) = j;
            ExpectedCells(2,iOverlapping) = 0;
            ExpectedCoords(0,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(1,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(2,iOverlapping) = 0.;
            ExpectedDestinations(0,iOverlapping) = i-2;
            ExpectedDestinations(1,iOverlapping) = 15;
            ExpectedDestinations(2,iOverlapping) = 0;
            ++iOverlapping;
          }
        }
      }

      EXPECT_THAT(OverlapM.Cells(), ElementsAreArray(ExpectedCells));
      EXPECT_THAT(OverlapM.Coords(), ElementsAreArray(ExpectedCoords));
      EXPECT_THAT(OverlapM.Destinations(), ElementsAreArray(ExpectedDestinations));
      EXPECT_THAT(OverlapM.DestinationRanks(), ElementsAreArray(ExpectedDestinationRanks));

      const ovk::overlap_n &OverlapN = OverlapComponent.OverlapN({2,1});

      ovk::field<bool> ExpectedMask(Grid.LocalRange());
      NumOverlapping = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedMask(i,j,0) = i > 2;
          if (ExpectedMask(i,j,0)) {
            ++NumOverlapping;
          }
        }
      }

      EXPECT_EQ(OverlapN.Size(), NumOverlapping);
      EXPECT_THAT(OverlapN.Mask(), ElementsAreArray(ExpectedMask));

      ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedSourceRanks({NumOverlapping});
      iOverlapping = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          if (ExpectedMask(i,j,0)) {
            ExpectedPoints(0,iOverlapping) = i;
            ExpectedPoints(1,iOverlapping) = j;
            ExpectedPoints(2,iOverlapping) = 0;
            if (j > 0) {
              ExpectedSources(0,iOverlapping) = i-3;
              ExpectedSources(1,iOverlapping) = j-1;
              ExpectedSources(2,iOverlapping) = 0;
              ExpectedSourceRanks(iOverlapping) = Comm.Rank() + Comm.Size()/2 +
                (j == LocalRange.Begin(1) ? -1 : 0);
            } else {
              ExpectedSources(0,iOverlapping) = i-3;
              ExpectedSources(1,iOverlapping) = 15;
              ExpectedSources(2,iOverlapping) = 0;
              ExpectedSourceRanks(iOverlapping) = Comm.Size()-1;
            }
            ++iOverlapping;
          }
        }
      }

      EXPECT_THAT(OverlapN.Points(), ElementsAreArray(ExpectedPoints));
      EXPECT_THAT(OverlapN.Sources(), ElementsAreArray(ExpectedSources));
      EXPECT_THAT(OverlapN.SourceRanks(), ElementsAreArray(ExpectedSourceRanks));

    }

    if (OuterIsLocal) {

      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      const ovk::range &CellGlobalRange = Grid.CellGlobalRange();
      const ovk::range &CellLocalRange = Grid.CellLocalRange();

      const ovk::overlap_m &OverlapM = OverlapComponent.OverlapM({2,1});

      ovk::range CellCoverRange = CreateCellCoverRange(Grid);
      ovk::range OverlapMRange = ovk::UnionRanges(CellGlobalRange, CellCoverRange);
      OverlapMRange.End(0) = 4;
      OverlapMRange = ovk::IntersectRanges(OverlapMRange, CellCoverRange);

      long long NumOverlapping = OverlapMRange.Count();

      EXPECT_EQ(OverlapM.Size(), NumOverlapping);

      ovk::array<int,2> ExpectedCells({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<Matcher<double>,2> ExpectedCoords({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedDestinations({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedDestinationRanks({NumOverlapping}, -1);
      long long iOverlapping = 0;
      // Must do these before the other cells because overlap M data is ordered by destination
      // point index
      for (int j = OverlapMRange.Begin(1); j < OverlapMRange.End(1); ++j) {
        for (int i = OverlapMRange.Begin(0); i < OverlapMRange.End(0); ++i) {
          if (j < 0 || j == CellGlobalRange.End(1)-1) {
            ExpectedCells(0,iOverlapping) = i;
            ExpectedCells(1,iOverlapping) = j;
            ExpectedCells(2,iOverlapping) = 0;
            ExpectedCoords(0,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(1,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(2,iOverlapping) = 0.;
            ExpectedDestinations(0,iOverlapping) = i+3;
            ExpectedDestinations(1,iOverlapping) = 0;
            ExpectedDestinations(2,iOverlapping) = 0;
            if (CellLocalRange.Contains({i,j,0})) {
              ExpectedDestinationRanks(iOverlapping) = 0;
            }
            ++iOverlapping;
          }
        }
      }
      for (int j = OverlapMRange.Begin(1); j < OverlapMRange.End(1); ++j) {
        for (int i = OverlapMRange.Begin(0); i < OverlapMRange.End(0); ++i) {
          if (j >= 0 && j != CellGlobalRange.End(1)-1) {
            ExpectedCells(0,iOverlapping) = i;
            ExpectedCells(1,iOverlapping) = j;
            ExpectedCells(2,iOverlapping) = 0;
            ExpectedCoords(0,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(1,iOverlapping) = DoubleNear(0.5, 1.e-2);
            ExpectedCoords(2,iOverlapping) = 0.;
            ExpectedDestinations(0,iOverlapping) = i+3;
            ExpectedDestinations(1,iOverlapping) = j+1;
            ExpectedDestinations(2,iOverlapping) = 0;
            if (CellLocalRange.Contains({i,j,0})) {
              ExpectedDestinationRanks(iOverlapping) = Comm.Rank() - Comm.Size()/2 +
                (j == CellLocalRange.End(1)-1 ? 1 : 0);
            }
            ++iOverlapping;
          }
        }
      }

      EXPECT_THAT(OverlapM.Cells(), ElementsAreArray(ExpectedCells));
      EXPECT_THAT(OverlapM.Coords(), ElementsAreArray(ExpectedCoords));
      EXPECT_THAT(OverlapM.Destinations(), ElementsAreArray(ExpectedDestinations));
      EXPECT_THAT(OverlapM.DestinationRanks(), ElementsAreArray(ExpectedDestinationRanks));

      const ovk::overlap_n &OverlapN = OverlapComponent.OverlapN({1,2});

      ovk::field<bool> ExpectedMask(Grid.LocalRange());
      NumOverlapping = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ExpectedMask(i,j,0) = i < 4;
          if (ExpectedMask(i,j,0)) {
            ++NumOverlapping;
          }
        }
      }

      EXPECT_EQ(OverlapN.Size(), NumOverlapping);
      EXPECT_THAT(OverlapN.Mask(), ElementsAreArray(ExpectedMask));

      ovk::array<int,2> ExpectedPoints({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int,2> ExpectedSources({{ovk::MAX_DIMS,NumOverlapping}});
      ovk::array<int> ExpectedSourceRanks({NumOverlapping});
      iOverlapping = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          if (ExpectedMask(i,j,0)) {
            ExpectedPoints(0,iOverlapping) = i;
            ExpectedPoints(1,iOverlapping) = j;
            ExpectedPoints(2,iOverlapping) = 0;
            ExpectedSources(0,iOverlapping) = i+2;
            ExpectedSources(1,iOverlapping) = j;
            ExpectedSources(2,iOverlapping) = 0;
            ExpectedSourceRanks(iOverlapping) = Comm.Rank() - Comm.Size()/2;
            ++iOverlapping;
          }
        }
      }

      EXPECT_THAT(OverlapN.Points(), ElementsAreArray(ExpectedPoints));
      EXPECT_THAT(OverlapN.Sources(), ElementsAreArray(ExpectedSources));
      EXPECT_THAT(OverlapN.SourceRanks(), ElementsAreArray(ExpectedSourceRanks));

    }

  }

  // Duplicated
  // Not supported yet
//   if (Comm) {
//   }

}

// TEST_F(AssemblerTests, BoundaryHoleCutting2D) {

//   // Cylinder in box case

// }

// TEST_F(AssemblerTests, BoundaryHoleCutting3D) {

//   // Sphere in box case (yin-yang spherical grids)

// }

// TEST_F(AssemblerTests, BoundaryHoleCuttingBlanked) {

//   // Box in box with boundary in the center + blanking inside

// }
