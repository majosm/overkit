// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/GeometryOps.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <type_traits>

using testing::ElementsAre;
using testing::ElementsAreArray;

class GeometryOpsTests : public tests::mpi_test {};

TEST_F(GeometryOpsTests, CartesianGridCell1D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::CartesianGridCell1D;

  // Above origin, cell interior
  {
    long long Cell = CartesianGridCell1D(1., 2., 6.);
    EXPECT_EQ(Cell, 2);
  }

  // Above origin, lower cell boundary
  {
    long long Cell = CartesianGridCell1D(1., 2., 5.);
    EXPECT_EQ(Cell, 2);
  }

  // Above origin, upper cell boundary
  {
    long long Cell = CartesianGridCell1D(1., 2., 7.);
    EXPECT_EQ(Cell, 3);
  }

  // Below origin, cell interior
  {
    long long Cell = CartesianGridCell1D(1., 2., -6.);
    EXPECT_EQ(Cell, -4);
  }

  // Below origin, lower cell boundary
  {
    long long Cell = CartesianGridCell1D(1., 2., -7.);
    EXPECT_EQ(Cell, -4);
  }

  // Below origin, upper cell boundary
  {
    long long Cell = CartesianGridCell1D(1., 2., -5.);
    EXPECT_EQ(Cell, -3);
  }

  // Explicit index type
  {
    auto Cell = CartesianGridCell1D<int>(1., 2., 6.);
    EXPECT_TRUE((std::is_same<decltype(Cell), int>::value));
    EXPECT_EQ(Cell, 2);
  }

}

TEST_F(GeometryOpsTests, CartesianGridCell2D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::CartesianGridCell2D;

  // Above origin, cell interior
  {
    ovk::tuple<long long> Cell = CartesianGridCell2D({1.,2.,0.}, {2.,4.,0.}, {4.,11.,0.});
    EXPECT_THAT(Cell, ElementsAre(1,2,0));
  }

  // Above origin, lower cell boundary
  {
    ovk::tuple<long long> Cell = CartesianGridCell2D({1.,2.,0.}, {2.,4.,0.}, {3.,10.,0.});
    EXPECT_THAT(Cell, ElementsAre(1,2,0));
  }

  // Above origin, upper cell boundary
  {
    ovk::tuple<long long> Cell = CartesianGridCell2D({1.,2.,0.}, {2.,4.,0.}, {5.,14.,0.});
    EXPECT_THAT(Cell, ElementsAre(2,3,0));
  }

  // Below origin, cell interior
  {
    ovk::tuple<long long> Cell = CartesianGridCell2D({1.,2.,0.}, {2.,4.,0.}, {-4.,-11.,0.});
    EXPECT_THAT(Cell, ElementsAre(-3,-4,0));
  }

  // Below origin, lower cell boundary
  {
    ovk::tuple<long long> Cell = CartesianGridCell2D({1.,2.,0.}, {2.,4.,0.}, {-5.,-14.,0.});
    EXPECT_THAT(Cell, ElementsAre(-3,-4,0));
  }

  // Below origin, upper cell boundary
  {
    ovk::tuple<long long> Cell = CartesianGridCell2D({1.,2.,0.}, {2.,4.,0.}, {-3.,-10.,0.});
    EXPECT_THAT(Cell, ElementsAre(-2,-3,0));
  }

  // Explicit index type
  {
    auto Cell = CartesianGridCell2D<int>({1.,2.,0.}, {2.,4.,0.}, {4.,11.,0.});
    EXPECT_TRUE((std::is_same<decltype(Cell), ovk::tuple<int>>::value));
    EXPECT_THAT(Cell, ElementsAre(1,2,0));
  }

}

TEST_F(GeometryOpsTests, CartesianGridCell3D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::CartesianGridCell3D;

  // Above origin, cell interior
  {
    ovk::tuple<long long> Cell = CartesianGridCell3D({1.,2.,3.}, {2.,4.,6.}, {4.,11.,22.});
    EXPECT_THAT(Cell, ElementsAre(1,2,3));
  }

  // Above origin, lower cell boundary
  {
    ovk::tuple<long long> Cell = CartesianGridCell3D({1.,2.,3.}, {2.,4.,6.}, {3.,10.,21.});
    EXPECT_THAT(Cell, ElementsAre(1,2,3));
  }

  // Above origin, upper cell boundary
  {
    ovk::tuple<long long> Cell = CartesianGridCell3D({1.,2.,3.}, {2.,4.,6.}, {5.,14.,27.});
    EXPECT_THAT(Cell, ElementsAre(2,3,4));
  }

  // Below origin, cell interior
  {
    ovk::tuple<long long> Cell = CartesianGridCell3D({1.,2.,3.}, {2.,4.,6.}, {-4.,-11.,-22.});
    EXPECT_THAT(Cell, ElementsAre(-3,-4,-5));
  }

  // Below origin, lower cell boundary
  {
    ovk::tuple<long long> Cell = CartesianGridCell3D({1.,2.,3.}, {2.,4.,6.}, {-5.,-14.,-27.});
    EXPECT_THAT(Cell, ElementsAre(-3,-4,-5));
  }

  // Below origin, upper cell boundary
  {
    ovk::tuple<long long> Cell = CartesianGridCell3D({1.,2.,3.}, {2.,4.,6.}, {-3.,-10.,-21.});
    EXPECT_THAT(Cell, ElementsAre(-2,-3,-4));
  }

  // Explicit index type
  {
    auto Cell = CartesianGridCell3D<int>({1.,2.,3.}, {2.,4.,6.}, {4.,11.,22.});
    EXPECT_TRUE((std::is_same<decltype(Cell), ovk::tuple<int>>::value));
    EXPECT_THAT(Cell, ElementsAre(1,2,3));
  }

}

TEST_F(GeometryOpsTests, ColumnDeterminant2D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::ColumnDeterminant2D;

  ovk::tuple<double> AI = {
    1.,
   -2.,
    0.
  };
  ovk::tuple<double> AJ = {
    2.,
    1.,
    0.
  };
  EXPECT_EQ(ColumnDeterminant2D(AI, AJ), 5.);

}

TEST_F(GeometryOpsTests, ColumnDeterminant3D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::ColumnDeterminant3D;

  ovk::tuple<double> AI = {
    1.,
   -3.,
   -2.
  };
  ovk::tuple<double> AJ = {
    2.,
    1.,
   -3.
  };
  ovk::tuple<double> AK = {
    3.,
    2.,
    1.
  };
  EXPECT_EQ(ColumnDeterminant3D(AI, AJ, AK), 38.);

}

TEST_F(GeometryOpsTests, ColumnSolve2D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::ColumnSolve2D;

  ovk::tuple<double> AI = {
    1.,
   -2.,
    0.
  };
  ovk::tuple<double> AJ = {
    2.,
    1.,
    0.
  };
  ovk::tuple<double> B = {
    8.,
   -1,
    0.
  };
  EXPECT_THAT(ColumnSolve2D(AI, AJ, B), ElementsAre(2.,3.,0.));

}

TEST_F(GeometryOpsTests, ColumnSolve3D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::ColumnSolve3D;

  ovk::tuple<double> AI = {
    1.,
   -3.,
   -2.
  };
  ovk::tuple<double> AJ = {
    2.,
    1.,
   -3.
  };
  ovk::tuple<double> AK = {
    3.,
    2.,
    1.
  };
  ovk::tuple<double> B = {
   20.,
    5.,
   -9.
  };
  EXPECT_THAT(ColumnSolve3D(AI, AJ, AK, B), ElementsAre(2.,3.,4.));

}

TEST_F(GeometryOpsTests, LagrangeInterpLinear) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::LagrangeInterpLinear;
  using ovk::core::LagrangeInterpLinearDeriv;

  // Point 1
  {
    ovk::elem<double,2> Interp = LagrangeInterpLinear(0.);
    ovk::elem<double,2> InterpDeriv = LagrangeInterpLinearDeriv(0.);
    EXPECT_THAT(Interp, ElementsAre(1.,0.));
    EXPECT_THAT(InterpDeriv, ElementsAre(-1.,1.));
  }

  // Between points 1 and 2
  {
    ovk::elem<double,2> Interp = LagrangeInterpLinear(0.5);
    ovk::elem<double,2> InterpDeriv = LagrangeInterpLinearDeriv(0.5);
    EXPECT_THAT(Interp, ElementsAre(0.5,0.5));
    EXPECT_THAT(InterpDeriv, ElementsAre(-1.,1.));
  }

  // Point 2
  {
    ovk::elem<double,2> Interp = LagrangeInterpLinear(1.);
    ovk::elem<double,2> InterpDeriv = LagrangeInterpLinearDeriv(1.);
    EXPECT_THAT(Interp, ElementsAre(0.,1.));
    EXPECT_THAT(InterpDeriv, ElementsAre(-1.,1.));
  }

}

TEST_F(GeometryOpsTests, LagrangeInterpCubic) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::LagrangeInterpCubic;
  using ovk::core::LagrangeInterpCubicDeriv;

  // Point 1
  {
    ovk::elem<double,4> Interp = LagrangeInterpCubic(-1.);
    ovk::elem<double,4> InterpDeriv = LagrangeInterpCubicDeriv(-1.);
    EXPECT_NEAR(Interp(0), 1., 1.e-12);
    EXPECT_NEAR(Interp(1), 0., 1.e-12);
    EXPECT_NEAR(Interp(2), 0., 1.e-12);
    EXPECT_NEAR(Interp(3), 0., 1.e-12);
    EXPECT_NEAR(InterpDeriv(0), -11./6., 1.e-12);
    EXPECT_NEAR(InterpDeriv(1), 3., 1.e-12);
    EXPECT_NEAR(InterpDeriv(2), -1.5, 1.e-12);
    EXPECT_NEAR(InterpDeriv(3), 1./3., 1.e-12);
  }

  // Between points 1 and 2
  {
    ovk::elem<double,4> Interp = LagrangeInterpCubic(-0.5);
    ovk::elem<double,4> InterpDeriv = LagrangeInterpCubicDeriv(-0.5);
    EXPECT_NEAR(Interp(0), 0.3125, 1.e-12);
    EXPECT_NEAR(Interp(1), 0.9375, 1.e-12);
    EXPECT_NEAR(Interp(2), -0.3125, 1.e-12);
    EXPECT_NEAR(Interp(3), 0.0625, 1.e-12);
    EXPECT_NEAR(InterpDeriv(0), -23./24., 1.e-12);
    EXPECT_NEAR(InterpDeriv(1), 0.875, 1.e-12);
    EXPECT_NEAR(InterpDeriv(2), 0.125, 1.e-12);
    EXPECT_NEAR(InterpDeriv(3), -1./24., 1.e-12);
  }

  // Point 2
  {
    ovk::elem<double,4> Interp = LagrangeInterpCubic(0.);
    ovk::elem<double,4> InterpDeriv = LagrangeInterpCubicDeriv(0.);
    EXPECT_NEAR(Interp(0), 0., 1.e-12);
    EXPECT_NEAR(Interp(1), 1., 1.e-12);
    EXPECT_NEAR(Interp(2), 0., 1.e-12);
    EXPECT_NEAR(Interp(3), 0., 1.e-12);
    EXPECT_NEAR(InterpDeriv(0), -1./3., 1.e-12);
    EXPECT_NEAR(InterpDeriv(1), -0.5, 1.e-12);
    EXPECT_NEAR(InterpDeriv(2), 1., 1.e-12);
    EXPECT_NEAR(InterpDeriv(3), -1./6., 1.e-12);
  }

  // Between points 2 and 3
  {
    ovk::elem<double,4> Interp = LagrangeInterpCubic(0.5);
    ovk::elem<double,4> InterpDeriv = LagrangeInterpCubicDeriv(0.5);
    EXPECT_NEAR(Interp(0), -0.0625, 1.e-12);
    EXPECT_NEAR(Interp(1), 0.5625, 1.e-12);
    EXPECT_NEAR(Interp(2), 0.5625, 1.e-12);
    EXPECT_NEAR(Interp(3), -0.0625, 1.e-12);
    EXPECT_NEAR(InterpDeriv(0), 1./24., 1.e-12);
    EXPECT_NEAR(InterpDeriv(1), -1.125, 1.e-12);
    EXPECT_NEAR(InterpDeriv(2), 1.125, 1.e-12);
    EXPECT_NEAR(InterpDeriv(3), -1./24., 1.e-12);
  }

  // Point 3
  {
    ovk::elem<double,4> Interp = LagrangeInterpCubic(1.);
    ovk::elem<double,4> InterpDeriv = LagrangeInterpCubicDeriv(1.);
    EXPECT_NEAR(Interp(0), 0., 1.e-12);
    EXPECT_NEAR(Interp(1), 0., 1.e-12);
    EXPECT_NEAR(Interp(2), 1., 1.e-12);
    EXPECT_NEAR(Interp(3), 0., 1.e-12);
    EXPECT_NEAR(InterpDeriv(0), 1./6., 1.e-12);
    EXPECT_NEAR(InterpDeriv(1), -1., 1.e-12);
    EXPECT_NEAR(InterpDeriv(2), 0.5, 1.e-12);
    EXPECT_NEAR(InterpDeriv(3), 1./3., 1.e-12);
  }

  // Between points 3 and 4
  {
    ovk::elem<double,4> Interp = LagrangeInterpCubic(1.5);
    ovk::elem<double,4> InterpDeriv = LagrangeInterpCubicDeriv(1.5);
    EXPECT_NEAR(Interp(0), 0.0625, 1.e-12);
    EXPECT_NEAR(Interp(1), -0.3125, 1.e-12);
    EXPECT_NEAR(Interp(2), 0.9375, 1.e-12);
    EXPECT_NEAR(Interp(3), 0.3125, 1.e-12);
    EXPECT_NEAR(InterpDeriv(0), 1./24., 1.e-12);
    EXPECT_NEAR(InterpDeriv(1), -0.125, 1.e-12);
    EXPECT_NEAR(InterpDeriv(2), -0.875, 1.e-12);
    EXPECT_NEAR(InterpDeriv(3),  23./24., 1.e-12);
  }

  // Point 4
  {
    ovk::elem<double,4> Interp = LagrangeInterpCubic(2.);
    ovk::elem<double,4> InterpDeriv = LagrangeInterpCubicDeriv(2.);
    EXPECT_NEAR(Interp(0), 0., 1.e-12);
    EXPECT_NEAR(Interp(1), 0., 1.e-12);
    EXPECT_NEAR(Interp(2), 0., 1.e-12);
    EXPECT_NEAR(Interp(3), 1., 1.e-12);
    EXPECT_NEAR(InterpDeriv(0), -1./3., 1.e-12);
    EXPECT_NEAR(InterpDeriv(1), 1.5, 1.e-12);
    EXPECT_NEAR(InterpDeriv(2), -3., 1.e-12);
    EXPECT_NEAR(InterpDeriv(3), 11./6., 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoLine2Node) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoLine2Node;

  // Line segment from 1 to 3

  auto CoordFunc = [](double U) -> double {
    return 1.+2.*U;
  };

  const double NodeCoords[] = {
    CoordFunc(0.),
    CoordFunc(1.)
  };

  // Sanity check
  EXPECT_EQ(NodeCoords[0], 1.);
  EXPECT_EQ(NodeCoords[1], 3.);

  const double LowerCoord = NodeCoords[0];
  const double UpperCoord = NodeCoords[1];

  // At nodes
  int iNode = 0;
  for (int i = 0; i < 2; ++i) {
    double LocalCoord = double(i);
    double Coord = IsoLine2Node(LowerCoord, UpperCoord, LocalCoord);
    EXPECT_EQ(Coord, NodeCoords[iNode]);
    ++iNode;
  }

  // Between nodes
  {
    double Coord = IsoLine2Node(LowerCoord, UpperCoord, 0.5);
    EXPECT_EQ(Coord, CoordFunc(0.5));
  }

}

TEST_F(GeometryOpsTests, IsoLine2NodeInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoLine2NodeInverse;

  // Line segment from 1 to 3

  auto CoordFunc = [](double U) -> double {
    return 1.+2.*U;
  };

  const double NodeCoords[] = {
    CoordFunc(0.),
    CoordFunc(1.)
  };

  // Sanity check
  EXPECT_EQ(NodeCoords[0], 1.);
  EXPECT_EQ(NodeCoords[1], 3.);

  const double LowerCoord = NodeCoords[0];
  const double UpperCoord = NodeCoords[1];

  // At nodes
  int iNode = 0;
  for (int i = 0; i < 2; ++i) {
    double LocalCoord = IsoLine2NodeInverse(LowerCoord, UpperCoord, NodeCoords[iNode]);
    EXPECT_EQ(LocalCoord, double(i));
    ++iNode;
  }

  // Between nodes
  {
    double LocalCoord = IsoLine2NodeInverse(LowerCoord, UpperCoord, CoordFunc(0.5));
    EXPECT_EQ(LocalCoord, 0.5);
  }

}

TEST_F(GeometryOpsTests, IsoLine4Node) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoLine4Node;

  // Line segment from 1 to 4, stretched

  auto CoordFunc = [](double U) -> double {
    auto Stretch = [](double X) -> double { return X*X*X; };
    return Stretch(2.+U);
  };

  const double NodeCoords[] = {
    CoordFunc(-1.),
    CoordFunc(0.),
    CoordFunc(1.),
    CoordFunc(2.)
  };

  // Sanity check
  EXPECT_EQ(NodeCoords[0], 1.);
  EXPECT_EQ(NodeCoords[1], 8.);
  EXPECT_EQ(NodeCoords[2], 27.);
  EXPECT_EQ(NodeCoords[3], 64.);

  // At nodes
  int iNode = 0;
  for (int i = 0; i < 4; ++i) {
    double LocalCoord = double(i)-1.;
    double Coord = IsoLine4Node(NodeCoords, LocalCoord);
    EXPECT_NEAR(Coord, NodeCoords[iNode], 1.e-12);
    ++iNode;
  }

  // Between nodes
  {
    double Coord = IsoLine4Node(NodeCoords, 0.5);
    EXPECT_NEAR(Coord, CoordFunc(0.5), 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoLine4NodeInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoLine4NodeInverse;

  // Line segment from 1 to 4, stretched

  auto CoordFunc = [](double U) -> double {
    auto Stretch = [](double X) -> double { return X*X*X; };
    return Stretch(2.+U);
  };

  const double NodeCoords[] = {
    CoordFunc(-1.),
    CoordFunc(0.),
    CoordFunc(1.),
    CoordFunc(2.)
  };

  // Sanity check
  EXPECT_EQ(NodeCoords[0], 1.);
  EXPECT_EQ(NodeCoords[1], 8.);
  EXPECT_EQ(NodeCoords[2], 27.);
  EXPECT_EQ(NodeCoords[3], 64.);

  // At nodes
  int iNode = 0;
  for (int i = 0; i < 4; ++i) {
    double LocalCoord = IsoLine4NodeInverse(NodeCoords, NodeCoords[iNode]);
    EXPECT_NEAR(LocalCoord, double(i)-1., 1.e-12);
    ++iNode;
  }

  // Between nodes
  {
    bool Success;
    double LocalCoord = IsoLine4NodeInverse(NodeCoords, CoordFunc(0.5), &Success);
    EXPECT_TRUE(Success);
    EXPECT_NEAR(LocalCoord, 0.5, 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoQuad4NodeUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoQuad4NodeUniform;

  // Rectangle of length (2,3) with lower corner at (1,2)

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    return {1.+2.*U, 2.+3.*V, 0.};
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,0.));

  const ovk::tuple<double> &LowerCoords = NodeCoords[0];
  const ovk::tuple<double> &UpperCoords = NodeCoords[3];

  // At nodes
  int iNode = 0;
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      ovk::tuple<double> LocalCoords = {double(i), double(j), 0.};
      ovk::tuple<double> Coords = IsoQuad4NodeUniform(LowerCoords, UpperCoords, LocalCoords);
      EXPECT_NEAR(Coords(0), NodeCoords[iNode](0), 1.e-12);
      EXPECT_NEAR(Coords(1), NodeCoords[iNode](1), 1.e-12);
      EXPECT_EQ(Coords(2), 0.);
      ++iNode;
    }
  }

  // Between nodes
  {
    ovk::tuple<double> Coords = IsoQuad4NodeUniform(LowerCoords, UpperCoords, {0.5,0.5,0.});
    EXPECT_NEAR(Coords(0), CoordFunc(0.5,0.5)(0), 1.e-12);
    EXPECT_NEAR(Coords(1), CoordFunc(0.5,0.5)(1), 1.e-12);
    EXPECT_EQ(Coords(2), 0.);
  }

}

TEST_F(GeometryOpsTests, IsoQuad4NodeUniformInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoQuad4NodeUniformInverse;

  // Rectangle of length (2,3) with lower corner at (1,2)

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    return {1.+2.*U, 2.+3.*V, 0.};
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,0.));

  const ovk::tuple<double> &LowerCoords = NodeCoords[0];
  const ovk::tuple<double> &UpperCoords = NodeCoords[3];

  // At nodes
  int iNode = 0;
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      ovk::tuple<double> LocalCoords = IsoQuad4NodeUniformInverse(LowerCoords, UpperCoords,
        NodeCoords[iNode]);
      EXPECT_NEAR(LocalCoords(0), double(i), 1.e-12);
      EXPECT_NEAR(LocalCoords(1), double(j), 1.e-12);
      EXPECT_EQ(LocalCoords(2), 0.);
      ++iNode;
    }
  }

  // Between nodes
  {
    ovk::tuple<double> LocalCoords = IsoQuad4NodeUniformInverse(LowerCoords, UpperCoords,
      CoordFunc(0.5,0.5));
    EXPECT_NEAR(LocalCoords(0), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(1), 0.5, 1.e-12);
    EXPECT_EQ(LocalCoords(2), 0.);
  }

}

TEST_F(GeometryOpsTests, IsoQuad4NodeOrientedUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoQuad4NodeOrientedUniform;

  // Rectangle of length (2,3) with lower corner at (1,2), "rotated" 45 degrees cw about
  // lower corner (actually stretched too, but only axially)

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Transform = [](double X, double Y) -> ovk::tuple<double> {
      double XRel = X-1.;
      double YRel = Y-2.;
      return {1.+XRel+YRel, 2.-XRel+YRel, 0.};
    };
    return Transform(1.+2.*U, 2.+3.*V);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,0.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(6.,3.,0.));

  // At nodes
  int iNode = 0;
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      ovk::tuple<double> LocalCoords = {double(i), double(j), 0.};
      ovk::tuple<double> Coords = IsoQuad4NodeOrientedUniform(NodeCoords, LocalCoords);
      EXPECT_NEAR(Coords(0), NodeCoords[iNode](0), 1.e-12);
      EXPECT_NEAR(Coords(1), NodeCoords[iNode](1), 1.e-12);
      EXPECT_EQ(Coords(2), 0.);
      ++iNode;
    }
  }

  // Between nodes
  {
    ovk::tuple<double> Coords = IsoQuad4NodeOrientedUniform(NodeCoords, {0.5,0.5,0.});
    EXPECT_NEAR(Coords(0), CoordFunc(0.5,0.5)(0), 1.e-12);
    EXPECT_NEAR(Coords(1), CoordFunc(0.5,0.5)(1), 1.e-12);
    EXPECT_EQ(Coords(2), 0.);
  }

}

TEST_F(GeometryOpsTests, IsoQuad4NodeOrientedUniformInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoQuad4NodeOrientedUniformInverse;

  // Rectangle of length (2,3) with lower corner at (1,2), "rotated" 45 degrees cw about
  // lower corner (actually stretched too, but only axially)

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Transform = [](double X, double Y) -> ovk::tuple<double> {
      double XRel = X-1.;
      double YRel = Y-2.;
      return {1.+XRel+YRel, 2.-XRel+YRel, 0.};
    };
    return Transform(1.+2.*U, 2.+3.*V);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,0.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(6.,3.,0.));

  // At nodes
  int iNode = 0;
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      ovk::tuple<double> LocalCoords = IsoQuad4NodeOrientedUniformInverse(NodeCoords,
        NodeCoords[iNode]);
      EXPECT_NEAR(LocalCoords(0), double(i), 1.e-12);
      EXPECT_NEAR(LocalCoords(1), double(j), 1.e-12);
      EXPECT_EQ(LocalCoords(2), 0.);
      ++iNode;
    }
  }

  // Between nodes
  {
    ovk::tuple<double> LocalCoords = IsoQuad4NodeOrientedUniformInverse(NodeCoords,
      CoordFunc(0.5,0.5));
    EXPECT_NEAR(LocalCoords(0), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(1), 0.5, 1.e-12);
    EXPECT_EQ(LocalCoords(2), 0.);
  }

}

TEST_F(GeometryOpsTests, IsoQuad4NodeNonUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoQuad4NodeNonUniform;

  // Unit-length square with lower corner at (1,2), stretched

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y) -> ovk::tuple<double> {
      return {2.*X+Y, X+2.*Y, 0.};
    };
    return Stretch(1.+U, 2.+V);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(6.,6.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(5.,7.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(7.,8.,0.));

  // At nodes
  int iNode = 0;
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      ovk::tuple<double> LocalCoords = {double(i), double(j), 0.};
      ovk::tuple<double> Coords = IsoQuad4NodeNonUniform(NodeCoords, LocalCoords);
      EXPECT_NEAR(Coords(0), NodeCoords[iNode](0), 1.e-12);
      EXPECT_NEAR(Coords(1), NodeCoords[iNode](1), 1.e-12);
      EXPECT_EQ(Coords(2), 0.);
      ++iNode;
    }
  }

  // Between nodes
  {
    ovk::tuple<double> Coords = IsoQuad4NodeNonUniform(NodeCoords, {0.5,0.5,0.});
    EXPECT_NEAR(Coords(0), CoordFunc(0.5,0.5)(0), 1.e-12);
    EXPECT_NEAR(Coords(1), CoordFunc(0.5,0.5)(1), 1.e-12);
    EXPECT_EQ(Coords(2), 0.);
  }

}

TEST_F(GeometryOpsTests, IsoQuad4NodeNonUniformInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoQuad4NodeNonUniformInverse;

  // Unit-length square with lower corner at (1,2), stretched

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y) -> ovk::tuple<double> {
      return {2.*X+Y, X+2.*Y, 0.};
    };
    return Stretch(1.+U, 2.+V);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(6.,6.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(5.,7.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(7.,8.,0.));

  // At nodes
  int iNode = 0;
  for (int j = 0; j < 2; ++j) {
    for (int i = 0; i < 2; ++i) {
      ovk::tuple<double> LocalCoords = IsoQuad4NodeNonUniformInverse(NodeCoords, NodeCoords[iNode]);
      EXPECT_NEAR(LocalCoords(0), double(i), 1.e-12);
      EXPECT_NEAR(LocalCoords(1), double(j), 1.e-12);
      EXPECT_EQ(LocalCoords(2), 0.);
      ++iNode;
    }
  }

  // Between nodes
  {
    bool Success;
    ovk::tuple<double> LocalCoords = IsoQuad4NodeNonUniformInverse(NodeCoords, CoordFunc(0.5,0.5),
      &Success);
    EXPECT_TRUE(Success);
    EXPECT_NEAR(LocalCoords(0), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(1), 0.5, 1.e-12);
    EXPECT_EQ(LocalCoords(2), 0.);
  }

}

TEST_F(GeometryOpsTests, IsoQuad16Node) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoQuad16Node;

  // Square of side length 3 with lower corner at (1,2), stretched

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y) -> ovk::tuple<double> {
      return {2.*X+Y, X+2.*Y, 0.};
    };
    return Stretch(2.+U, 3.+V);
  };

  ovk::tuple<double> NodeCoords[16];
  int iNode = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      NodeCoords[iNode] = CoordFunc(double(i)-1.,double(j)-1.);
      ++iNode;
    }
  }

  // Sanity check (too many nodes -- just pick a few)
  EXPECT_THAT(NodeCoords[0], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(10.,8.,0.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(7.,8.,0.));
  EXPECT_THAT(NodeCoords[12], ElementsAre(7.,11.,0.));
  EXPECT_THAT(NodeCoords[15], ElementsAre(13.,14.,0.));

  // At nodes
  iNode = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      ovk::tuple<double> LocalCoords = {double(i)-1., double(j)-1., 0.};
      ovk::tuple<double> Coords = IsoQuad16Node(NodeCoords, LocalCoords);
      EXPECT_NEAR(Coords(0), NodeCoords[iNode](0), 1.e-12);
      EXPECT_NEAR(Coords(1), NodeCoords[iNode](1), 1.e-12);
      EXPECT_EQ(Coords(2), 0.);
      ++iNode;
    }
  }

  // Between nodes
  {
    ovk::tuple<double> Coords = IsoQuad16Node(NodeCoords, {0.5,0.5,0.});
    EXPECT_NEAR(Coords(0), CoordFunc(0.5,0.5)(0), 1.e-12);
    EXPECT_NEAR(Coords(1), CoordFunc(0.5,0.5)(1), 1.e-12);
    EXPECT_EQ(Coords(2), 0.);
  }

}

TEST_F(GeometryOpsTests, IsoQuad16NodeInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoQuad16NodeInverse;

  // Square of side length 3 with lower corner at (1,2), stretched

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y) -> ovk::tuple<double> {
      return {2.*X+Y, X+2.*Y, 0.};
    };
    return Stretch(2.+U, 3.+V);
  };

  ovk::tuple<double> NodeCoords[16];
  int iNode = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      NodeCoords[iNode] = CoordFunc(double(i)-1.,double(j)-1.);
      ++iNode;
    }
  }

  // Sanity check (too many nodes -- just pick a few)
  EXPECT_THAT(NodeCoords[0], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(10.,8.,0.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(7.,8.,0.));
  EXPECT_THAT(NodeCoords[12], ElementsAre(7.,11.,0.));
  EXPECT_THAT(NodeCoords[15], ElementsAre(13.,14.,0.));

  // At nodes
  iNode = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      ovk::tuple<double> LocalCoords = IsoQuad16NodeInverse(NodeCoords, NodeCoords[iNode]);
      EXPECT_NEAR(LocalCoords(0), double(i)-1., 1.e-12);
      EXPECT_NEAR(LocalCoords(1), double(j)-1., 1.e-12);
      EXPECT_EQ(LocalCoords(2), 0.);
      ++iNode;
    }
  }

  // Between nodes
  {
    bool Success;
    ovk::tuple<double> LocalCoords = IsoQuad16NodeInverse(NodeCoords, CoordFunc(0.5,0.5), &Success);
    EXPECT_TRUE(Success);
    EXPECT_NEAR(LocalCoords(0), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(1), 0.5, 1.e-12);
    EXPECT_EQ(LocalCoords(2), 0.);
  }

}

TEST_F(GeometryOpsTests, IsoHex8NodeUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoHex8NodeUniform;

  // Cuboid of length (2,3,4) with lower corner at (1,2,3)

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    return {1.+2.*U, 2.+3.*V, 3.+4.*W};
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,3.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,3.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,3.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,3.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(1.,2.,7.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(3.,2.,7.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(1.,5.,7.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(3.,5.,7.));

  const ovk::tuple<double> &LowerCoords = NodeCoords[0];
  const ovk::tuple<double> &UpperCoords = NodeCoords[7];

  // At nodes
  int iNode = 0;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        ovk::tuple<double> LocalCoords = {double(i), double(j), double(k)};
        ovk::tuple<double> Coords = IsoHex8NodeUniform(LowerCoords, UpperCoords, LocalCoords);
        EXPECT_NEAR(Coords(0), NodeCoords[iNode](0), 1.e-12);
        EXPECT_NEAR(Coords(1), NodeCoords[iNode](1), 1.e-12);
        EXPECT_NEAR(Coords(2), NodeCoords[iNode](2), 1.e-12);
        ++iNode;
      }
    }
  }

  // Between nodes
  {
    ovk::tuple<double> Coords = IsoHex8NodeUniform(LowerCoords, UpperCoords, {0.5,0.5,0.5});
    EXPECT_NEAR(Coords(0), CoordFunc(0.5,0.5,0.5)(0), 1.e-12);
    EXPECT_NEAR(Coords(1), CoordFunc(0.5,0.5,0.5)(1), 1.e-12);
    EXPECT_NEAR(Coords(2), CoordFunc(0.5,0.5,0.5)(2), 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoHex8NodeUniformInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoHex8NodeUniformInverse;

  // Cuboid of length (2,3,4) with lower corner at (1,2,3)

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    return {1.+2.*U, 2.+3.*V, 3.+4.*W};
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,3.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,3.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,3.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,3.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(1.,2.,7.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(3.,2.,7.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(1.,5.,7.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(3.,5.,7.));

  const ovk::tuple<double> &LowerCoords = NodeCoords[0];
  const ovk::tuple<double> &UpperCoords = NodeCoords[7];

  // At nodes
  int iNode = 0;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        ovk::tuple<double> LocalCoords = IsoHex8NodeUniformInverse(LowerCoords, UpperCoords,
          NodeCoords[iNode]);
        EXPECT_NEAR(LocalCoords(0), double(i), 1.e-12);
        EXPECT_NEAR(LocalCoords(1), double(j), 1.e-12);
        EXPECT_NEAR(LocalCoords(2), double(k), 1.e-12);
        ++iNode;
      }
    }
  }

  // Between nodes
  {
    ovk::tuple<double> LocalCoords = IsoHex8NodeUniformInverse(LowerCoords, UpperCoords,
      CoordFunc(0.5,0.5,0.5));
    EXPECT_NEAR(LocalCoords(0), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(1), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(2), 0.5, 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoHex8NodeOrientedUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoHex8NodeOrientedUniform;

  // Cuboid of length (2,3,4) with lower corner at (1,2,3), "rotated" 45 degrees cw about x through
  // lower corner (actually stretched too, but only axially)

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Transform = [](double X, double Y, double Z) -> ovk::tuple<double> {
      double YRel = Y-2.;
      double ZRel = Z-3.;
      return {X, 2.+YRel+ZRel, 3.-YRel+ZRel};
    };
    return Transform(1.+2.*U, 2.+3.*V, 3.+4.*W);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,3.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,3.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,0.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(1.,6.,7.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(3.,6.,7.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(1.,9.,4.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(3.,9.,4.));

  // At nodes
  int iNode = 0;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        ovk::tuple<double> LocalCoords = {double(i), double(j), double(k)};
        ovk::tuple<double> Coords = IsoHex8NodeOrientedUniform(NodeCoords, LocalCoords);
        EXPECT_NEAR(Coords(0), NodeCoords[iNode](0), 1.e-12);
        EXPECT_NEAR(Coords(1), NodeCoords[iNode](1), 1.e-12);
        EXPECT_NEAR(Coords(2), NodeCoords[iNode](2), 1.e-12);
        ++iNode;
      }
    }
  }

  // Between nodes
  {
    ovk::tuple<double> Coords = IsoHex8NodeOrientedUniform(NodeCoords, {0.5,0.5,0.5});
    EXPECT_NEAR(Coords(0), CoordFunc(0.5,0.5,0.5)(0), 1.e-12);
    EXPECT_NEAR(Coords(1), CoordFunc(0.5,0.5,0.5)(1), 1.e-12);
    EXPECT_NEAR(Coords(2), CoordFunc(0.5,0.5,0.5)(2), 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoHex8NodeOrientedUniformInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoHex8NodeOrientedUniformInverse;

  // Cuboid of length (2,3,4) with lower corner at (1,2,3), "rotated" 45 degrees cw about x through
  // lower corner (actually stretched too, but only axially)

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Transform = [](double X, double Y, double Z) -> ovk::tuple<double> {
      double YRel = Y-2.;
      double ZRel = Z-3.;
      return {X, 2.+YRel+ZRel, 3.-YRel+ZRel};
    };
    return Transform(1.+2.*U, 2.+3.*V, 3.+4.*W);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,3.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,3.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,0.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(1.,6.,7.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(3.,6.,7.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(1.,9.,4.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(3.,9.,4.));

  // At nodes
  int iNode = 0;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        ovk::tuple<double> LocalCoords = IsoHex8NodeOrientedUniformInverse(NodeCoords,
          NodeCoords[iNode]);
        EXPECT_NEAR(LocalCoords(0), double(i), 1.e-12);
        EXPECT_NEAR(LocalCoords(1), double(j), 1.e-12);
        EXPECT_NEAR(LocalCoords(2), double(k), 1.e-12);
        ++iNode;
      }
    }
  }

  // Between nodes
  {
    ovk::tuple<double> LocalCoords = IsoHex8NodeOrientedUniformInverse(NodeCoords,
      CoordFunc(0.5,0.5,0.5));
    EXPECT_NEAR(LocalCoords(0), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(1), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(2), 0.5, 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoHex8NodeNonUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoHex8NodeNonUniform;

  // Unit-length cube with lower corner at (1,2,3), stretched

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y, double Z) -> ovk::tuple<double> {
      return {2.*X+Y+Z, X+2.*Y+Z, X+Y+2.*Z};
    };
    return Stretch(1.+U, 2.+V, 3.+W);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(7.,8.,9.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(9.,9.,10.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(8.,10.,10.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(10.,11.,11.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(8.,9.,11.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(10.,10.,12.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(9.,11.,12.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(11.,12.,13.));

  // At nodes
  int iNode = 0;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        ovk::tuple<double> LocalCoords = {double(i), double(j), double(k)};
        ovk::tuple<double> Coords = IsoHex8NodeNonUniform(NodeCoords, LocalCoords);
        EXPECT_NEAR(Coords(0), NodeCoords[iNode](0), 1.e-12);
        EXPECT_NEAR(Coords(1), NodeCoords[iNode](1), 1.e-12);
        EXPECT_NEAR(Coords(2), NodeCoords[iNode](2), 1.e-12);
        ++iNode;
      }
    }
  }

  // Between nodes
  {
    ovk::tuple<double> Coords = IsoHex8NodeNonUniform(NodeCoords, {0.5,0.5,0.5});
    EXPECT_NEAR(Coords(0), CoordFunc(0.5,0.5,0.5)(0), 1.e-12);
    EXPECT_NEAR(Coords(1), CoordFunc(0.5,0.5,0.5)(1), 1.e-12);
    EXPECT_NEAR(Coords(2), CoordFunc(0.5,0.5,0.5)(2), 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoHex8NodeNonUniformInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoHex8NodeNonUniformInverse;

  // Unit-length cube with lower corner at (1,2,3), stretched

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y, double Z) -> ovk::tuple<double> {
      return {2.*X+Y+Z, X+2.*Y+Z, X+Y+2.*Z};
    };
    return Stretch(1.+U, 2.+V, 3.+W);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(7.,8.,9.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(9.,9.,10.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(8.,10.,10.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(10.,11.,11.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(8.,9.,11.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(10.,10.,12.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(9.,11.,12.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(11.,12.,13.));

  // At nodes
  int iNode = 0;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        ovk::tuple<double> LocalCoords = IsoHex8NodeNonUniformInverse(NodeCoords,
          NodeCoords[iNode]);
        EXPECT_NEAR(LocalCoords(0), double(i), 1.e-12);
        EXPECT_NEAR(LocalCoords(1), double(j), 1.e-12);
        EXPECT_NEAR(LocalCoords(2), double(k), 1.e-12);
        ++iNode;
      }
    }
  }

  // Between nodes
  {
    bool Success;
    ovk::tuple<double> LocalCoords = IsoHex8NodeNonUniformInverse(NodeCoords,
      CoordFunc(0.5,0.5,0.5), &Success);
    EXPECT_TRUE(Success);
    EXPECT_NEAR(LocalCoords(0), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(1), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(2), 0.5, 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoHex64Node) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoHex64Node;

  // Cube of side length 3 with lower corner at (1,2,3), stretched

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y, double Z) -> ovk::tuple<double> {
      return {2.*X+Y+Z, X+2.*Y+Z, X+Y+2.*Z};
    };
    return Stretch(2.+U, 3.+V, 4.+W);
  };

  ovk::tuple<double> NodeCoords[64];
  int iNode = 0;
  for (int k = 0; k < 4; ++k) {
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 4; ++i) {
        NodeCoords[iNode] = CoordFunc(double(i)-1.,double(j)-1.,double(k)-1.);
        ++iNode;
      }
    }
  }

  // Sanity check (too many nodes -- just pick a few)
  EXPECT_THAT(NodeCoords[0], ElementsAre(7.,8.,9.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(13.,11.,12.));
  EXPECT_THAT(NodeCoords[12], ElementsAre(10.,14.,12.));
  EXPECT_THAT(NodeCoords[15], ElementsAre(16.,17.,15.));
  EXPECT_THAT(NodeCoords[21], ElementsAre(11.,12.,13.));
  EXPECT_THAT(NodeCoords[48], ElementsAre(10.,11.,15.));
  EXPECT_THAT(NodeCoords[51], ElementsAre(16.,14.,18.));
  EXPECT_THAT(NodeCoords[60], ElementsAre(13.,17.,18.));
  EXPECT_THAT(NodeCoords[63], ElementsAre(19.,20.,21.));

  // At nodes
  iNode = 0;
  for (int k = 0; k < 4; ++k) {
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 4; ++i) {
        ovk::tuple<double> LocalCoords = {double(i)-1., double(j)-1., double(k)-1.};
        ovk::tuple<double> Coords = IsoHex64Node(NodeCoords, LocalCoords);
        EXPECT_NEAR(Coords(0), NodeCoords[iNode](0), 1.e-12);
        EXPECT_NEAR(Coords(1), NodeCoords[iNode](1), 1.e-12);
        EXPECT_NEAR(Coords(2), NodeCoords[iNode](2), 1.e-12);
        ++iNode;
      }
    }
  }

  // Between nodes
  {
    ovk::tuple<double> Coords = IsoHex64Node(NodeCoords, {0.5,0.5,0.5});
    EXPECT_NEAR(Coords(0), CoordFunc(0.5,0.5,0.5)(0), 1.e-12);
    EXPECT_NEAR(Coords(1), CoordFunc(0.5,0.5,0.5)(1), 1.e-12);
    EXPECT_NEAR(Coords(2), CoordFunc(0.5,0.5,0.5)(2), 1.e-12);
  }

}

TEST_F(GeometryOpsTests, IsoHex64NodeInverse) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::IsoHex64NodeInverse;

  // Cube of side length 3 with lower corner at (1,2,3), stretched

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y, double Z) -> ovk::tuple<double> {
      return {2.*X+Y+Z, X+2.*Y+Z, X+Y+2.*Z};
    };
    return Stretch(2.+U, 3.+V, 4.+W);
  };

  ovk::tuple<double> NodeCoords[64];
  int iNode = 0;
  for (int k = 0; k < 4; ++k) {
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 4; ++i) {
        NodeCoords[iNode] = CoordFunc(double(i)-1.,double(j)-1.,double(k)-1.);
        ++iNode;
      }
    }
  }

  // Sanity check (too many nodes -- just pick a few)
  EXPECT_THAT(NodeCoords[0], ElementsAre(7.,8.,9.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(13.,11.,12.));
  EXPECT_THAT(NodeCoords[12], ElementsAre(10.,14.,12.));
  EXPECT_THAT(NodeCoords[15], ElementsAre(16.,17.,15.));
  EXPECT_THAT(NodeCoords[21], ElementsAre(11.,12.,13.));
  EXPECT_THAT(NodeCoords[48], ElementsAre(10.,11.,15.));
  EXPECT_THAT(NodeCoords[51], ElementsAre(16.,14.,18.));
  EXPECT_THAT(NodeCoords[60], ElementsAre(13.,17.,18.));
  EXPECT_THAT(NodeCoords[63], ElementsAre(19.,20.,21.));

  // At nodes
  iNode = 0;
  for (int k = 0; k < 4; ++k) {
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 4; ++i) {
        ovk::tuple<double> LocalCoords = IsoHex64NodeInverse(NodeCoords, NodeCoords[iNode]);
        EXPECT_NEAR(LocalCoords(0), double(i)-1., 1.e-12);
        EXPECT_NEAR(LocalCoords(1), double(j)-1., 1.e-12);
        EXPECT_NEAR(LocalCoords(2), double(k)-1., 1.e-12);
        ++iNode;
      }
    }
  }

  // Between nodes
  {
    bool Success;
    ovk::tuple<double> LocalCoords = IsoHex64NodeInverse(NodeCoords, CoordFunc(0.5,0.5,0.5),
      &Success);
    EXPECT_TRUE(Success);
    EXPECT_NEAR(LocalCoords(0), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(1), 0.5, 1.e-12);
    EXPECT_NEAR(LocalCoords(2), 0.5, 1.e-12);
  }

}

TEST_F(GeometryOpsTests, OverlapsLine) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::OverlapsLine;

  // Line segment from 1 to 3

  auto CoordFunc = [](double U) -> double {
    return 1.+2.*U;
  };

  const double NodeCoords[] = {
    CoordFunc(0.),
    CoordFunc(1.)
  };

  // Sanity check
  EXPECT_EQ(NodeCoords[0], 1.);
  EXPECT_EQ(NodeCoords[1], 3.);

  const double LowerCoord = NodeCoords[0];
  const double UpperCoord = NodeCoords[1];

  // Inside
  EXPECT_TRUE(OverlapsLine(LowerCoord, UpperCoord, CoordFunc(0.5)));

  // Outside, lower
  EXPECT_FALSE(OverlapsLine(LowerCoord, UpperCoord, CoordFunc(-0.5)));

  // Outside, upper
  EXPECT_FALSE(OverlapsLine(LowerCoord, UpperCoord, CoordFunc(1.5)));

  // Boundary, lower
  EXPECT_TRUE(OverlapsLine(LowerCoord, UpperCoord, CoordFunc(0.)));

  // Boundary, upper
  EXPECT_TRUE(OverlapsLine(LowerCoord, UpperCoord, CoordFunc(1.)));

  // Within tolerance, lower
  EXPECT_TRUE(OverlapsLine(LowerCoord, UpperCoord, CoordFunc(-0.01), 0.02));

  // Within tolerance, upper
  EXPECT_TRUE(OverlapsLine(LowerCoord, UpperCoord, CoordFunc(1.01), 0.02));

}

TEST_F(GeometryOpsTests, OverlapsQuadUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::OverlapsQuadUniform;

  // Rectangle of length (2,3) with lower corner at (1,2)

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    return {1.+2.*U, 2.+3.*V, 0.};
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,0.));

  const ovk::tuple<double> &LowerCoords = NodeCoords[0];
  const ovk::tuple<double> &UpperCoords = NodeCoords[3];

  // Inside
  EXPECT_TRUE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.5)));

  // Outside, x, lower
  EXPECT_FALSE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(-0.5,0.5)));

  // Outside, x, upper
  EXPECT_FALSE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(1.5,0.5)));

  // Outside, y, lower
  EXPECT_FALSE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(0.5,-0.5)));

  // Outside, y, upper
  EXPECT_FALSE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(0.5,1.5)));

  // Boundary, x, lower
  EXPECT_TRUE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(0.,0.5)));

  // Boundary, x, upper
  EXPECT_TRUE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(1.,0.5)));

  // Boundary, y, lower
  EXPECT_TRUE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.)));

  // Boundary, y, upper
  EXPECT_TRUE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(0.5,1.)));

  // Within tolerance, x, lower
  EXPECT_TRUE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(-0.01,0.5), 0.02));

  // Within tolerance, x, upper
  EXPECT_TRUE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(1.01,0.5), 0.02));

  // Within tolerance, y, lower
  EXPECT_TRUE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(0.5,-0.01), 0.02));

  // Within tolerance, y, upper
  EXPECT_TRUE(OverlapsQuadUniform(LowerCoords, UpperCoords, CoordFunc(0.5,1.01), 0.02));

}

TEST_F(GeometryOpsTests, OverlapsQuadOrientedUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::OverlapsQuadOrientedUniform;

  // Rectangle of length (2,3) with lower corner at (1,2), "rotated" 45 degrees cw about
  // lower corner (actually stretched too, but only axially)

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Transform = [](double X, double Y) -> ovk::tuple<double> {
      double XRel = X-1.;
      double YRel = Y-2.;
      return {1.+XRel+YRel, 2.-XRel+YRel, 0.};
    };
    return Transform(1.+2.*U, 2.+3.*V);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,0.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(6.,3.,0.));

  // Inside
  EXPECT_TRUE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(0.5,0.5)));

  // Outside, I, lower
  EXPECT_FALSE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(-0.5,0.5)));

  // Outside, I, upper
  EXPECT_FALSE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(1.5,0.5)));

  // Outside, J, lower
  EXPECT_FALSE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(0.5,-0.5)));

  // Outside, J, upper
  EXPECT_FALSE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(0.5,1.5)));

  // Boundary, I, lower
  EXPECT_TRUE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(0.,0.5)));

  // Boundary, I, upper
  EXPECT_TRUE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(1.,0.5)));

  // Boundary, J, lower
  EXPECT_TRUE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(0.5,0.)));

  // Boundary, J, upper
  EXPECT_TRUE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(0.5,1.)));

  // Within tolerance, I, lower
  EXPECT_TRUE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(-0.01,0.5), 0.02));

  // Within tolerance, I, upper
  EXPECT_TRUE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(1.01,0.5), 0.02));

  // Within tolerance, J, lower
  EXPECT_TRUE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(0.5,-0.01), 0.02));

  // Within tolerance, J, upper
  EXPECT_TRUE(OverlapsQuadOrientedUniform(NodeCoords, CoordFunc(0.5,1.01), 0.02));

}

TEST_F(GeometryOpsTests, OverlapsQuadNonUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::OverlapsQuadNonUniform;

  // Unit-length square with lower corner at (1,2), stretched

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y) -> ovk::tuple<double> {
      return {2.*X+Y, X+2.*Y, 0.};
    };
    return Stretch(1.+U, 2.+V);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(6.,6.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(5.,7.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(7.,8.,0.));

  // Inside
  EXPECT_TRUE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(0.5,0.5)));

  // Outside, I, lower
  EXPECT_FALSE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(-0.5,0.5)));

  // Outside, I, upper
  EXPECT_FALSE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(1.5,0.5)));

  // Outside, J, lower
  EXPECT_FALSE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(0.5,-0.5)));

  // Outside, J, upper
  EXPECT_FALSE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(0.5,1.5)));

  // Boundary, I, lower
  EXPECT_TRUE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(0.,0.5)));

  // Boundary, I, upper
  EXPECT_TRUE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(1.,0.5)));

  // Boundary, J, lower
  EXPECT_TRUE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(0.5,0.)));

  // Boundary, J, upper
  EXPECT_TRUE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(0.5,1.)));

  // Within tolerance, I, lower
  EXPECT_TRUE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(-0.01,0.5), 0.02));

  // Within tolerance, I, upper
  EXPECT_TRUE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(1.01,0.5), 0.02));

  // Within tolerance, J, lower
  EXPECT_TRUE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(0.5,-0.01), 0.02));

  // Within tolerance, J, upper
  EXPECT_TRUE(OverlapsQuadNonUniform(NodeCoords, CoordFunc(0.5,1.01), 0.02));

}

TEST_F(GeometryOpsTests, OverlapsHexUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::OverlapsHexUniform;

  // Cuboid of length (2,3,4) with lower corner at (1,2,3)

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    return {1.+2.*U, 2.+3.*V, 3.+4.*W};
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,3.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,3.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,3.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,3.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(1.,2.,7.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(3.,2.,7.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(1.,5.,7.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(3.,5.,7.));

  const ovk::tuple<double> &LowerCoords = NodeCoords[0];
  const ovk::tuple<double> &UpperCoords = NodeCoords[7];

  // Inside
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.5,0.5)));

  // Outside, x, lower
  EXPECT_FALSE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(-0.5,0.5,0.5)));

  // Outside, x, upper
  EXPECT_FALSE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(1.5,0.5,0.5)));

  // Outside, y, lower
  EXPECT_FALSE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,-0.5,0.5)));

  // Outside, y, upper
  EXPECT_FALSE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,1.5,0.5)));

  // Outside, z, lower
  EXPECT_FALSE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.5,-0.5)));

  // Outside, z, upper
  EXPECT_FALSE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.5,1.5)));

  // Boundary, x, lower
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.,0.5,0.5)));

  // Boundary, x, upper
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(1.,0.5,0.5)));

  // Boundary, y, lower
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.,0.5)));

  // Boundary, y, upper
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,1.,0.5)));

  // Boundary, z, lower
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.5,0.)));

  // Boundary, z, upper
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.5,1.)));

  // Within tolerance, x, lower
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(-0.01,0.5,0.5), 0.02));

  // Within tolerance, x, upper
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(1.01,0.5,0.5), 0.02));

  // Within tolerance, y, lower
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,-0.01,0.5), 0.02));

  // Within tolerance, y, upper
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,1.01,0.5), 0.02));

  // Within tolerance, z, lower
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.5,-0.01), 0.02));

  // Within tolerance, z, upper
  EXPECT_TRUE(OverlapsHexUniform(LowerCoords, UpperCoords, CoordFunc(0.5,0.5,1.01), 0.02));

}

TEST_F(GeometryOpsTests, OverlapsHexOrientedUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::OverlapsHexOrientedUniform;

  // Cuboid of length (2,3,4) with lower corner at (1,2,3), "rotated" 45 degrees cw about x through
  // lower corner (actually stretched too, but only axially)

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Transform = [](double X, double Y, double Z) -> ovk::tuple<double> {
      double YRel = Y-2.;
      double ZRel = Z-3.;
      return {X, 2.+YRel+ZRel, 3.-YRel+ZRel};
    };
    return Transform(1.+2.*U, 2.+3.*V, 3.+4.*W);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,3.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,3.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,0.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(1.,6.,7.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(3.,6.,7.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(1.,9.,4.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(3.,9.,4.));

  // Inside
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,0.5,0.5)));

  // Outside, I, lower
  EXPECT_FALSE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(-0.5,0.5,0.5)));

  // Outside, I, upper
  EXPECT_FALSE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(1.5,0.5,0.5)));

  // Outside, J, lower
  EXPECT_FALSE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,-0.5,0.5)));

  // Outside, J, upper
  EXPECT_FALSE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,1.5,0.5)));

  // Outside, K, lower
  EXPECT_FALSE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,0.5,-0.5)));

  // Outside, K, upper
  EXPECT_FALSE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,0.5,1.5)));

  // Boundary, I, lower
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.,0.5,0.5)));

  // Boundary, I, upper
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(1.,0.5,0.5)));

  // Boundary, J, lower
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,0.,0.5)));

  // Boundary, J, upper
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,1.,0.5)));

  // Boundary, K, lower
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,0.5,0.)));

  // Boundary, K, upper
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,0.5,1.)));

  // Within tolerance, I, lower
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(-0.01,0.5,0.5), 0.02));

  // Within tolerance, I, upper
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(1.01,0.5,0.5), 0.02));

  // Within tolerance, J, lower
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,-0.01,0.5), 0.02));

  // Within tolerance, J, upper
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,1.01,0.5), 0.02));

  // Within tolerance, K, lower
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,0.5,-0.01), 0.02));

  // Within tolerance, K, upper
  EXPECT_TRUE(OverlapsHexOrientedUniform(NodeCoords, CoordFunc(0.5,0.5,1.01), 0.02));

}

TEST_F(GeometryOpsTests, OverlapsHexNonUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::OverlapsHexNonUniform;

  // Unit-length cube with lower corner at (1,2,3), stretched

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y, double Z) -> ovk::tuple<double> {
      return {2.*X+Y+Z, X+2.*Y+Z, X+Y+2.*Z};
    };
    return Stretch(1.+U, 2.+V, 3.+W);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(7.,8.,9.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(9.,9.,10.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(8.,10.,10.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(10.,11.,11.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(8.,9.,11.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(10.,10.,12.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(9.,11.,12.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(11.,12.,13.));

  // Inside
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,0.5,0.5)));

  // Outside, I, lower
  EXPECT_FALSE(OverlapsHexNonUniform(NodeCoords, CoordFunc(-0.5,0.5,0.5)));

  // Outside, I, upper
  EXPECT_FALSE(OverlapsHexNonUniform(NodeCoords, CoordFunc(1.5,0.5,0.5)));

  // Outside, J, lower
  EXPECT_FALSE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,-0.5,0.5)));

  // Outside, J, upper
  EXPECT_FALSE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,1.5,0.5)));

  // Outside, K, lower
  EXPECT_FALSE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,0.5,-0.5)));

  // Outside, K, upper
  EXPECT_FALSE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,0.5,1.5)));

  // Boundary, I, lower
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.,0.5,0.5)));

  // Boundary, I, upper
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(1.,0.5,0.5)));

  // Boundary, J, lower
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,0.,0.5)));

  // Boundary, J, upper
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,1.,0.5)));

  // Boundary, K, lower
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,0.5,0.)));

  // Boundary, K, upper
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,0.5,1.)));

  // Within tolerance, I, lower
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(-0.01,0.5,0.5), 0.02));

  // Within tolerance, I, upper
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(1.01,0.5,0.5), 0.02));

  // Within tolerance, J, lower
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,-0.01,0.5), 0.02));

  // Within tolerance, J, upper
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,1.01,0.5), 0.02));

  // Within tolerance, K, lower
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,0.5,-0.01), 0.02));

  // Within tolerance, K, upper
  EXPECT_TRUE(OverlapsHexNonUniform(NodeCoords, CoordFunc(0.5,0.5,1.01), 0.02));

}






















TEST_F(GeometryOpsTests, VolumeLine) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::VolumeLine;

  // Line segment from 1 to 3

  auto CoordFunc = [](double U) -> double {
    return 1.+2.*U;
  };

  const double NodeCoords[] = {
    CoordFunc(0.),
    CoordFunc(1.)
  };

  // Sanity check
  EXPECT_EQ(NodeCoords[0], 1.);
  EXPECT_EQ(NodeCoords[1], 3.);

  const double LowerCoord = NodeCoords[0];
  const double UpperCoord = NodeCoords[1];

  EXPECT_EQ(VolumeLine(LowerCoord, UpperCoord), 2.);

}

TEST_F(GeometryOpsTests, VolumeQuadUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::VolumeQuadUniform;

  // Rectangle of length (2,3) with lower corner at (1,2)

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    return {1.+2.*U, 2.+3.*V, 0.};
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,0.));

  const ovk::tuple<double> &LowerCoords = NodeCoords[0];
  const ovk::tuple<double> &UpperCoords = NodeCoords[3];

  EXPECT_EQ(VolumeQuadUniform(LowerCoords, UpperCoords), 6.);

}

TEST_F(GeometryOpsTests, VolumeQuadOrientedUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::VolumeQuadOrientedUniform;

  // Rectangle of length (2,3) with lower corner at (1,2), "rotated" 45 degrees cw about
  // lower corner (actually stretched too, but only axially)

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Transform = [](double X, double Y) -> ovk::tuple<double> {
      double XRel = X-1.;
      double YRel = Y-2.;
      return {1.+XRel+YRel, 2.-XRel+YRel, 0.};
    };
    return Transform(1.+2.*U, 2.+3.*V);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,0.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(6.,3.,0.));

  EXPECT_NEAR(VolumeQuadOrientedUniform(NodeCoords), 12., 1.e-12);

}

TEST_F(GeometryOpsTests, VolumeQuadNonUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::VolumeQuadNonUniform;

  // Unit-length square with lower corner at (1,2), stretched

  auto CoordFunc = [](double U, double V) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y) -> ovk::tuple<double> {
      return {2.*X+Y, X+2.*Y, 0.};
    };
    return Stretch(1.+U, 2.+V);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.),
    CoordFunc(1.,0.),
    CoordFunc(0.,1.),
    CoordFunc(1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(4.,5.,0.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(6.,6.,0.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(5.,7.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(7.,8.,0.));

  EXPECT_NEAR(VolumeQuadNonUniform(NodeCoords), 3., 1.e-12);

}

TEST_F(GeometryOpsTests, VolumeHexUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::VolumeHexUniform;

  // Cuboid of length (2,3,4) with lower corner at (1,2,3)

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    return {1.+2.*U, 2.+3.*V, 3.+4.*W};
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,3.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,3.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,3.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,3.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(1.,2.,7.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(3.,2.,7.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(1.,5.,7.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(3.,5.,7.));

  const ovk::tuple<double> &LowerCoords = NodeCoords[0];
  const ovk::tuple<double> &UpperCoords = NodeCoords[7];

  EXPECT_EQ(VolumeHexUniform(LowerCoords, UpperCoords), 24.);

}

TEST_F(GeometryOpsTests, VolumeHexOrientedUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::VolumeHexOrientedUniform;

  // Cuboid of length (2,3,4) with lower corner at (1,2,3), "rotated" 45 degrees cw about x through
  // lower corner (actually stretched too, but only axially)

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Transform = [](double X, double Y, double Z) -> ovk::tuple<double> {
      double YRel = Y-2.;
      double ZRel = Z-3.;
      return {X, 2.+YRel+ZRel, 3.-YRel+ZRel};
    };
    return Transform(1.+2.*U, 2.+3.*V, 3.+4.*W);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(1.,2.,3.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(3.,2.,3.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(1.,5.,0.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(3.,5.,0.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(1.,6.,7.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(3.,6.,7.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(1.,9.,4.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(3.,9.,4.));

  EXPECT_NEAR(VolumeHexOrientedUniform(NodeCoords), 48., 1.e-12);

}

TEST_F(GeometryOpsTests, VolumeHexNonUniform) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::VolumeHexNonUniform;

  // Unit-length cube with lower corner at (1,2,3), stretched

  auto CoordFunc = [](double U, double V, double W) -> ovk::tuple<double> {
    auto Stretch = [](double X, double Y, double Z) -> ovk::tuple<double> {
      return {2.*X+Y+Z, X+2.*Y+Z, X+Y+2.*Z};
    };
    return Stretch(1.+U, 2.+V, 3.+W);
  };

  const ovk::tuple<double> NodeCoords[] = {
    CoordFunc(0.,0.,0.),
    CoordFunc(1.,0.,0.),
    CoordFunc(0.,1.,0.),
    CoordFunc(1.,1.,0.),
    CoordFunc(0.,0.,1.),
    CoordFunc(1.,0.,1.),
    CoordFunc(0.,1.,1.),
    CoordFunc(1.,1.,1.)
  };

  // Sanity check
  EXPECT_THAT(NodeCoords[0], ElementsAre(7.,8.,9.));
  EXPECT_THAT(NodeCoords[1], ElementsAre(9.,9.,10.));
  EXPECT_THAT(NodeCoords[2], ElementsAre(8.,10.,10.));
  EXPECT_THAT(NodeCoords[3], ElementsAre(10.,11.,11.));
  EXPECT_THAT(NodeCoords[4], ElementsAre(8.,9.,11.));
  EXPECT_THAT(NodeCoords[5], ElementsAre(10.,10.,12.));
  EXPECT_THAT(NodeCoords[6], ElementsAre(9.,11.,12.));
  EXPECT_THAT(NodeCoords[7], ElementsAre(11.,12.,13.));

  EXPECT_NEAR(VolumeHexNonUniform(NodeCoords), 4., 1.e-12);

}
