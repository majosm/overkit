// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Math.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>
#include <ovk/core/Elem.hpp>

#include <mpi.h>

using testing::ElementsAre;

class MathTests : public tests::mpi_test {};

TEST_F(MathTests, ColumnDeterminant2D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::ColumnDeterminant2D;

  ovk::elem<double,2> AI = {
    1.,
   -2.
  };
  ovk::elem<double,2> AJ = {
    2.,
    1.
  };
  EXPECT_EQ(ColumnDeterminant2D(AI, AJ), 5.);

}

TEST_F(MathTests, ColumnDeterminant3D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::ColumnDeterminant3D;

  ovk::elem<double,3> AI = {
    1.,
   -3.,
   -2.
  };
  ovk::elem<double,3> AJ = {
    2.,
    1.,
   -3.
  };
  ovk::elem<double,3> AK = {
    3.,
    2.,
    1.
  };
  EXPECT_EQ(ColumnDeterminant3D(AI, AJ, AK), 38.);

}

TEST_F(MathTests, ColumnSolve2D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::ColumnSolve2D;

  ovk::elem<double,2> AI = {
    1.,
   -2.
  };
  ovk::elem<double,2> AJ = {
    2.,
    1.
  };
  ovk::elem<double,2> B = {
    8.,
   -1
  };
  EXPECT_THAT(ColumnSolve2D(AI, AJ, B), ElementsAre(2.,3.));

}

TEST_F(MathTests, ColumnSolve3D) {

  if (TestComm().Rank() != 0) return;

  using ovk::core::ColumnSolve3D;

  ovk::elem<double,3> AI = {
    1.,
   -3.,
   -2.
  };
  ovk::elem<double,3> AJ = {
    2.,
    1.,
   -3.
  };
  ovk::elem<double,3> AK = {
    3.,
    2.,
    1.
  };
  ovk::elem<double,3> B = {
   20.,
    5.,
   -9.
  };
  EXPECT_THAT(ColumnSolve3D(AI, AJ, AK, B), ElementsAre(2.,3.,4.));

}

TEST_F(MathTests, LagrangeInterpLinear) {

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

TEST_F(MathTests, LagrangeInterpCubic) {

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
