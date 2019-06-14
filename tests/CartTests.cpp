// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Cart.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

using testing::ElementsAre;

class CartTests : public tests::mpi_test {};

TEST_F(CartTests, PeriodicAdjust) {

  if (TestComm().Rank() != 0) return;

  auto CreateCart = [](int NumDims, bool IsPeriodic, bool Duplicated) -> ovk::cart {
    ovk::range Range = ovk::MakeEmptyRange(NumDims);
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      Range.Begin(iDim) = 1+iDim;
      Range.End(iDim) = 5+iDim;
    }
    ovk::tuple<bool> Periodic = ovk::MakeUniformTuple<bool>(NumDims, false);
    if (IsPeriodic) Periodic[NumDims-1] = true;
    ovk::periodic_storage PeriodicStorage = Duplicated ? ovk::periodic_storage::DUPLICATED :
      ovk::periodic_storage::UNIQUE;
    return {NumDims, Range, Periodic, PeriodicStorage};
  };

  // Non-periodic, unique, 2D
  {
    ovk::cart Cart = CreateCart(2, false, false);
    // Inside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({3,4,0});
      EXPECT_THAT(Tuple, ElementsAre(3,4,0));
    }
    // Outside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({5,6,0});
      EXPECT_THAT(Tuple, ElementsAre(5,6,0));
    }
  }

  // Non-periodic, unique, 3D
  {
    ovk::cart Cart = CreateCart(3, false, false);
    // Inside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({3,4,5});
      EXPECT_THAT(Tuple, ElementsAre(3,4,5));
    }
    // Outside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({5,6,7});
      EXPECT_THAT(Tuple, ElementsAre(5,6,7));
    }
  }

  // Non-periodic, duplicated, 2D
  {
    ovk::cart Cart = CreateCart(2, false, true);
    // Inside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({3,4,0});
      EXPECT_THAT(Tuple, ElementsAre(3,4,0));
    }
    // Outside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({5,6,0});
      EXPECT_THAT(Tuple, ElementsAre(5,6,0));
    }
  }

  // Non-periodic, duplicated, 3D
  {
    ovk::cart Cart = CreateCart(3, false, true);
    // Inside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({3,4,5});
      EXPECT_THAT(Tuple, ElementsAre(3,4,5));
    }
    // Outside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({5,6,7});
      EXPECT_THAT(Tuple, ElementsAre(5,6,7));
    }
  }

  // Periodic, unique, 2D
  {
    ovk::cart Cart = CreateCart(2, true, false);
    // Inside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({3,4,0});
      EXPECT_THAT(Tuple, ElementsAre(3,4,0));
    }
    // Outside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({5,6,0});
      EXPECT_THAT(Tuple, ElementsAre(5,2,0));
    }
  }

  // Periodic, unique, 3D
  {
    ovk::cart Cart = CreateCart(3, true, false);
    // Inside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({3,4,5});
      EXPECT_THAT(Tuple, ElementsAre(3,4,5));
    }
    // Outside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({5,6,7});
      EXPECT_THAT(Tuple, ElementsAre(5,6,3));
    }
  }

  // Periodic, duplicated, 2D
  {
    ovk::cart Cart = CreateCart(2, true, true);
    // Inside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({3,4,0});
      EXPECT_THAT(Tuple, ElementsAre(3,4,0));
    }
    // Outside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({5,6,0});
      EXPECT_THAT(Tuple, ElementsAre(5,3,0));
    }
  }

  // Periodic, duplicated, 3D
  {
    ovk::cart Cart = CreateCart(3, true, true);
    // Inside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({3,4,5});
      EXPECT_THAT(Tuple, ElementsAre(3,4,5));
    }
    // Outside
    {
      ovk::tuple<int> Tuple = Cart.PeriodicAdjust({5,6,7});
      EXPECT_THAT(Tuple, ElementsAre(5,6,4));
    }
  }

}

// TODO: Add remaining tests
