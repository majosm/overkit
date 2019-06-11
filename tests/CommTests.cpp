// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Comm.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Tuple.hpp>

#include <mpi.h>

using testing::ElementsAre;

namespace ovk {
namespace core {
template <> class test_helper<comm_view> {
public:
  static MPI_Comm GetMPIComm(const comm_view &Comm) { return Comm.Comm_; }
};
template <> class test_helper<comm> {
public:
  static MPI_Comm GetMPIComm(const comm &Comm) { return Comm.View_.Get(); }
};
}}

TEST(CommViewTests, Create) {

  using helper = ovk::core::test_helper<ovk::comm_view>;

  // Default
  {
    ovk::comm_view Comm;
    EXPECT_EQ(helper::GetMPIComm(Comm), MPI_COMM_NULL);
  }

  // Null
  {
    ovk::comm_view Comm(MPI_COMM_NULL);
    EXPECT_EQ(helper::GetMPIComm(Comm), MPI_COMM_NULL);
  }

  // Not null
  {
    ovk::comm_view Comm(MPI_COMM_WORLD);
    EXPECT_EQ(helper::GetMPIComm(Comm), MPI_COMM_WORLD);
  }

}

TEST(CommViewTests, Copy) {

  using helper = ovk::core::test_helper<ovk::comm_view>;

  // Copy construct
  {
    ovk::comm_view Comm1(MPI_COMM_WORLD);
    ovk::comm_view Comm2(Comm1);
    EXPECT_EQ(helper::GetMPIComm(Comm2), MPI_COMM_WORLD);
  }

  // Copy assign
  {
    ovk::comm_view Comm1(MPI_COMM_WORLD);
    ovk::comm_view Comm2;
    Comm2 = Comm1;
    EXPECT_EQ(helper::GetMPIComm(Comm2), MPI_COMM_WORLD);
  }

}

TEST(CommViewTests, Reset) {

  using helper = ovk::core::test_helper<ovk::comm_view>;

  // Initially null
  {
    ovk::comm_view Comm;
    Comm.Reset();
    EXPECT_EQ(helper::GetMPIComm(Comm), MPI_COMM_NULL);
  }

  // Initially not null
  {
    ovk::comm_view Comm(MPI_COMM_WORLD);
    Comm.Reset();
    EXPECT_EQ(helper::GetMPIComm(Comm), MPI_COMM_NULL);
  }

}

TEST(CommViewTests, Equality) {

  int WorldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

  MPI_Comm SplitComm;
  MPI_Comm_split(MPI_COMM_WORLD, WorldRank % 2, WorldRank, &SplitComm);

  ovk::comm_view Comm1(MPI_COMM_WORLD);
  ovk::comm_view Comm2(MPI_COMM_WORLD);
  ovk::comm_view Comm3(SplitComm);

  // Self
  EXPECT_EQ(Comm1 == Comm1, true);
  EXPECT_EQ(Comm1 != Comm1, false);

  // Other with same comm
  EXPECT_EQ(Comm1 == Comm2, true);
  EXPECT_EQ(Comm1 != Comm2, false);

  // Different comms
  EXPECT_EQ(Comm1 == Comm3, false);
  EXPECT_EQ(Comm1 != Comm3, true);

}

TEST(CommViewTests, ConvertToBool) {

  // Null
  {
    ovk::comm_view Comm;
    EXPECT_EQ(static_cast<bool>(Comm), false);
  }

  // Not null
  {
    ovk::comm_view Comm(MPI_COMM_WORLD);
    EXPECT_EQ(static_cast<bool>(Comm), true);
  }

}

TEST(CommViewTests, Get) {

  ovk::comm_view Comm(MPI_COMM_WORLD);
  EXPECT_EQ(Comm.Get(), MPI_COMM_WORLD);

}

TEST(CommViewTests, ConvertToMPIComm) {

  auto ImplicitConvertFunc = [](MPI_Comm Comm) -> MPI_Comm { return Comm; };

  ovk::comm_view Comm(MPI_COMM_WORLD);
  EXPECT_EQ(ImplicitConvertFunc(Comm), MPI_COMM_WORLD);

}

TEST(CommViewTests, Size) {

  int WorldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &WorldSize);

  ovk::comm_view Comm(MPI_COMM_WORLD);
  EXPECT_EQ(Comm.Size(), WorldSize);

}

TEST(CommViewTests, Rank) {

  int WorldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

  ovk::comm_view Comm(MPI_COMM_WORLD);
  EXPECT_EQ(Comm.Rank(), WorldRank);

}

TEST(CommTests, CreateDestroy) {

  using helper = ovk::core::test_helper<ovk::comm>;

//   // Duplicate MPI_COMM_WORLD so we can modify error handler
//   MPI_Comm ErrorComm;
//   MPI_Comm_dup(MPI_COMM_WORLD, &ErrorComm);
//   MPI_Comm_set_errhandler(ErrorComm, MPI_ERRORS_RETURN);

  // Null
  {
    ovk::comm Comm(MPI_COMM_NULL);
    EXPECT_EQ(helper::GetMPIComm(Comm), MPI_COMM_NULL);
  }

  // Not null
  {
    MPI_Comm CommRaw;
    MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);
    {
      ovk::comm Comm(CommRaw);
      EXPECT_EQ(helper::GetMPIComm(Comm), CommRaw);
    }
//     // Make sure communicator was cleaned up
//     int Compare;
//     int Error = MPI_Comm_compare(ErrorComm, CommRaw, &Compare);
//     EXPECT_NE(Error, MPI_SUCCESS);
  }

//   MPI_Comm_free(&ErrorComm);

}

TEST(CommTests, Move) {

  using helper = ovk::core::test_helper<ovk::comm>;

//   // Duplicate MPI_COMM_WORLD so we can modify error handler
//   MPI_Comm ErrorComm;
//   MPI_Comm_dup(MPI_COMM_WORLD, &ErrorComm);
//   MPI_Comm_set_errhandler(ErrorComm, MPI_ERRORS_RETURN);

  // Move construction
  {
    MPI_Comm CommRaw;
    MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);
    ovk::comm Comm1(CommRaw);
    ovk::comm Comm2(std::move(Comm1));
    EXPECT_EQ(helper::GetMPIComm(Comm2), CommRaw);
    EXPECT_EQ(helper::GetMPIComm(Comm1), MPI_COMM_NULL);
  }

  // Move assignment
  {
    MPI_Comm Comm1Raw;
    MPI_Comm_dup(MPI_COMM_WORLD, &Comm1Raw);
    int WorldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);
    MPI_Comm Comm2Raw;
    MPI_Comm_split(MPI_COMM_WORLD, WorldRank % 2, WorldRank, &Comm2Raw);
    ovk::comm Comm1(Comm1Raw);
    ovk::comm Comm2(Comm2Raw);
    Comm2 = std::move(Comm1);
    EXPECT_EQ(helper::GetMPIComm(Comm2), Comm1Raw);
    EXPECT_EQ(helper::GetMPIComm(Comm1), MPI_COMM_NULL);
//     // Make sure original communicator in Comm2 was cleaned up
//     int Error = MPI_Comm_compare(ErrorComm, Comm2Raw, &Compare);
//     EXPECT_NE(Error, MPI_SUCCESS);
  }

//   MPI_Comm_free(&ErrorComm);

}

TEST(CommTests, Reset) {

  using helper = ovk::core::test_helper<ovk::comm>;

//   // Duplicate MPI_COMM_WORLD so we can modify error handler
//   MPI_Comm ErrorComm;
//   MPI_Comm_dup(MPI_COMM_WORLD, &ErrorComm);
//   MPI_Comm_set_errhandler(ErrorComm, MPI_ERRORS_RETURN);

  {
    MPI_Comm CommRaw;
    MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);
    ovk::comm Comm(CommRaw);
//     MPI_Comm CommRaw = helper::GetMPIComm(Comm);
    Comm.Reset();
    EXPECT_EQ(helper::GetMPIComm(Comm), MPI_COMM_NULL);
//     // Make sure communicator was cleaned up
//     int Compare;
//     int Error = MPI_Comm_compare(ErrorComm, CommRaw, &Compare);
//     EXPECT_NE(Error, MPI_SUCCESS);
  }

//   MPI_Comm_free(&ErrorComm);

}

TEST(CommTests, ConvertToBool) {

  // Null
  {
    ovk::comm Comm;
    EXPECT_EQ(static_cast<bool>(Comm), false);
  }

  // Not null
  {
    MPI_Comm CommRaw;
    MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);
    ovk::comm Comm(CommRaw);
    EXPECT_EQ(static_cast<bool>(Comm), true);
  }

}

TEST(CommTests, Get) {

  MPI_Comm CommRaw;
  MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);

  ovk::comm Comm(CommRaw);
  EXPECT_EQ(Comm.Get(), CommRaw);

}

TEST(CommTests, ConvertToMPIComm) {

  auto ImplicitConvertFunc = [](MPI_Comm Comm) -> MPI_Comm { return Comm; };

  MPI_Comm CommRaw;
  MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);

  ovk::comm Comm(CommRaw);
  EXPECT_EQ(ImplicitConvertFunc(Comm), CommRaw);

}

TEST(CommTests, Size) {

  int WorldSize;
  MPI_Comm_size(MPI_COMM_WORLD, &WorldSize);

  MPI_Comm CommRaw;
  MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);

  ovk::comm Comm(CommRaw);
  EXPECT_EQ(Comm.Size(), WorldSize);

}

TEST(CommTests, Rank) {

  int WorldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

  MPI_Comm CommRaw;
  MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);

  ovk::comm Comm(CommRaw);
  EXPECT_EQ(Comm.Rank(), WorldRank);

}

TEST(CommTests, Duplicate) {

  MPI_Comm CommRaw;
  MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);

  ovk::comm Comm(CommRaw);
  ovk::comm DuplicatedComm = ovk::DuplicateComm(Comm);

  int Compare;
  MPI_Comm_compare(Comm, DuplicatedComm, &Compare);

  EXPECT_EQ(Compare, MPI_CONGRUENT);

}

TEST(CommTests, Subset) {

  MPI_Comm CommRaw;
  MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);

  ovk::comm Comm(CommRaw);
  ASSERT_GE(Comm.Size(), 2);

  bool InSubset = Comm.Rank() % 2 == 0;
  ovk::comm SubsetComm = ovk::CreateSubsetComm(Comm, InSubset);

  if (InSubset) {
    EXPECT_TRUE(static_cast<bool>(SubsetComm));
  } else {
    EXPECT_FALSE(static_cast<bool>(SubsetComm));
  }

}

TEST(CommTests, Cartesian) {

  MPI_Comm CommRaw;
  MPI_Comm_dup(MPI_COMM_WORLD, &CommRaw);

  ovk::comm Comm(CommRaw);
  ASSERT_GE(Comm.Size(), 6);

  ovk::comm SubsetComm = ovk::CreateSubsetComm(Comm, Comm.Rank() < 6);

  if (SubsetComm) {

    ovk::tuple<int> Dims = {2,3,1};
    ovk::tuple<bool> Periodic = {false,true,false};
    ovk::comm CartComm = ovk::CreateCartComm(SubsetComm, 2, Dims, Periodic);

    EXPECT_TRUE(ovk::IsCartComm(CartComm));
    EXPECT_EQ(ovk::GetCartCommDimension(CartComm), 2);
    EXPECT_THAT(ovk::GetCartCommDims(CartComm), ElementsAre(2,3,1));
    EXPECT_THAT(ovk::GetCartCommPeriodic(CartComm), ElementsAre(false,true,false));
    switch (CartComm.Rank()) {
    // Lower corner
    case 0:
      EXPECT_THAT(ovk::GetCartCommCoords(CartComm), ElementsAre(0,0,0));
      break;
    // Middle
    case 1:
      EXPECT_THAT(ovk::GetCartCommCoords(CartComm), ElementsAre(0,1,0));
      break;
    // Upper corner
    case 5:
      EXPECT_THAT(ovk::GetCartCommCoords(CartComm), ElementsAre(1,2,0));
      break;
    default:
      break;
    }

  }

}
