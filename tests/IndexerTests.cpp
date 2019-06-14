// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/Indexer.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <ovk/core/Comm.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Interval.hpp>

#include <mpi.h>

#include <type_traits>

using testing::ElementsAre;

class IndexerTests : public tests::mpi_test {};

namespace ovk {
namespace core {
template <array_layout Layout> class test_helper<indexer<long long, int, 3, Layout>> {
public:
  using indexer_type = indexer<long long, int, 3, Layout>;
  static const elem<int,3> &GetBegin(const indexer_type &Indexer) { return Indexer.Begin_; }
  static const elem<long long,3> &GetStride(const indexer_type &Indexer) { return Indexer.Stride_; }
};
}}

TEST_F(IndexerTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using indexer = ovk::indexer<long long, int, 3>;

  EXPECT_TRUE((std::is_same<typename indexer::index_type, long long>::value));
  EXPECT_TRUE((std::is_same<typename indexer::tuple_element_type, int>::value));
  EXPECT_EQ(int(indexer::Rank), 3);
  EXPECT_EQ(ovk::array_layout(indexer::Layout), ovk::array_layout::ROW_MAJOR);
  EXPECT_TRUE((std::is_same<typename indexer::tuple_type, ovk::elem<int,3>>::value));
  EXPECT_TRUE((std::is_same<typename indexer::stride_type, ovk::elem<long long,3>>::value));
  EXPECT_TRUE((std::is_same<typename indexer::interval_type, ovk::interval<int,3>>::value));

  using indexer_col = ovk::indexer<long long, int, 3, ovk::array_layout::COLUMN_MAJOR>;

  EXPECT_EQ(ovk::array_layout(indexer_col::Layout), ovk::array_layout::COLUMN_MAJOR);

}

TEST_F(IndexerTests, Create) {

  if (TestComm().Rank() != 0) return;

  using indexer_row = ovk::indexer<long long, int, 3, ovk::array_layout::ROW_MAJOR>;
  using indexer_col = ovk::indexer<long long, int, 3, ovk::array_layout::COLUMN_MAJOR>;
  using helper_row = ovk::core::test_helper<indexer_row>;
  using helper_col = ovk::core::test_helper<indexer_col>;

  // Row major, default
  {
    indexer_row Indexer;
    EXPECT_THAT(helper_row::GetBegin(Indexer), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetStride(Indexer), ElementsAre(0,0,1));
  }

  // Row major, size
  {
    indexer_row Indexer({1,2,3});
    EXPECT_THAT(helper_row::GetBegin(Indexer), ElementsAre(0,0,0));
    EXPECT_THAT(helper_row::GetStride(Indexer), ElementsAre(6,3,1));
  }

  // Row major, interval
  {
    indexer_row Indexer({{1,2,3}, {4,5,6}});
    EXPECT_THAT(helper_row::GetBegin(Indexer), ElementsAre(1,2,3));
    EXPECT_THAT(helper_row::GetStride(Indexer), ElementsAre(9,3,1));
  }

  // Column major, default
  {
    indexer_col Indexer;
    EXPECT_THAT(helper_col::GetBegin(Indexer), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetStride(Indexer), ElementsAre(1,0,0));
  }

  // Column major, size
  {
    indexer_col Indexer({1,2,3});
    EXPECT_THAT(helper_col::GetBegin(Indexer), ElementsAre(0,0,0));
    EXPECT_THAT(helper_col::GetStride(Indexer), ElementsAre(1,1,2));
  }

  // Column major, interval
  {
    indexer_col Indexer({{1,2,3}, {4,5,6}});
    EXPECT_THAT(helper_col::GetBegin(Indexer), ElementsAre(1,2,3));
    EXPECT_THAT(helper_col::GetStride(Indexer), ElementsAre(1,3,9));
  }

}

TEST_F(IndexerTests, Equality) {

  if (TestComm().Rank() != 0) return;

  using indexer = ovk::indexer<long long, int, 3>;

  // Self
  {
    indexer Indexer({{1,2,3}, {4,5,6}});
    EXPECT_TRUE(Indexer == Indexer);
    EXPECT_FALSE(Indexer != Indexer);
  }

  // Other with same values
  {
    indexer Indexer1({{1,2,3}, {4,5,6}});
    indexer Indexer2({{1,2,3}, {4,5,6}});
    EXPECT_TRUE(Indexer1 == Indexer2);
    EXPECT_FALSE(Indexer1 != Indexer2);
  }

  // Different begin
  {
    indexer Indexer1({{1,2,3}, {4,5,6}});
    indexer Indexer2({{1,3,2}, {4,6,5}});
    EXPECT_FALSE(Indexer1 == Indexer2);
    EXPECT_TRUE(Indexer1 != Indexer2);
  }

  // Different stride
  {
    indexer Indexer1({{1,2,3}, {4,5,6}});
    indexer Indexer2({{1,2,3}, {4,6,5}});
    EXPECT_FALSE(Indexer1 == Indexer2);
    EXPECT_TRUE(Indexer1 != Indexer2);
  }

}

TEST_F(IndexerTests, Begin) {

  if (TestComm().Rank() != 0) return;

  using indexer = ovk::indexer<long long, int, 3>;

  indexer Indexer({{1,2,3}, {4,5,6}});
  EXPECT_THAT(Indexer.Begin(), ElementsAre(1,2,3));

}

TEST_F(IndexerTests, Stride) {

  if (TestComm().Rank() != 0) return;

  using indexer = ovk::indexer<long long, int, 3>;

  indexer Indexer({{1,2,3}, {4,5,6}});
  EXPECT_THAT(Indexer.Stride(), ElementsAre(9,3,1));

}

TEST_F(IndexerTests, TupleToIndex) {

  if (TestComm().Rank() != 0) return;

  using indexer_row = ovk::indexer<long long, int, 3, ovk::array_layout::ROW_MAJOR>;
  using indexer_col = ovk::indexer<long long, int, 3, ovk::array_layout::COLUMN_MAJOR>;

  // Row major, iterator
  {
    indexer_row Indexer({{1,1,1}, {4,5,6}});
    ovk::elem<int,3> Tuple = {2,3,4};
    long long Index = Indexer.ToIndex(Tuple.Data());
    EXPECT_EQ(Index, 33);
  }

  // Row major, array
  {
    indexer_row Indexer({{1,1,1}, {4,5,6}});
    long long Index = Indexer.ToIndex(ovk::elem<int,3>(2,3,4));
    EXPECT_EQ(Index, 33);
  }

  // Row major, separate
  {
    indexer_row Indexer({{1,1,1}, {4,5,6}});
    long long Index = Indexer.ToIndex(2,3,4);
    EXPECT_EQ(Index, 33);
  }

  // Column major, iterator
  {
    indexer_col Indexer({{1,1,1}, {4,5,6}});
    ovk::elem<int,3> Tuple = {2,3,4};
    long long Index = Indexer.ToIndex(Tuple.Data());
    EXPECT_EQ(Index, 43);
  }

  // Column major, array
  {
    indexer_col Indexer({{1,1,1}, {4,5,6}});
    long long Index = Indexer.ToIndex(ovk::elem<int,3>(2,3,4));
    EXPECT_EQ(Index, 43);
  }

  // Column major, separate
  {
    indexer_col Indexer({{1,1,1}, {4,5,6}});
    long long Index = Indexer.ToIndex(2,3,4);
    EXPECT_EQ(Index, 43);
  }

}

TEST_F(IndexerTests, IndexToTuple) {

  if (TestComm().Rank() != 0) return;

  using indexer_row = ovk::indexer<long long, int, 3, ovk::array_layout::ROW_MAJOR>;
  using indexer_col = ovk::indexer<long long, int, 3, ovk::array_layout::COLUMN_MAJOR>;

  // Row major
  {
    indexer_row Indexer({{1,1,1}, {4,5,6}});
    ovk::elem<int,3> Tuple = Indexer.ToTuple(33);
    EXPECT_THAT(Tuple, ElementsAre(2,3,4));
  }

  // Column major
  {
    indexer_col Indexer({{1,1,1}, {4,5,6}});
    ovk::elem<int,3> Tuple = Indexer.ToTuple(43);
    EXPECT_THAT(Tuple, ElementsAre(2,3,4));
  }

}
