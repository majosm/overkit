// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <ovk/core/DistributedField.hpp>

#include "tests/MPITest.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "support/Decomp.hpp"

#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/Partition.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <type_traits>
#include <utility>

using testing::ElementsAre;
using testing::ElementsAreArray;

class DistributedFieldTests : public tests::mpi_test {};

using support::CartesianDecomp;

namespace ovk {
namespace core {
template <typename T> class test_helper<distributed_field<T>> {
public:
  static const std::shared_ptr<const partition> &GetPartition(const distributed_field<T> &Field) {
    return Field.Partition_;
  }
  static const field<T> &GetValues(const distributed_field<T> &Field) { return Field.Values_; }
  static field<T> &GetValues(distributed_field<T> &Field) { return Field.Values_; }
};
}}

namespace {
// Have to also return cart comm at the moment because partition only stores comm_view
std::shared_ptr<const ovk::partition> CreatePartition(ovk::comm_view CommOfSize4, bool Duplicated,
  ovk::comm &CartComm) {

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(CommOfSize4)
  ));

  ovk::periodic_storage PeriodicStorage = Duplicated ? ovk::periodic_storage::DUPLICATED :
    ovk::periodic_storage::UNIQUE;

  ovk::cart Cart(2, {{1,2,0}, {13,14,1}}, {false,true,false}, PeriodicStorage);

  CartComm = ovk::CreateCartComm(CommOfSize4, Cart.Dimension(), {2,2,1}, Cart.Periodic());

  ovk::range LocalRange = CartesianDecomp(Cart.Dimension(), Cart.Range(), CartComm);
  ovk::range ExtendedRange = ovk::core::ExtendLocalRange(Cart, LocalRange, 1);

  ovk::core::decomp_hash DecompHash = ovk::core::CreateDecompHash(Cart.Dimension(), CartComm,
    LocalRange);

  ovk::array<int> NeighborRanks = ovk::core::DetectNeighbors(Cart, CartComm, LocalRange,
    DecompHash);

  return std::make_shared<ovk::partition>(std::move(Context), Cart, CartComm, LocalRange,
    ExtendedRange, 1, NeighborRanks);

}
}

TEST_F(DistributedFieldTests, Meta) {

  if (TestComm().Rank() != 0) return;

  using dfield = ovk::distributed_field<int>;

  EXPECT_TRUE((std::is_same<typename dfield::value_type, int>::value));
  EXPECT_TRUE((std::is_same<typename dfield::index_type, long long>::value));
  EXPECT_TRUE((std::is_same<typename dfield::tuple_element_type, int>::value));
  EXPECT_TRUE((std::is_same<typename dfield::tuple_type, ovk::tuple<int>>::value));
  EXPECT_TRUE((std::is_same<typename dfield::interval_type, ovk::range>::value));
  EXPECT_TRUE((std::is_same<typename dfield::indexer_type, ovk::range_indexer_c<long long>>::value));
  EXPECT_TRUE((std::is_same<typename dfield::iterator::pointer, int *>::value));
  EXPECT_TRUE((std::is_same<typename dfield::const_iterator::pointer, const int *>::value));

}

TEST_F(DistributedFieldTests, Create) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto SourcePartition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &GlobalRange = SourcePartition->GlobalRange();
    const ovk::range &ExtendedRange = SourcePartition->ExtendedRange();
    ovk::range_indexer_c<long long> GlobalIndexer(GlobalRange);

    // Default
    {
      ovk::distributed_field<int> Field;
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_FALSE(static_cast<bool>(Partition));
      EXPECT_TRUE(Values.Empty());
    }

    // Partition
    {
      ovk::distributed_field<int> Field(SourcePartition);
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
    }

    // Partition and value
    {
      ovk::distributed_field<int> Field(SourcePartition, 1);
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      ovk::field<int> ExpectedValues(ExtendedRange, 1);
      EXPECT_THAT(Values, ElementsAreArray(ExpectedValues));
    }

    // Partition and iterator
    {
      ovk::field<int> SourceValues(ExtendedRange);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            SourceValues(i,j,k) = GlobalIndexer.ToIndex(i,j,k);
          }
        }
      }
      ovk::distributed_field<int> Field(SourcePartition, SourceValues.Begin());
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      EXPECT_THAT(Values, ElementsAreArray(SourceValues));
    }

    // Partition and field view
    {
      ovk::field<int> SourceValues(ExtendedRange);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            SourceValues(i,j,k) = GlobalIndexer.ToIndex(i,j,k);
          }
        }
      }
      ovk::field_view<int> View(SourceValues);
      ovk::distributed_field<int> Field(SourcePartition, View);
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      EXPECT_THAT(Values, ElementsAreArray(SourceValues));
    }

    // Partition and field
    {
      ovk::field<int> SourceValues(ExtendedRange);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            SourceValues(i,j,k) = GlobalIndexer.ToIndex(i,j,k);
          }
        }
      }
      ovk::distributed_field<int> Field(SourcePartition, SourceValues);
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      EXPECT_THAT(Values, ElementsAreArray(SourceValues));
    }

  }

}

TEST_F(DistributedFieldTests, Copy) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto SourcePartition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &ExtendedRange = SourcePartition->ExtendedRange();

    // Copy construct
    {
      ovk::distributed_field<int> Field1(SourcePartition, CartComm.Rank());
      ovk::distributed_field<int> Field2 = Field1;
      auto &Partition = helper::GetPartition(Field2);
      auto &Values = helper::GetValues(Field2);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      ovk::field<int> ExpectedValues(ExtendedRange, CartComm.Rank());
      EXPECT_THAT(Values, ElementsAreArray(ExpectedValues));
    }

    // Copy assign, operator=
    {
      ovk::distributed_field<int> Field1(SourcePartition, CartComm.Rank());
      ovk::distributed_field<int> Field2;
      Field2 = Field1;
      auto &Partition = helper::GetPartition(Field2);
      auto &Values = helper::GetValues(Field2);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      ovk::field<int> ExpectedValues(ExtendedRange, CartComm.Rank());
      EXPECT_THAT(Values, ElementsAreArray(ExpectedValues));
    }

    // Copy assign, Assign
    {
      ovk::distributed_field<int> Field1(SourcePartition, CartComm.Rank());
      ovk::distributed_field<int> Field2;
      Field2.Assign(Field1);
      auto &Partition = helper::GetPartition(Field2);
      auto &Values = helper::GetValues(Field2);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      ovk::field<int> ExpectedValues(ExtendedRange, CartComm.Rank());
      EXPECT_THAT(Values, ElementsAreArray(ExpectedValues));
    }

  }

}

TEST_F(DistributedFieldTests, Move) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto SourcePartition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &ExtendedRange = SourcePartition->ExtendedRange();

    // Move construct
    {
      ovk::distributed_field<int> Field1(SourcePartition, CartComm.Rank());
      ovk::distributed_field<int> Field2 = std::move(Field1);
      auto &Partition1 = helper::GetPartition(Field1);
      auto &Values1 = helper::GetValues(Field1);
      auto &Partition2 = helper::GetPartition(Field2);
      auto &Values2 = helper::GetValues(Field2);
      EXPECT_EQ(Partition2, SourcePartition);
      EXPECT_THAT(Values2.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values2.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      ovk::field<int> ExpectedValues(ExtendedRange, CartComm.Rank());
      EXPECT_THAT(Values2, ElementsAreArray(ExpectedValues));
      EXPECT_FALSE(static_cast<bool>(Partition1));
      EXPECT_TRUE(Values1.Empty());
    }

    // Move assign, operator=
    {
      ovk::distributed_field<int> Field1(SourcePartition, CartComm.Rank());
      ovk::distributed_field<int> Field2;
      Field2 = std::move(Field1);
      auto &Partition1 = helper::GetPartition(Field1);
      auto &Values1 = helper::GetValues(Field1);
      auto &Partition2 = helper::GetPartition(Field2);
      auto &Values2 = helper::GetValues(Field2);
      EXPECT_EQ(Partition2, SourcePartition);
      EXPECT_THAT(Values2.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values2.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      ovk::field<int> ExpectedValues(ExtendedRange, CartComm.Rank());
      EXPECT_THAT(Values2, ElementsAreArray(ExpectedValues));
      EXPECT_FALSE(static_cast<bool>(Partition1));
      EXPECT_TRUE(Values1.Empty());
    }

    // Move assign, Assign
    {
      ovk::distributed_field<int> Field1(SourcePartition, CartComm.Rank());
      ovk::distributed_field<int> Field2;
      Field2.Assign(std::move(Field1));
      auto &Partition1 = helper::GetPartition(Field1);
      auto &Values1 = helper::GetValues(Field1);
      auto &Partition2 = helper::GetPartition(Field2);
      auto &Values2 = helper::GetValues(Field2);
      EXPECT_EQ(Partition2, SourcePartition);
      EXPECT_THAT(Values2.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values2.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      ovk::field<int> ExpectedValues(ExtendedRange, CartComm.Rank());
      EXPECT_THAT(Values2, ElementsAreArray(ExpectedValues));
      EXPECT_FALSE(static_cast<bool>(Partition1));
      EXPECT_TRUE(Values1.Empty());
    }

  }

}

TEST_F(DistributedFieldTests, Assign) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto SourcePartition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &GlobalRange = SourcePartition->GlobalRange();
    const ovk::range &ExtendedRange = SourcePartition->ExtendedRange();
    ovk::range_indexer_c<long long> GlobalIndexer(GlobalRange);

    // Partition
    {
      ovk::distributed_field<int> Field;
      Field.Assign(SourcePartition);
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
    }

    // Partition and value
    {
      ovk::distributed_field<int> Field;
      Field.Assign(SourcePartition, 1);
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      ovk::field<int> ExpectedValues(ExtendedRange, 1);
      EXPECT_THAT(Values, ElementsAreArray(ExpectedValues));
    }

    // Partition and iterator
    {
      ovk::field<int> SourceValues(ExtendedRange);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            SourceValues(i,j,k) = GlobalIndexer.ToIndex(i,j,k);
          }
        }
      }
      ovk::distributed_field<int> Field;
      Field.Assign(SourcePartition, SourceValues.Begin());
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      EXPECT_THAT(Values, ElementsAreArray(SourceValues));
    }

    // Partition and field view
    {
      ovk::field<int> SourceValues(ExtendedRange);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            SourceValues(i,j,k) = GlobalIndexer.ToIndex(i,j,k);
          }
        }
      }
      ovk::field_view<int> View(SourceValues);
      ovk::distributed_field<int> Field;
      Field.Assign(SourcePartition, View);
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      EXPECT_THAT(Values, ElementsAreArray(SourceValues));
    }

    // Partition and field
    {
      ovk::field<int> SourceValues(ExtendedRange);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            SourceValues(i,j,k) = GlobalIndexer.ToIndex(i,j,k);
          }
        }
      }
      ovk::distributed_field<int> Field;
      Field.Assign(SourcePartition, SourceValues);
      auto &Partition = helper::GetPartition(Field);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Partition, SourcePartition);
      EXPECT_THAT(Values.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
      EXPECT_THAT(Values.Extents().End(), ElementsAreArray(ExtendedRange.End()));
      EXPECT_THAT(Values, ElementsAreArray(SourceValues));
    }

  }

}

TEST_F(DistributedFieldTests, Exchange) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &LocalRange = Partition->LocalRange();

    ovk::distributed_field<int> Field(Partition, -1);

    auto &Values = helper::GetValues(Field);

    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          Values(i,j,k) = CartComm.Rank();
        }
      }
    }

    // Sanity check
    int MinValue = CartComm.Size();
    for (auto &Value : Values) {
      MinValue = ovk::Min(MinValue, Value);
    }
    EXPECT_EQ(MinValue, -1);

    Field.Exchange();

    // Simple check since halo is tested elsewhere
    MinValue = CartComm.Size();
    for (auto &Value : Values) {
      MinValue = ovk::Min(MinValue, Value);
    }
    EXPECT_GE(MinValue, 0);

  }

}

TEST_F(DistributedFieldTests, Fill) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    // Constant value
    {
      ovk::comm CartComm;
      auto Partition = CreatePartition(CommOfSize4, false, CartComm);
      const ovk::range &ExtendedRange = Partition->ExtendedRange();
      ovk::distributed_field<int> Field(Partition, 0);
      Field.Fill(1);
      auto &Values = helper::GetValues(Field);
      ovk::field<int> ExpectedValues(ExtendedRange, 1);
      EXPECT_THAT(Values, ElementsAreArray(ExpectedValues));
    }

    // Range and constant value, interior
    {
      ovk::comm CartComm;
      auto Partition = CreatePartition(CommOfSize4, false, CartComm);
      const ovk::range &GlobalRange = Partition->GlobalRange();
      const ovk::range &ExtendedRange = Partition->ExtendedRange();
      ovk::distributed_field<int> Field(Partition, 0);
      ovk::range FillRange = ovk::MakeEmptyRange(2);
      for (int iDim = 0; iDim < 2; ++iDim) {
        FillRange.Begin(iDim) = GlobalRange.Begin(iDim)+2;
        FillRange.End(iDim) = GlobalRange.End(iDim)-2;
      }
      Field.Fill(FillRange, 1);
      auto &Values = helper::GetValues(Field);
      ovk::range LocalFillRange = ovk::IntersectRanges(ExtendedRange, FillRange);
      ovk::field<int> ExpectedValues(ExtendedRange, 0);
      ExpectedValues.Fill(LocalFillRange, 1);
      EXPECT_THAT(Values, ElementsAreArray(ExpectedValues));
    }

    // Range and constant value, periodic boundary, unique
    {
      ovk::comm CartComm;
      auto Partition = CreatePartition(CommOfSize4, false, CartComm);
      const ovk::range &GlobalRange = Partition->GlobalRange();
      const ovk::range &ExtendedRange = Partition->ExtendedRange();
      ovk::distributed_field<int> Field(Partition, 0);
      ovk::range FillRange = {
        {GlobalRange.Begin(0),GlobalRange.End(1)-2,0},
        {GlobalRange.End(0),GlobalRange.End(1)+2,1}
      };
      Field.Fill(FillRange, 1);
      auto &Values = helper::GetValues(Field);
      ovk::field<int> ExpectedValues(ExtendedRange, 0);
      if (CartComm.Rank() == 0 || CartComm.Rank() == 2) {
        ExpectedValues.Fill({{ExtendedRange.Begin(0),ExtendedRange.Begin(1),0},
          {ExtendedRange.End(0),GlobalRange.Begin(1)+2,1}}, 1);
      } else {
        ExpectedValues.Fill({{ExtendedRange.Begin(0),GlobalRange.End(1)-2,0},
          {ExtendedRange.End(0),ExtendedRange.End(1),1}}, 1);
      }
      EXPECT_THAT(Values, ElementsAreArray(ExpectedValues));
    }

    // Range and constant value, periodic boundary, duplicated
    {
      ovk::comm CartComm;
      auto Partition = CreatePartition(CommOfSize4, true, CartComm);
      const ovk::range &GlobalRange = Partition->GlobalRange();
      const ovk::range &ExtendedRange = Partition->ExtendedRange();
      ovk::distributed_field<int> Field(Partition, 0);
      ovk::range FillRange = {
        {GlobalRange.Begin(0),GlobalRange.End(1)-2,0},
        {GlobalRange.End(0),GlobalRange.End(1)+2,1}
      };
      Field.Fill(FillRange, 1);
      auto &Values = helper::GetValues(Field);
      ovk::field<int> ExpectedValues(ExtendedRange, 0);
      if (CartComm.Rank() == 0 || CartComm.Rank() == 2) {
        ExpectedValues.Fill({{ExtendedRange.Begin(0),ExtendedRange.Begin(1),0},
          {ExtendedRange.End(0),GlobalRange.Begin(1)+3,1}}, 1);
      } else {
        ExpectedValues.Fill({{ExtendedRange.Begin(0),GlobalRange.End(1)-2,0},
          {ExtendedRange.End(0),ExtendedRange.End(1),1}}, 1);
      }
      EXPECT_THAT(Values, ElementsAreArray(ExpectedValues));
    }

    // Iterator
    {
      ovk::comm CartComm;
      auto Partition = CreatePartition(CommOfSize4, false, CartComm);
      const ovk::range &GlobalRange = Partition->GlobalRange();
      const ovk::range &ExtendedRange = Partition->ExtendedRange();
      ovk::range_indexer_c<long long> GlobalIndexer(GlobalRange);
      ovk::distributed_field<int> Field(Partition, 0);
      ovk::field<int> SourceValues(ExtendedRange);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            SourceValues(i,j,k) = GlobalIndexer.ToIndex(i,j,k);
          }
        }
      }
      Field.Fill(SourceValues.Begin());
      auto &Values = helper::GetValues(Field);
      EXPECT_THAT(Values, ElementsAreArray(SourceValues));
    }

    // Field view
    {
      ovk::comm CartComm;
      auto Partition = CreatePartition(CommOfSize4, false, CartComm);
      const ovk::range &GlobalRange = Partition->GlobalRange();
      const ovk::range &ExtendedRange = Partition->ExtendedRange();
      ovk::range_indexer_c<long long> GlobalIndexer(GlobalRange);
      ovk::distributed_field<int> Field(Partition, 0);
      ovk::field<int> SourceValues(ExtendedRange);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            SourceValues(i,j,k) = GlobalIndexer.ToIndex(i,j,k);
          }
        }
      }
      ovk::field_view<int> View(SourceValues);
      Field.Fill(View);
      auto &Values = helper::GetValues(Field);
      EXPECT_THAT(Values, ElementsAreArray(SourceValues));
    }

    // Field
    {
      ovk::comm CartComm;
      auto Partition = CreatePartition(CommOfSize4, false, CartComm);
      const ovk::range &GlobalRange = Partition->GlobalRange();
      const ovk::range &ExtendedRange = Partition->ExtendedRange();
      ovk::range_indexer_c<long long> GlobalIndexer(GlobalRange);
      ovk::distributed_field<int> Field(Partition, 0);
      ovk::field<int> SourceValues(ExtendedRange);
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            SourceValues(i,j,k) = GlobalIndexer.ToIndex(i,j,k);
          }
        }
      }
      Field.Fill(SourceValues);
      auto &Values = helper::GetValues(Field);
      EXPECT_THAT(Values, ElementsAreArray(SourceValues));
    }

  }

}

TEST_F(DistributedFieldTests, ConvertToBool) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);

    // Doesn't have partition
    {
      ovk::distributed_field<int> Field;
      EXPECT_FALSE(static_cast<bool>(Field));
    }

    // Has partition
    {
      ovk::distributed_field<int> Field(Partition);
      EXPECT_TRUE(static_cast<bool>(Field));
    }

  }

}

TEST_F(DistributedFieldTests, Partition) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);

    ovk::distributed_field<int> Field(Partition);

    EXPECT_EQ(&Field.Partition(), Partition.get());
    EXPECT_EQ(Field.SharedPartition(), Partition);

  }

}

TEST_F(DistributedFieldTests, Cart) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);

    ovk::distributed_field<int> Field(Partition);
    EXPECT_EQ(&Field.Cart(), &Partition->Cart());

  }

}

TEST_F(DistributedFieldTests, Comm) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);

    ovk::distributed_field<int> Field(Partition);

    EXPECT_EQ(Field.Comm().Get(), Partition->Comm().Get());

  }

}

TEST_F(DistributedFieldTests, Ranges) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);

    ovk::distributed_field<int> Field(Partition);

    EXPECT_EQ(&Field.GlobalRange(), &Partition->GlobalRange());
    EXPECT_EQ(&Field.LocalRange(), &Partition->LocalRange());
    EXPECT_EQ(&Field.ExtendedRange(), &Partition->ExtendedRange());

  }

}

TEST_F(DistributedFieldTests, Extents) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &ExtendedRange = Partition->ExtendedRange();

    ovk::distributed_field<int> Field(Partition);

    EXPECT_THAT(Field.Extents().Begin(), ElementsAreArray(ExtendedRange.Begin()));
    EXPECT_THAT(Field.Extents().End(), ElementsAreArray(ExtendedRange.End()));

  }

}

TEST_F(DistributedFieldTests, Size) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &ExtendedRange = Partition->ExtendedRange();

    ovk::distributed_field<int> Field(Partition);

    EXPECT_THAT(Field.Size(), ElementsAreArray(ExtendedRange.Size()));
    EXPECT_EQ(Field.Size(0), ExtendedRange.Size(0));
    EXPECT_EQ(Field.Size(1), ExtendedRange.Size(1));
    EXPECT_EQ(Field.Size(2), ExtendedRange.Size(2));

  }

}

TEST_F(DistributedFieldTests, Count) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &ExtendedRange = Partition->ExtendedRange();

    ovk::distributed_field<int> Field(Partition);

    EXPECT_EQ(Field.Count(), ExtendedRange.Count());

  }

}

// TEST_F(DistributedFieldTests, Indexer) {

//   ASSERT_GE(TestComm().Size(), 4);

//   ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

//   if (CommOfSize4) {

//     ovk::comm CartComm;
//     auto Partition = CreatePartition(CommOfSize4, false, CartComm);
//     const ovk::range &ExtendedRange = Partition->ExtendedRange();

//     ovk::distributed_field<int> Field(Partition);

//     ovk::range_indexer_c<long long> ExpectedIndexer(ExtendedRange);

//     EXPECT_EQ(Field.Indexer(), ExpectedIndexer);

//   }

// }

TEST_F(DistributedFieldTests, Values) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);

    // Const
    {
      const ovk::distributed_field<int> Field(Partition, -1);
      auto &Values = helper::GetValues(Field);
      EXPECT_TRUE((std::is_same<decltype(Field.Values()), const ovk::field<int> &>::value));
      EXPECT_EQ(&Field.Values(), &Values);
    }

    // Non-const
    {
      ovk::distributed_field<int> Field(Partition, -1);
      auto &Values = helper::GetValues(Field);
      EXPECT_TRUE((std::is_same<decltype(Field.Values()), ovk::field<int> &>::value));
      EXPECT_EQ(&Field.Values(), &Values);
    }

  }

}

TEST_F(DistributedFieldTests, ParenthesisOperator) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &LocalRange = Partition->LocalRange();

    // Tuple
    {
      ovk::distributed_field<int> Field(Partition);
      auto &Values = helper::GetValues(Field);
      ovk::tuple<int> Point = LocalRange.Begin();
      EXPECT_EQ(&Field(Point), Values.Data(Point));
    }

    // Separate elements
    {
      ovk::distributed_field<int> Field(Partition);
      auto &Values = helper::GetValues(Field);
      int i = LocalRange.Begin(0);
      int j = LocalRange.Begin(1);
      int k = LocalRange.Begin(2);
      EXPECT_EQ(&Field(i,j,k), Values.Data(i,j,k));
    }

  }

}

TEST_F(DistributedFieldTests, BracketOperator) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &LocalRange = Partition->LocalRange();
    const ovk::range &ExtendedRange = Partition->ExtendedRange();

    ovk::range_indexer_c<long long> Indexer(ExtendedRange);

    ovk::distributed_field<int> Field(Partition);

    auto &Values = helper::GetValues(Field);
    long long Index = Indexer.ToIndex(LocalRange.Begin());

    EXPECT_EQ(&Field[Index], Values.Data()+Index);

  }

}

TEST_F(DistributedFieldTests, Data) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &LocalRange = Partition->LocalRange();

    // No argument
    {
      ovk::distributed_field<int> Field(Partition);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Field.Data(), Values.Data());
    }

    // Tuple
    {
      ovk::distributed_field<int> Field(Partition);
      auto &Values = helper::GetValues(Field);
      ovk::tuple<int> Point = LocalRange.Begin();
      EXPECT_EQ(Field.Data(Point), Values.Data(Point));
    }

    // Separate elements
    {
      ovk::distributed_field<int> Field(Partition);
      auto &Values = helper::GetValues(Field);
      int i = LocalRange.Begin(0);
      int j = LocalRange.Begin(1);
      int k = LocalRange.Begin(2);
      EXPECT_EQ(Field.Data(i,j,k), Values.Data(i,j,k));
    }

  }

}

TEST_F(DistributedFieldTests, BeginEnd) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using helper = ovk::core::test_helper<ovk::distributed_field<int>>;

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &ExtendedRange = Partition->ExtendedRange();

    // Const
    {
      const ovk::distributed_field<int> Field(Partition, 1);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Field.Begin().Pointer(), Values.Data());
      EXPECT_EQ(Field.End().Pointer(), Values.Data()+ExtendedRange.Count());
      long long Sum = 0;
      for (auto &Value : Field) Sum += Value;
      EXPECT_EQ(Sum, ExtendedRange.Count());
    }

    // Non-const
    {
      ovk::distributed_field<int> Field(Partition, 0);
      auto &Values = helper::GetValues(Field);
      EXPECT_EQ(Field.Begin().Pointer(), Values.Data());
      EXPECT_EQ(Field.End().Pointer(), Values.Data()+ExtendedRange.Count());
      for (auto &Value : Field) Value = 1;
      long long Sum = 0;
      for (auto &Value : Values) Sum += Value;
      EXPECT_EQ(Sum, ExtendedRange.Count());
    }

  }

}

TEST_F(DistributedFieldTests, Traits) {

  ASSERT_GE(TestComm().Size(), 4);

  ovk::comm CommOfSize4 = CreateSubsetComm(TestComm(), TestComm().Rank() < 4);

  if (CommOfSize4) {

    using dfield = ovk::distributed_field<int>;

    ovk::comm CartComm;
    auto Partition = CreatePartition(CommOfSize4, false, CartComm);
    const ovk::range &ExtendedRange = Partition->ExtendedRange();

    EXPECT_TRUE(ovk::core::IsArray<dfield>());
    EXPECT_TRUE((std::is_same<ovk::core::array_value_type<dfield>, int>::value));
    EXPECT_EQ(ovk::core::ArrayRank<dfield>(), OVK_MAX_DIMS);
    EXPECT_EQ(ovk::core::ArrayLayout<dfield>(), ovk::array_layout::COLUMN_MAJOR);
    EXPECT_TRUE(ovk::core::ArrayHasRuntimeExtents<dfield>());

    dfield Field(Partition);
    EXPECT_THAT(ovk::core::ArrayExtents(Field).Begin(), ElementsAreArray(ExtendedRange.Begin()));
    EXPECT_THAT(ovk::core::ArrayExtents(Field).End(), ElementsAreArray(ExtendedRange.End()));
    EXPECT_EQ(ovk::core::ArrayData(Field), Field.Data());

  }

}
