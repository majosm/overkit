// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISPERSE_OVERWRITE_HPP_INCLUDED
#define OVK_CORE_DISPERSE_OVERWRITE_HPP_INCLUDED

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Context.hpp"
#include "ovk/core/DisperseBase.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Range.hpp"
#include "ovk/core/Tuple.hpp"

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {
namespace disperse_internal {

template <typename T, array_layout Layout> class disperse_overwrite : public disperse_base_for_type<
  T, Layout> {

protected:

  using parent_type = disperse_base_for_type<T, Layout>;

  using parent_type::Context_;
  using parent_type::Points_;
  using parent_type::Count_;
  using parent_type::FieldValuesRange_;
  using parent_type::FieldValuesIndexer_;
  using parent_type::PackedValues_;
  using parent_type::FieldValues_;

public:

  using value_type = T;

  disperse_overwrite(std::shared_ptr<context> &&Context, const array<int,2> &Points, int Count,
    const range &FieldValuesRange):
    parent_type(std::move(Context), Points, Count, FieldValuesRange)
  {}

  disperse_overwrite(const disperse_overwrite &Other) = delete;
  disperse_overwrite(disperse_overwrite &&Other) noexcept = default;

  disperse_overwrite &operator=(const disperse_overwrite &Other) = delete;
  disperse_overwrite &operator=(disperse_overwrite &&Other) noexcept = default;

  void Disperse(const void *PackedValuesVoid, void *FieldValuesVoid) {

    parent_type::SetBufferViews(PackedValuesVoid, FieldValuesVoid);

    long long NumPoints = Points_.Size(1);

    for (long long iPoint = 0; iPoint < NumPoints; ++iPoint) {
      tuple<int> Point = {
        Points_(0,iPoint),
        Points_(1,iPoint),
        Points_(2,iPoint)
      };
      long long iFieldValue = FieldValuesIndexer_.ToIndex(Point);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        FieldValues_(iCount)(iFieldValue) = PackedValues_(iCount)(iPoint);
      }
    }

  }

};

}}}

#endif
