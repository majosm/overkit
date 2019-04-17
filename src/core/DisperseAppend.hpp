// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISPERSE_APPEND_HPP_INCLUDED
#define OVK_CORE_DISPERSE_APPEND_HPP_INCLUDED

#include "ovk/core/Array.hpp"
#include "ovk/core/ArrayView.hpp"
#include "ovk/core/DisperseBase.hpp"
#include "ovk/core/Comm.hpp"
#include "ovk/core/ConnectivityR.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Elem.hpp"
#include "ovk/core/Exchange.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Profiler.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

namespace ovk {
namespace core {
namespace disperse_internal {

template <typename T, array_layout Layout> class disperse_append : public disperse_base_for_type<T,
  Layout> {

protected:

  using parent_type = disperse_base_for_type<T, Layout>;

  using parent_type::Receivers_;
  using parent_type::Grid_;
  using parent_type::Count_;
  using parent_type::NumReceivers_;
  using parent_type::GridValuesRange_;
  using parent_type::GridValuesIndexer_;
  using parent_type::Points_;
  using parent_type::ReceiverValues_;
  using parent_type::GridValues_;

public:

  using value_type = T;

  disperse_append() = default;
  disperse_append(const disperse_append &Other) = delete;
  disperse_append(disperse_append &&Other) noexcept = default;

  disperse_append &operator=(const disperse_append &Other) = delete;
  disperse_append &operator=(disperse_append &&Other) noexcept = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {

    parent_type::Initialize(Exchange, Count, GridValuesRange);

  }

  void Disperse(const void * const *ReceiverValuesVoid, void **GridValuesVoid) {

    parent_type::SetBufferViews(ReceiverValuesVoid, GridValuesVoid);

    for (long long iReceiver = 0; iReceiver < NumReceivers_; ++iReceiver) {
      elem<int,MAX_DIMS> Point = {
        Points_(0,iReceiver),
        Points_(1,iReceiver),
        Points_(2,iReceiver)
      };
      long long iPoint = GridValuesIndexer_.ToIndex(Point);
      for (int iCount = 0; iCount < Count_; ++iCount) {
        GridValues_(iCount)(iPoint) += ReceiverValues_(iCount)(iReceiver);
      }
    }

  }

};

}}}

#endif
