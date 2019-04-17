// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_COLLECT_HPP_INCLUDED
#define OVK_CORE_COLLECT_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Connectivity.hpp>
#include <ovk/core/ConnectivityD.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Exchange.hpp>
#include <ovk/core/Global.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

class collect {

public:

  collect() = default;

  template <typename T> collect(T &&Collect):
    Collect_(new model<T>(std::forward<T>(Collect)))
  {}

  collect(const collect &Other) = delete;
  collect(collect &&Other) noexcept = default;

  template <typename T> collect &operator=(T &&Collect) {
    Collect_.reset(new model<T>(std::forward<T>(Collect)));
    return *this;
  }

  collect &operator=(const collect &Other) = delete;
  collect &operator=(collect &&Other) noexcept = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {
    Collect_->Initialize(Exchange, Count, GridValuesRange);
  }

  void Collect(const void * const *GridValues, void **DonorValues) {
    Collect_->Collect(GridValues, DonorValues);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) = 0;
    virtual void Collect(const void * const *GridValues, void **DonorValues) = 0;
  };

  template <typename T> class model : public concept {
  public:
    using value_type = typename T::value_type;
    model(T Collect):
      Collect_(std::move(Collect))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange)
      override {
      Collect_.Initialize(Exchange, Count, GridValuesRange);
      GridValuesRange_ = GridValuesRange;
      const connectivity_d *Donors;
      GetConnectivityDonorSide(*Exchange.Connectivity_, Donors);
      GetConnectivityDonorSideCount(*Donors, NumDonors_);
      GridValues_.Resize({Count});
      DonorValues_.Resize({Count});
    }
    virtual void Collect(const void * const *GridValuesVoid, void **DonorValuesVoid) override {
      OVK_DEBUG_ASSERT(GridValuesVoid || GridValues_.Count() == 0, "Invalid grid values pointer.");
      OVK_DEBUG_ASSERT(DonorValuesVoid || DonorValues_.Count() == 0, "Invalid donor values pointer.");
      for (int iCount = 0; iCount < GridValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(GridValuesVoid[iCount] || NumDonors_ == 0, "Invalid grid values pointer.");
        GridValues_(iCount) = {static_cast<const value_type *>(GridValuesVoid[iCount]),
          {GridValuesRange_.Count()}};
      }
      for (int iCount = 0; iCount < DonorValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(DonorValuesVoid[iCount] || NumDonors_ == 0, "Invalid donor values "
          "pointer.");
        DonorValues_(iCount) = {static_cast<value_type *>(DonorValuesVoid[iCount]),
          {NumDonors_}};
      }
      Collect_.Collect(GridValues_, DonorValues_);
    }
  private:
    T Collect_;
    range GridValuesRange_;
    long long NumDonors_;
    array<array_view<const value_type>> GridValues_;
    array<array_view<value_type>> DonorValues_;
  };

  std::unique_ptr<concept> Collect_;

};

namespace collect_internal {
collect MakeCollectRow(collect_op CollectOp, data_type ValueType);
collect MakeCollectCol(collect_op CollectOp, data_type ValueType);
}

inline collect MakeCollect(collect_op CollectOp, data_type ValueType, array_layout Layout) {

  switch (Layout) {
  case array_layout::ROW_MAJOR:
    return collect_internal::MakeCollectRow(CollectOp, ValueType);
  case array_layout::COLUMN_MAJOR:
    return collect_internal::MakeCollectCol(CollectOp, ValueType);
  }

  return {};

}

}}

#endif
