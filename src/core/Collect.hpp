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
    }
    virtual void Collect(const void * const *GridValues, void **DonorValues) override {
      Collect_.Collect(GridValues, DonorValues);
    }
  private:
    T Collect_;
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
