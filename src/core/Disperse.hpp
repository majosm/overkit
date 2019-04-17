// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISPERSE_HPP_INCLUDED
#define OVK_CORE_DISPERSE_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Connectivity.hpp>
#include <ovk/core/ConnectivityR.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Exchange.hpp>
#include <ovk/core/Global.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

class disperse {

public:

  disperse() = default;

  template <typename T> disperse(T &&Disperse):
    Disperse_(new model<T>(std::forward<T>(Disperse)))
  {}

  disperse(const disperse &Other) = delete;
  disperse(disperse &&Other) noexcept = default;

  template <typename T> disperse &operator=(T &&Disperse) {
    Disperse_.reset(new model<T>(std::forward<T>(Disperse)));
    return *this;
  }

  disperse &operator=(const disperse &Other) = delete;
  disperse &operator=(disperse &&Other) noexcept = default;

  void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) {
    Disperse_->Initialize(Exchange, Count, GridValuesRange);
  }

  void Disperse(const void * const *ReceiverValues, void **GridValues) {
    Disperse_->Disperse(ReceiverValues, GridValues);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange) = 0;
    virtual void Disperse(const void * const *ReceiverValues, void **GridValues) = 0;
  };

  template <typename T> class model : public concept {
  public:
    using value_type = typename T::value_type;
    explicit model(T Disperse):
      Disperse_(std::move(Disperse))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, const range &GridValuesRange)
      override {
      Disperse_.Initialize(Exchange, Count, GridValuesRange);
    }
    virtual void Disperse(const void * const *ReceiverValues, void **GridValues) override {
      Disperse_.Disperse(ReceiverValues, GridValues);
    }
  private:
    T Disperse_;
  };

  std::unique_ptr<concept> Disperse_;

};

disperse MakeDisperse(disperse_op DisperseOp, data_type ValueType, array_layout Layout);

}}

#endif
