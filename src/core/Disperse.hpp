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
      const connectivity_r *Receivers;
      GetConnectivityReceiverSide(*Exchange.Connectivity_, Receivers);
      GetConnectivityReceiverSideCount(*Receivers, NumReceivers_);
      GridValuesRange_ = GridValuesRange;
      ReceiverValues_.Resize({Count});
      GridValues_.Resize({Count});
    }
    virtual void Disperse(const void * const *ReceiverValuesVoid, void **GridValuesVoid) override {
      OVK_DEBUG_ASSERT(ReceiverValuesVoid || ReceiverValues_.Count() == 0, "Invalid receiver values "
        "pointer.");
      OVK_DEBUG_ASSERT(GridValuesVoid || GridValues_.Count() == 0, "Invalid grid values pointer.");
      for (int iCount = 0; iCount < ReceiverValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(ReceiverValuesVoid[iCount] || NumReceivers_ == 0, "Invalid receiver "
          "values pointer.");
        ReceiverValues_(iCount) = {static_cast<const value_type *>(ReceiverValuesVoid[iCount]),
          {NumReceivers_}};
      }
      for (int iCount = 0; iCount < GridValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(GridValuesVoid[iCount] || NumReceivers_ == 0, "Invalid grid values "
          "pointer.");
        GridValues_(iCount) = {static_cast<value_type *>(GridValuesVoid[iCount]),
          {GridValuesRange_.Count()}};
      }
      Disperse_.Disperse(ReceiverValues_, GridValues_);
    }
  private:
    T Disperse_;
    long long NumReceivers_;
    range GridValuesRange_;
    array<array_view<const value_type>> ReceiverValues_;
    array<array_view<value_type>> GridValues_;
  };

  std::unique_ptr<concept> Disperse_;

};

disperse MakeDisperse(disperse_op DisperseOp, data_type ValueType, array_layout Layout);

}}

#endif
