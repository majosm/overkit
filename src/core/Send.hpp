// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_SEND_HPP_INCLUDED
#define OVK_CORE_SEND_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Connectivity.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Exchange.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Profiler.hpp>
#include <ovk/core/Request.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

class send {

public:

  send() = default;

  template <typename T> send(T &&Send):
    Send_(new model<T>(std::forward<T>(Send)))
  {}

  send(const send &Other) = delete;
  send(send &&Other) noexcept = default;

  template <typename T> send &operator=(T &&Send) {
    Send_.reset(new model<T>(std::forward<T>(Send)));
    return *this;
  }

  send &operator=(const send &Other) = delete;
  send &operator=(send &&Other) noexcept = default;

  void Initialize(const exchange &Exchange, int Count, int Tag) {
    Send_->Initialize(Exchange, Count, Tag);
  }

  request Send(const void * const *DonorValues) {
    return Send_->Send(DonorValues);
  }

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Initialize(const exchange &Exchange, int Count, int Tag) = 0;
    virtual request Send(const void * const *DonorValues) = 0;
  };

  template <typename T> class model : public concept {
  public:
    using value_type = typename T::value_type;
    explicit model(T Send):
      Send_(std::move(Send))
    {}
    virtual void Initialize(const exchange &Exchange, int Count, int Tag) override {
      Send_.Initialize(Exchange, Count, Tag);
      const connectivity_d *Donors;
      GetConnectivityDonorSide(*Exchange.Connectivity_, Donors);
      GetConnectivityDonorSideCount(*Donors, NumDonors_);
      DonorValues_.Resize({Count});
    }
    virtual request Send(const void * const *DonorValuesVoid) override {
      OVK_DEBUG_ASSERT(DonorValuesVoid || DonorValues_.Count() == 0, "Invalid donor values "
        "pointer.");
      for (int iCount = 0; iCount < DonorValues_.Count(); ++iCount) {
        OVK_DEBUG_ASSERT(DonorValuesVoid[iCount] || NumDonors_ == 0, "Invalid donor values "
          "pointer.");
        DonorValues_(iCount) = {static_cast<const value_type *>(DonorValuesVoid[iCount]),
          {NumDonors_}};
      }
      return Send_.Send(DonorValues_);
    }
  private:
    T Send_;
    long long NumDonors_;
    array<array_view<const value_type>> DonorValues_;
  };

  std::unique_ptr<concept> Send_;

};

send MakeSend(data_type ValueType);

}}

#endif
