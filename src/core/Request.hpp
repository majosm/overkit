// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_REQUEST_HPP_INCLUDED
#define OVK_CORE_REQUEST_HPP_INCLUDED

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {

class request {

public:

  request() = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<core::remove_cvref<T>, request>::value)>
    request(T &&Request):
    Request_(new model<T>(std::forward<T>(Request)))
  {}

  ~request() noexcept {
    Wait();
  }

  request(const request &Other) = delete;
  request(request &&Other) noexcept = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<core::remove_cvref<T>, request>::value)>
    request &operator=(T &&Request) {
    Request_.reset(new model<T>(std::forward<T>(Request)));
    return *this;
  }

  request &operator=(const request &Other) = delete;
  request &operator=(request &&Other) noexcept = default;

  explicit operator bool() { return static_cast<bool>(Request_); }

  void Wait();

  static void internal_WaitAll(array_view<request> Requests);
  static void internal_WaitAny(array_view<request> Requests, int &Index);
  // Needed for C API
  static void internal_WaitAll(array_view<request *> Requests);
  static void internal_WaitAny(array_view<request *> Requests, int &Index);

private:

  class concept {
  public:
    virtual ~concept() noexcept {}
    virtual array_view<MPI_Request> MPIRequests() = 0;
    virtual void OnMPIRequestComplete(int iMPIRequest) = 0;
    virtual void OnComplete() = 0;
    virtual void StartWaitTime() const = 0;
    virtual void StopWaitTime() const = 0;
    virtual void StartMPITime() const = 0;
    virtual void StopMPITime() const = 0;
  };

  template <typename T> class model final : public concept {
  public:
    explicit model(T Request):
      Request_(std::move(Request))
    {}
    virtual array_view<MPI_Request> MPIRequests() override {
      return Request_.MPIRequests();
    }
    virtual void OnMPIRequestComplete(int iMPIRequest) override {
      Request_.OnMPIRequestComplete(iMPIRequest);
    }
    virtual void OnComplete() override {
      Request_.OnComplete();
    }
    virtual void StartWaitTime() const override {
      Request_.StartWaitTime();
    }
    virtual void StopWaitTime() const override {
      Request_.StopWaitTime();
    }
    virtual void StartMPITime() const override {
      Request_.StartMPITime();
    }
    virtual void StopMPITime() const override {
      Request_.StopMPITime();
    }
  private:
    T Request_;
  };

  std::unique_ptr<concept> Request_;

  array_view<MPI_Request> MPIRequests_() {
    return Request_->MPIRequests();
  }

  void OnMPIRequestComplete_(int iMPIRequest) {
    Request_->OnMPIRequestComplete(iMPIRequest);
  }

  void OnComplete_() {
    Request_->OnComplete();
  }

  void StartWaitTime_() {
    Request_->StartWaitTime();
  }

  void StopWaitTime_() {
    Request_->StopWaitTime();
  }

  void StartMPITime_() {
    Request_->StartMPITime();
  }

  void StopMPITime_() {
    Request_->StopMPITime();
  }

  void Reset_() {
    Request_.reset();
  }

};

void WaitAll(array_view<request> Requests);
void WaitAny(array_view<request> Requests, int &Index);
// Needed for C API
void WaitAll(array_view<request *> Requests);
void WaitAny(array_view<request *> Requests, int &Index);

}

#endif
