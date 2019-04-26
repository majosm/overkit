// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
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

  ~request() {
    if (Request_) Request_->Wait();
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

  void Wait() {
    Request_->Wait();
    Request_.reset();
  }

  static void internal_WaitAll(array_view<request> Requests);
  static void internal_WaitAny(array_view<request> Requests, int &Index);
  // Needed for C API
  static void internal_WaitAll(array_view<request *> Requests);
  static void internal_WaitAny(array_view<request *> Requests, int &Index);

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual array_view<MPI_Request> MPIRequests() = 0;
    virtual void Finish(int iMPIRequest) = 0;
    virtual void Wait() = 0;
    virtual void StartProfileMemAlloc() const = 0;
    virtual void EndProfileMemAlloc() const = 0;
    virtual void StartProfileMPI() const = 0;
    virtual void EndProfileMPI() const = 0;
  };

  template <typename T> class model : public concept {
  public:
    explicit model(T Request):
      Request_(std::move(Request))
    {}
    virtual array_view<MPI_Request> MPIRequests() override {
      return Request_.MPIRequests();
    }
    virtual void Finish(int iMPIRequest) override {
      Request_.Finish(iMPIRequest);
    }
    virtual void Wait() override {
      Request_.Wait();
    }
    virtual void StartProfileMemAlloc() const override {
      Request_.StartProfileMemAlloc();
    }
    virtual void EndProfileMemAlloc() const override {
      Request_.EndProfileMemAlloc();
    }
    virtual void StartProfileMPI() const override {
      Request_.StartProfileMPI();
    }
    virtual void EndProfileMPI() const override {
      Request_.EndProfileMPI();
    }
  private:
    T Request_;
  };

  std::unique_ptr<concept> Request_;

  array_view<MPI_Request> MPIRequests() {
    return Request_->MPIRequests();
  }

  void Finish(int iMPIRequest) {
    Request_->Finish(iMPIRequest);
  }

  void StartProfileMemAlloc() {
    Request_->StartProfileMemAlloc();
  }

  void EndProfileMemAlloc() {
    Request_->EndProfileMemAlloc();
  }

  void StartProfileMPI() {
    Request_->StartProfileMPI();
  }

  void EndProfileMPI() {
    Request_->EndProfileMPI();
  }

};

void RequestWaitAll(array_view<request> Requests);
void RequestWaitAny(array_view<request> Requests, int &Index);
// Needed for C API
void RequestWaitAll(array_view<request *> Requests);
void RequestWaitAny(array_view<request *> Requests, int &Index);

}

#endif
