// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_REQUEST_HPP_INCLUDED
#define OVK_CORE_REQUEST_HPP_INCLUDED

#include <ovk/core/Global.hpp>

#include <mpi.h>

#include <memory>
#include <utility>
#include <vector>

namespace ovk {

class request {

public:

  request() = default;
  request(const request &) = delete;
  request(request &&) = default;

  request &operator=(const request &) = delete;
  request &operator=(request &&) = default;

  explicit operator bool() { return static_cast<bool>(Request_); }

  template <typename T> explicit request(T &&Request):
    Request_(new model<T>(std::forward<T>(Request)))
  {}

  template <typename T> request &operator=(T &&Request) {
    Request_.reset(new model<T>(std::forward<T>(Request)));
    return *this;
  }

  void Wait() {
    Request_->Wait();
  }

  static void core_WaitAll(int NumRequests, request **Requests);
  static void core_WaitAny(int NumRequests, request **Requests, int &Index);

private:

  class concept {
  public:
    virtual ~concept() {}
    virtual void Wait() = 0;
    virtual int NumMPIRequests() const = 0;
    virtual MPI_Request *MPIRequests() = 0;
  };

  template <typename T> class model : public concept {
  public:
    explicit model(T Request):
      Request_(std::move(Request))
    {}
    virtual void Wait() override {
      Request_.Wait();
    }
    virtual int NumMPIRequests() const override {
      return Request_.NumMPIRequests();
    }
    virtual MPI_Request *MPIRequests() override {
      return Request_.MPIRequests();
    }
  private:
    T Request_;
  };

  std::unique_ptr<concept> Request_;

  int NumMPIRequests() const {
    return Request_->NumMPIRequests();
  }

  MPI_Request *MPIRequests() {
    return Request_->MPIRequests();
  }

};

}

#endif
