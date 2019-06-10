// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_RECV_HPP_INCLUDED
#define OVK_CORE_RECV_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/RecvMap.hpp>
#include <ovk/core/Request.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

class recv {

public:

  recv() = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, recv>::value)>
    recv(T &&Recv):
    Recv_(new model<T>(std::forward<T>(Recv)))
  {}

  recv(const recv &Other) = delete;
  recv(recv &&Other) noexcept = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, recv>::value)>
    recv &operator=(T &&Recv) {
    Recv_.reset(new model<T>(std::forward<T>(Recv)));
    return *this;
  }

  recv &operator=(const recv &Other) = delete;
  recv &operator=(recv &&Other) noexcept = default;

  request Recv(void **ReceiverValues) {
    return Recv_->Recv(ReceiverValues);
  }

private:

  class concept {
  public:
    virtual ~concept() noexcept {}
    virtual request Recv(void **ReceiverValues) = 0;
  };

  template <typename T> class model final : public concept {
  public:
    explicit model(T Recv):
      Recv_(std::move(Recv))
    {}
    virtual request Recv(void **ReceiverValues) override {
      return Recv_.Recv(ReceiverValues);
    }
  private:
    T Recv_;
  };

  std::unique_ptr<concept> Recv_;

};

recv CreateRecv(std::shared_ptr<context> Context, comm_view Comm, const recv_map &RecvMap, data_type
  ValueType, int Count, int Tag);

}}

#endif
