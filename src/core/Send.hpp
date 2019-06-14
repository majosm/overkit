// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_SEND_HPP_INCLUDED
#define OVK_CORE_SEND_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/Context.hpp>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Request.hpp>
#include <ovk/core/SendMap.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <mpi.h>

#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

class send {

public:

  send() = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, send>::value)>
    send(T &&Send):
    Send_(new model<T>(std::forward<T>(Send)))
  {}

  send(const send &Other) = delete;
  send(send &&Other) noexcept = default;

  template <typename T, OVK_FUNCTION_REQUIRES(!std::is_same<remove_cvref<T>, send>::value)>
    send &operator=(T &&Send) {
    Send_.reset(new model<T>(std::forward<T>(Send)));
    return *this;
  }

  send &operator=(const send &Other) = delete;
  send &operator=(send &&Other) noexcept = default;

  request Send(const void *Values) {
    return Send_->Send(Values);
  }

private:

  class concept {
  public:
    virtual ~concept() noexcept {}
    virtual request Send(const void *Values) = 0;
  };

  template <typename T> class model final : public concept {
  public:
    explicit model(T Send):
      Send_(std::move(Send))
    {}
    virtual request Send(const void *Values) override {
      return Send_.Send(Values);
    }
  private:
    T Send_;
  };

  std::unique_ptr<concept> Send_;

};

send CreateSend(std::shared_ptr<context> Context, comm_view Comm, const send_map &SendMap, data_type
  ValueType, int Count, int Tag);

}}

#endif
