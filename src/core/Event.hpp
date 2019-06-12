// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_EVENT_HPP_INCLUDED
#define OVK_CORE_EVENT_HPP_INCLUDED

#include <ovk/core/FloatingRef.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/IDMap.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <functional>
#include <memory>
#include <utility>

namespace ovk {

class event_listener_handle;

template <typename FuncSignature> class event;

template <typename... Args> class event<void(Args...)> {

public:

  event();

  template <typename F, OVK_FUNCDECL_REQUIRES(core::IsCallableWith<F, Args...>())>
    event_listener_handle AddListener(F Listener);

  void Trigger(Args... Arguments);

private:

  using listener = std::function<void(Args...)>;

  floating_ref_generator<event> FloatingRefGenerator_;

  // Set non-contiguous because std::function is not noexcept movable until C++20
  id_map<1,listener,false> Listeners_;

  friend class event_listener_handle;

};

class event_listener_handle {

public:

  event_listener_handle() = default;

  event_listener_handle(const event_listener_handle &Other) = delete;
  event_listener_handle(event_listener_handle &&Other) noexcept;

  event_listener_handle &operator=(event_listener_handle Other) noexcept;

  ~event_listener_handle() noexcept;

  void Reset();

private:

  template <typename FuncSignature> event_listener_handle(event<FuncSignature> &Event, int ID);

  int ID_ = -1;
  // std::function is not noexcept movable until C++20
  std::unique_ptr<std::function<void(int)>> RemoveFunc_;

  template <typename FuncSignature> friend class event;

};

}

#include <ovk/core/Event.inl>

#endif
