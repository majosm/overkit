// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_HANDLE_HPP_INCLUDED
#define OVK_CORE_HANDLE_HPP_INCLUDED

#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <functional>
#include <utility>

namespace ovk {
namespace core {

template <typename T> class handle {

public:

  handle() = default;
  template <typename F, OVK_FUNCDECL_REQUIRES(IsCallableWith<F, T>())> handle(T Handle, F Delete);
  template <typename F, OVK_FUNCDECL_REQUIRES(!IsCallableWith<F, T>() && IsCallableWith<F, T *>())>
    handle(T Handle, F Delete);

  handle(const handle &Other) = delete;
  handle(handle &&Other) noexcept = default;

  handle &operator=(const handle &Other) = delete;
  handle &operator=(handle &&Other) noexcept = default;

  ~handle() noexcept;

  explicit operator bool() const;

  operator T() const;

  const T &Get() const;
  T &Get();

  void Reset();

  T Release();

private:

  optional<T> Handle_;
  // std::function is not noexcept movable until C++20
  std::unique_ptr<std::function<void(T &)>> Delete_;

};

template <typename T, typename F, OVK_FUNCDECL_REQUIRES(IsCallableWith<F, T>())> handle<T>
  MakeHandle(T Handle, F Delete);
template <typename T, typename F, OVK_FUNCDECL_REQUIRES(!IsCallableWith<F, T>() &&
  IsCallableWith<F, T *>())> handle<T> MakeHandle(T Handle, F Delete);

}}

#include <ovk/core/Handle.inl>

#endif
