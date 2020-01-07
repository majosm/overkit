// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_HANDLE_HPP_INCLUDED
#define OVK_CORE_HANDLE_HPP_INCLUDED

#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/ScopeGuard.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <memory>
#include <utility>

namespace ovk {
namespace core {

template <typename T> class handle {

public:

  handle() = default;
  template <typename F, OVK_FUNCDECL_REQUIRES(IsCallableWith<F, T *>())> handle(T Handle, F Delete);
  template <typename F, OVK_FUNCDECL_REQUIRES(!IsCallableWith<F, T *>() && IsCallableWith<F, T>())>
    handle(T Handle, F Delete);

  explicit operator bool() const;

  operator T() const;

  T Get() const;

  void Reset();

private:

  std::shared_ptr<T> Ptr_;

};

template <typename T> bool operator==(const handle<T> &Left, const handle<T> &Right);
template <typename T> bool operator!=(const handle<T> &Left, const handle<T> &Right);

}}

#include <ovk/core/Handle.inl>

#endif
