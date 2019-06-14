// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_SCOPE_GUARD_HPP_LOADED
#define OVK_CORE_SCOPE_GUARD_HPP_LOADED

#include <ovk/core/Global.hpp>

#include <utility>

namespace ovk {
namespace core {

template <typename F> class scope_guard {

public:

  scope_guard(F ScopeExitFunc);

  ~scope_guard() noexcept;

  scope_guard(const scope_guard &) = delete;
  scope_guard(scope_guard &&Other) noexcept;

  scope_guard &operator=(const scope_guard &Other) = delete;
  scope_guard &operator=(scope_guard &&Other) noexcept;

  scope_guard &Dismiss();

  bool Dismissed() const;

private:

  F ScopeExitFunc_;
  bool Dismissed_;

};

template <typename F> scope_guard<typename std::decay<F>::type> OnScopeExit(F &&ScopeExitFunc) {
  return {std::forward<F>(ScopeExitFunc)};
}

}}

#include <ovk/core/ScopeGuard.inl>

#endif
