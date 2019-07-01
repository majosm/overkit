// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename F> scope_guard<F>::scope_guard(F ScopeExitFunc):
  ScopeExitFunc_(std::move(ScopeExitFunc)),
  Dismissed_(false)
{}

template <typename F> scope_guard<F>::scope_guard(scope_guard &&Other) noexcept:
  ScopeExitFunc_(std::move(Other.ScopeExitFunc_)),
  Dismissed_(Other.Dismissed_)
{
  Other.Dismissed_ = true;
}

template <typename F> scope_guard<F>::~scope_guard() noexcept {

  if (!Dismissed_) {
    ScopeExitFunc_();
  }

}

template <typename F> scope_guard<F> &scope_guard<F>::operator=(scope_guard &&Other) noexcept {

  if (!Dismissed_) {
    ScopeExitFunc_();
  }

  ScopeExitFunc_ = std::move(Other.ScopeExitFunc_);
  Dismissed_ = Other.Dismissed_;
  Other.Dismissed_ = true;

  return *this;

}

template <typename F> scope_guard<F> &scope_guard<F>::Dismiss() {

  Dismissed_ = true;

  return *this;

}

template <typename F> bool scope_guard<F>::Dismissed() const {

  return Dismissed_;

}

}}
