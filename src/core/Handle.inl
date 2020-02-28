// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

template <typename T> template <typename F, OVK_FUNCDEF_REQUIRES(IsCallableWith<F, T>())>
  handle<T>::handle(T Handle, F Delete):
  Handle_(Handle),
  Delete_(new std::function<void(T &)>(std::move(Delete)))
{}

template <typename T> template <typename F, OVK_FUNCDEF_REQUIRES(!IsCallableWith<F, T>() &&
  IsCallableWith<F, T *>())> handle<T>::handle(T Handle, F Delete):
  Handle_(Handle),
  Delete_(new std::function<void(T &)>([Delete](T &Handle) { Delete(&Handle); }))
{}

template <typename T> handle<T>::~handle() noexcept {

  Reset();

}

template <typename T> handle<T>::operator bool() const {

  return static_cast<bool>(Handle_);

}

template <typename T> handle<T>::operator T() const {

  OVK_DEBUG_ASSERT(Handle_, "Invalid handle.");

  return *Handle_;

}

template <typename T> const T &handle<T>::Get() const {

  OVK_DEBUG_ASSERT(Handle_, "Invalid handle.");

  return *Handle_;

}

template <typename T> T &handle<T>::Get() {

  OVK_DEBUG_ASSERT(Handle_, "Invalid handle.");

  return *Handle_;

}

template <typename T> void handle<T>::Reset() {

  if (Handle_) {
    std::move(*Delete_)(*Handle_);
    Handle_.Reset();
    Delete_ = nullptr;
  }

}

template <typename T> T handle<T>::Release() {

  Delete_ = nullptr;

  return Handle_.Release();

}

template <typename T, typename F, OVK_FUNCDEF_REQUIRES(IsCallableWith<F, T>())> handle<T>
  MakeHandle(T Handle, F Delete) {
  return {std::move(Handle), std::move(Delete)};
}

template <typename T, typename F, OVK_FUNCDEF_REQUIRES(!IsCallableWith<F, T>() &&
  IsCallableWith<F, T *>())> handle<T> MakeHandle(T Handle, F Delete) {
  return {std::move(Handle), std::move(Delete)};
}

}}
