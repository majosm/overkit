namespace ovk {
namespace core {

template <typename T> template <typename F, OVK_FUNCDEF_REQUIRES(IsCallableWith<F, T *>())>
  handle<T>::handle(T Handle, F Delete) {

  if (Handle == handle_traits<T>::NullValue) return;

  auto CleanUpHandle = core::OnScopeExit([&] {
    Delete(&Handle);
  });

  T *HandlePtr = new T(Handle);

  CleanUpHandle.Dismiss();

  Ptr_ = std::shared_ptr<T>(HandlePtr, [Delete](T *HandlePtr) {
    Delete(HandlePtr);
    delete HandlePtr;
  });

}

template <typename T> handle<T>::operator bool() const {

  return Ptr_ != nullptr;

}

template <typename T> handle<T>::operator T() const {

  return Ptr_ ? *Ptr_ : handle_traits<T>::NullValue;

}

template <typename T> T handle<T>::Get() const {

  return Ptr_ ? *Ptr_ : handle_traits<T>::NullValue;

}

template <typename T> void handle<T>::Reset() {

  Ptr_.reset();

}

template <typename T> bool operator==(const handle<T> &Left, const handle<T> &Right) {

  return Left.Get() == Right.Get();

}

template <typename T> bool operator!=(const handle<T> &Left, const handle<T> &Right) {

  return Left.Get() != Right.Get();

}

}}
