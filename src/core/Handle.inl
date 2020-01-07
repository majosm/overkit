namespace ovk {
namespace core {

template <typename T> template <typename F, OVK_FUNCDEF_REQUIRES(IsCallableWith<F, T *>())>
  handle<T>::handle(T Handle, F Delete) {

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

  OVK_DEBUG_ASSERT(Ptr_, "Invalid handle.");

  return *Ptr_;

}

template <typename T> T handle<T>::Get() const {

  OVK_DEBUG_ASSERT(Ptr_, "Invalid handle.");

  return *Ptr_;

}

template <typename T> void handle<T>::Reset() {

  Ptr_.reset();

}

template <typename T> bool operator==(const handle<T> &Left, const handle<T> &Right) {

  return (!Left && !Right) || ((Left && Right) && Left.Get() == Right.Get());

}

template <typename T> bool operator!=(const handle<T> &Left, const handle<T> &Right) {

  return !(Left == Right);

}

}}
