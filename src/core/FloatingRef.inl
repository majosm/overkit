// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

inline floating_ref_generator::floating_ref_generator():
  ReferenceLoc_(new void *(this))
{}

inline floating_ref_generator::floating_ref_generator(const floating_ref_generator &Other):
  ReferenceLoc_(new void *(this))
{}

inline floating_ref_generator::floating_ref_generator(floating_ref_generator &&Other) noexcept:
  ReferenceLoc_(std::move(Other.ReferenceLoc_))
{
  *ReferenceLoc_ = this;
}

inline floating_ref_generator &floating_ref_generator::operator=(const floating_ref_generator
  &Other) {

  return *this;

}

inline floating_ref_generator &floating_ref_generator::operator=(floating_ref_generator &&Other)
  noexcept {

  ReferenceLoc_ = std::move(Other.ReferenceLoc_);

  *ReferenceLoc_ = this;

  return *this;

}

template <typename T> floating_ref<typename std::remove_reference<T>::type>
  floating_ref_generator::Generate(T &&Target) const {

  return {ReferenceLoc_.get(), std::forward<T>(Target)};

}

template <typename T> floating_ref<T>::floating_ref(void * const *ReferenceLoc, T &Target):
  ReferenceLoc_(ReferenceLoc),
  Offset_(reinterpret_cast<byte_ptr>(&Target) - reinterpret_cast<byte_ptr>(*ReferenceLoc_))
{}

template <typename T> template <typename U, OVK_FUNCDEF_REQUIRES(!std::is_same<U, T>::value &&
  std::is_convertible<U *, T *>::value)> floating_ref<T>::floating_ref(const floating_ref<U>
  &Other):
  ReferenceLoc_(Other.ReferenceLoc_)
{

  using other_byte_ptr = typename floating_ref<U>::byte_ptr;

  auto OtherTarget = reinterpret_cast<U *>(reinterpret_cast<other_byte_ptr>(*ReferenceLoc_) +
    Other.Offset_);
  auto Target = static_cast<T *>(OtherTarget);

  Offset_ = Other.Offset_ + (reinterpret_cast<byte_ptr>(Target) -
    reinterpret_cast<other_byte_ptr>(OtherTarget));

}

template <typename T> floating_ref<T>::floating_ref(floating_ref &&Other) noexcept:
  ReferenceLoc_(Other.ReferenceLoc_),
  Offset_(Other.Offset_)
{
  Other.ReferenceLoc_ = nullptr;
}

template <typename T> floating_ref<T> &floating_ref<T>::operator=(floating_ref &&Other) noexcept {

  ReferenceLoc_ = Other.ReferenceLoc_;
  Offset_ = Other.Offset_;

  Other.ReferenceLoc_ = nullptr;

  return *this;

}

template <typename T> T &floating_ref<T>::operator*() const {

  return *reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(*ReferenceLoc_) + Offset_);

}

template <typename T> T *floating_ref<T>::operator->() const {

  return reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(*ReferenceLoc_) + Offset_);

}

template <typename T> T &floating_ref<T>::Get() const {

  return *reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(*ReferenceLoc_) + Offset_);

}

template <typename T> template <typename U, OVK_FUNCDEF_REQUIRES(!std::is_same<U, T>::value &&
  !std::is_convertible<T *, U *>::value && std::is_base_of<T, U>::value)> floating_ref<T>::operator
  floating_ref<U>() const {

  using other_byte_ptr = typename floating_ref<U>::byte_ptr;

  auto Target = reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(*ReferenceLoc_) + Offset_);
  auto OtherTarget = static_cast<U *>(Target);

  floating_ref<U> Other;
  Other.ReferenceLoc_ = ReferenceLoc_;
  Other.Offset_ = Offset_ + (reinterpret_cast<other_byte_ptr>(OtherTarget) -
    reinterpret_cast<byte_ptr>(Target));

  return Other;

}

template <typename U, typename T, OVK_FUNCDEF_REQUIRES(std::is_convertible<T *, U *>::value ||
  std::is_base_of<T, U>::value)> floating_ref<U> FloatingRefCast(const floating_ref<T> &FloatingRef)
  {

  return static_cast<floating_ref<U>>(FloatingRef);

}

}
