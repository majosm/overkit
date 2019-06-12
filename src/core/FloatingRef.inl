// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

template <typename T> floating_ref_generator<T>::floating_ref_generator(T &Target):
  Resource_(new resource(&Target))
{}

template <typename T> floating_ref_generator<T>::floating_ref_generator(const floating_ref_generator
  &Other) {

  if (Other.Resource_) {

    auto OtherTarget = static_cast<T *>(Other.Resource_->Target);

    std::ptrdiff_t Offset = reinterpret_cast<byte_ptr>(OtherTarget) -
      reinterpret_cast<byte_ptr>(&Other);

    auto Target = reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(this) + Offset);

    Resource_.reset(new resource(Target));

  } else {

    Resource_.reset();

  }

}

template <typename T> floating_ref_generator<T>::floating_ref_generator(floating_ref_generator
  &&Other) noexcept {

  if (Other.Resource_) {

    auto OtherTarget = static_cast<T *>(Other.Resource_->Target);

    std::ptrdiff_t Offset = reinterpret_cast<byte_ptr>(OtherTarget) -
      reinterpret_cast<byte_ptr>(&Other);

    auto Target = reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(this) + Offset);

    Resource_ = std::move(Other.Resource_);
    Resource_->Target = Target;

  } else {

    Resource_.reset();

  }

}

template <typename T> floating_ref_generator<T> &floating_ref_generator<T>::operator=(const
  floating_ref_generator &Other) {

  if (Other.Resource_) {

    auto OtherTarget = static_cast<T *>(Other.Resource_->Target);

    std::ptrdiff_t Offset = reinterpret_cast<byte_ptr>(OtherTarget) -
      reinterpret_cast<byte_ptr>(&Other);

    auto Target = reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(this) + Offset);

    Resource_.reset(new resource(Target));

  } else {

    Resource_.reset();

  }

  return *this;

}

template <typename T> floating_ref_generator<T> &floating_ref_generator<T>::operator=(
  floating_ref_generator &&Other) noexcept {

  if (Other.Resource_) {

    auto OtherTarget = static_cast<T *>(Other.Resource_->Target);

    std::ptrdiff_t Offset = reinterpret_cast<byte_ptr>(OtherTarget) -
      reinterpret_cast<byte_ptr>(&Other);

    auto Target = reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(this) + Offset);

    Resource_ = std::move(Other.Resource_);
    Resource_->Target = Target;

  } else {

    Resource_.reset();

  }

  return *this;

}

template <typename T> template <typename U, OVK_FUNCDEF_REQUIRES(std::is_convertible<T *,
  U *>::value || std::is_base_of<T, U>::value)> floating_ref<U> floating_ref_generator<T>::
  Generate() const {

  return static_cast<floating_ref<U>>(floating_ref<T>(*Resource_));

}

template <typename T> template <typename U, OVK_FUNCDEF_REQUIRES(std::is_convertible<U *,
  T *>::value)> floating_ref<T>::floating_ref(const floating_ref<U> &Other):
  Resource_(Other.Resource_)
{

  using other_byte_ptr = typename floating_ref<U>::byte_ptr;

  auto OtherTarget = reinterpret_cast<U *>(reinterpret_cast<other_byte_ptr>(Resource_->Target) +
    Other.TypeOffset_);
  auto Target = static_cast<T *>(OtherTarget);

  TypeOffset_ = Other.TypeOffset_ + (reinterpret_cast<byte_ptr>(Target) -
    reinterpret_cast<other_byte_ptr>(OtherTarget));

}

template <typename T> floating_ref<T>::floating_ref(floating_ref &&Other) noexcept:
  Resource_(Other.Resource_),
  TypeOffset_(Other.TypeOffset_)
{
  Other.Resource_ = nullptr;
  Other.TypeOffset_ = 0;
}

template <typename T> floating_ref<T> &floating_ref<T>::operator=(floating_ref &&Other) noexcept {

  Resource_ = Other.Resource_;
  TypeOffset_ = Other.TypeOffset_;

  Other.Resource_ = nullptr;
  Other.TypeOffset_ = 0;

  return *this;

}

template <typename T> T &floating_ref<T>::operator*() const {

  return *reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(Resource_->Target) + TypeOffset_);

}

template <typename T> T *floating_ref<T>::operator->() const {

  return reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(Resource_->Target) + TypeOffset_);

}

template <typename T> T &floating_ref<T>::Get() const {

  return *reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(Resource_->Target) + TypeOffset_);

}

template <typename T> template <typename U, OVK_FUNCDEF_REQUIRES(std::is_base_of<T, U>::value)>
  floating_ref<T>::operator floating_ref<U>() const {

  using other_byte_ptr = typename floating_ref<U>::byte_ptr;

  auto Target = reinterpret_cast<T *>(reinterpret_cast<byte_ptr>(Resource_->Target) + TypeOffset_);
  auto OtherTarget = static_cast<U *>(Target);

  floating_ref<U> Other;
  Other.Resource_ = Resource_;
  Other.TypeOffset_ = TypeOffset_ + (reinterpret_cast<other_byte_ptr>(OtherTarget) -
    reinterpret_cast<byte_ptr>(Target));

  return Other;

}

template <typename U, typename T, OVK_FUNCDEF_REQUIRES(std::is_convertible<T *, U *>::value ||
  std::is_base_of<T, U>::value)> floating_ref<U> floating_ref_cast(const floating_ref<T>
  &FloatingRef) {

  return static_cast<floating_ref<U>>(FloatingRef);

}

}
