// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_FLOATING_REF_HPP_INCLUDED
#define OVK_CORE_FLOATING_REF_HPP_INCLUDED

#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <cstdint>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace ovk {

namespace floating_ref_internal {

struct resource {
  void *Target;
  explicit resource(void *Target_):
    Target(Target_)
  {}
};

}

template <typename T> class floating_ref;

template <typename T> class floating_ref_generator {

public:

  floating_ref_generator(T &Target);

  floating_ref_generator(const floating_ref_generator &Other);
  floating_ref_generator(floating_ref_generator &&Other) noexcept;

  floating_ref_generator &operator=(const floating_ref_generator &Other);
  floating_ref_generator &operator=(floating_ref_generator &&Other) noexcept;

  T &Target() const { return *static_cast<T *>(Resource_->Target); }

  template <typename U=T, OVK_FUNCDECL_REQUIRES(std::is_convertible<T *, U *>::value ||
    std::is_base_of<T, U>::value)> floating_ref<U> Generate() const;

private:

  using resource = floating_ref_internal::resource;

  using byte_ptr = core::mimic_cvref<T, unsigned char> *;

  std::unique_ptr<resource> Resource_;

};

template <typename T> class floating_ref {

public:

  floating_ref() = default;

  template <typename U, OVK_FUNCDECL_REQUIRES(std::is_convertible<U *, T *>::value)>
    floating_ref(const floating_ref<U> &Other);

  floating_ref(const floating_ref &Other) = default;
  floating_ref(floating_ref &&Other) noexcept;

  floating_ref &operator=(const floating_ref &Other) = default;
  floating_ref &operator=(floating_ref &&Other) noexcept;

  explicit operator bool() const { return Resource_ != nullptr; }

  T &operator*() const;

  T *operator->() const;

  T &Get() const;

  void Reset() { Resource_ = nullptr; }

  template <typename U, OVK_FUNCDECL_REQUIRES(std::is_base_of<T, U>::value)> explicit operator
    floating_ref<U>() const;

private:

  using resource = floating_ref_internal::resource;

  using byte_ptr = core::mimic_cvref<T, unsigned char> *;

  resource *Resource_ = nullptr;
  std::ptrdiff_t TypeOffset_ = 0;

  explicit floating_ref(resource &Resource):
    Resource_(&Resource)
  {}

  template <typename U> friend class floating_ref;
  template <typename U> friend class floating_ref_generator;

};

template <typename U, typename T, OVK_FUNCDECL_REQUIRES(std::is_convertible<T *, U *>::value ||
  std::is_base_of<T, U>::value)> floating_ref<U> floating_ref_cast(const floating_ref<T>
  &FloatingRef);

}

#include <ovk/core/FloatingRef.inl>

#endif
