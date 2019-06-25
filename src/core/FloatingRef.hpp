// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
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

template <typename T> class floating_ref;

class floating_ref_generator {

public:

  floating_ref_generator();

  floating_ref_generator(const floating_ref_generator &Other);
  floating_ref_generator(floating_ref_generator &&Other) noexcept;

  floating_ref_generator &operator=(const floating_ref_generator &Other);
  floating_ref_generator &operator=(floating_ref_generator &&Other) noexcept;

  template <typename T> floating_ref<T> Generate(T &Target) const;

private:

  std::unique_ptr<void *> ReferenceLoc_;

};

template <typename T> class floating_ref {

public:

  floating_ref() = default;

  template <typename U, OVK_FUNCDECL_REQUIRES(!std::is_same<U, T>::value && std::is_convertible<
    U *, T *>::value)> floating_ref(const floating_ref<U> &Other);

  floating_ref(const floating_ref &Other) = default;
  floating_ref(floating_ref &&Other) noexcept;

  floating_ref &operator=(const floating_ref &Other) = default;
  floating_ref &operator=(floating_ref &&Other) noexcept;

  explicit operator bool() const { return ReferenceLoc_ != nullptr; }

  T &operator*() const;

  T *operator->() const;

  T &Get() const;

  template <typename U, OVK_FUNCDECL_REQUIRES(!std::is_same<U, T>::value && !std::is_convertible<
    T *, U *>::value && std::is_base_of<T, U>::value)> explicit operator floating_ref<U>() const;

  void Reset() { ReferenceLoc_ = nullptr; }

private:

  using byte_ptr = core::mimic_cvref<T, unsigned char> *;

  void * const *ReferenceLoc_ = nullptr;
  std::ptrdiff_t Offset_ = 0;

  floating_ref(void * const *ReferenceLoc, T &Target);

  friend class floating_ref_generator;
  template <typename U> friend class floating_ref;

};

template <typename U, typename T, OVK_FUNCDECL_REQUIRES(std::is_convertible<T *, U *>::value ||
  std::is_base_of<T, U>::value)> floating_ref<U> floating_ref_cast(const floating_ref<T>
  &FloatingRef);

}

#include <ovk/core/FloatingRef.inl>

#endif
