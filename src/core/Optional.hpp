// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_OPTIONAL_HPP_LOADED
#define OVK_CORE_OPTIONAL_HPP_LOADED

#include <ovk/core/Global.hpp>
#include <ovk/core/Requires.hpp>
#include <ovk/core/TypeTraits.hpp>

#include <type_traits>
#include <utility>

namespace ovk {

template <typename T> class optional {

public:

  static_assert(std::is_same<core::remove_cvref<T>, T>::value, "Type parameter in optional cannot "
    "be cv-qualified or a reference type.");

  using value_type = T;

  constexpr optional() = default;

  optional(const value_type &Value);
  optional(value_type &&Value);
  // No forwarding constructor (would potentially cause ambiguities with default constructor)

  optional(const optional &Other);
  optional(optional &&Other) noexcept;

  ~optional() noexcept;

  optional &operator=(const optional &Other);
  optional &operator=(optional &&Other) noexcept;

  explicit operator bool() const;
  bool Present() const;

  optional &Assign(const value_type &Value);
  optional &Assign(value_type &&Value);
  template <typename... Args, OVK_FUNCDECL_REQUIRES(std::is_constructible<value_type, Args &&...
    >::value && !core::IsCopyOrMoveArgument<value_type, Args &&...>())> optional &Assign(Args &&...
    Arguments);

  const value_type *operator->() const;
  value_type *operator->();

  const value_type &operator*() const;
  value_type &operator*();

  const value_type &Get() const;
  value_type &Get();

  value_type Release();

  void Clear();

private:

  class value_storage {
  public:
    template <typename... Args> void Create(Args &&... Arguments);
    void Destroy();
    const value_type &Get() const;
    value_type &Get();
  private:
    typename std::aligned_storage<sizeof(value_type), alignof(value_type)>::type Storage_;
  };

  value_storage ValueStorage_;
  bool Present_ = false;

  friend class core::test_helper<optional>;

};

}

#include <ovk/core/Optional.inl>

#endif
