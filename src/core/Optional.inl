// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {

template <typename T> optional<T>::optional(const value_type &Value):
  Present_(true)
{
  ValueStorage_.Create(Value);
}

template <typename T> optional<T>::optional(value_type &&Value):
  Present_(true)
{
  ValueStorage_.Create(std::move(Value));
}

template <typename T> optional<T>::optional(const optional &Other):
  Present_(Other.Present_)
{
  if (Other.Present_) {
    ValueStorage_.Create(Other.ValueStorage_.Get());
  }
}

template <typename T> optional<T>::optional(optional &&Other) noexcept:
  Present_(Other.Present_)
{
  if (Other.Present_) {
    ValueStorage_.Create(std::move(Other.ValueStorage_.Get()));
    Other.ValueStorage_.Destroy();
  }
  Other.Present_ = false;
}

template <typename T> optional<T>::~optional() noexcept {

  if (Present_) {
    ValueStorage_.Destroy();
  }

}

template <typename T> optional<T> &optional<T>::operator=(const optional &Other) {

  if (&Other == this) { return *this; }

  if (Other.Present_) {
    // Create first in case constructor throws
    value_type TempValue = Other.ValueStorage_.Get();
    if (Present_) {
      ValueStorage_.Destroy();
    }
    ValueStorage_.Create(std::move(TempValue));
  } else {
    if (Present_) {
      ValueStorage_.Destroy();
    }
  }

  Present_ = Other.Present_;

  return *this;

}

template <typename T> optional<T> &optional<T>::operator=(optional &&Other) noexcept {

  if (Present_) {
    ValueStorage_.Destroy();
  }

  if (Other.Present_) {
    ValueStorage_.Create(std::move(Other.ValueStorage_.Get()));
  }

  Present_ = Other.Present_;
  Other.Present_ = false;

  return *this;

}

template <typename T> optional<T>::operator bool() const {

  return Present_;

}

template <typename T> bool optional<T>::Present() const {

  return Present_;

}

template <typename T> optional<T> &optional<T>::Assign(const value_type &Value) {

  // Create first in case constructor throws
  value_type TempValue = Value;

  if (Present_) {
    ValueStorage_.Destroy();
  }

  ValueStorage_.Create(std::move(TempValue));

  Present_ = true;

  return *this;

}

template <typename T> optional<T> &optional<T>::Assign(value_type &&Value) {

  if (Present_) {
    ValueStorage_.Destroy();
  }

  ValueStorage_.Create(std::move(Value));

  Present_ = true;

  return *this;

}

template <typename T> template <typename... Args, OVK_FUNCDEF_REQUIRES(std::is_constructible<T,
  Args &&...>::value && !core::IsCopyOrMoveArgument<T, Args &&...>())> optional<T> &optional<T>::
  Assign(Args &&... Arguments) {

  // Create first in case constructor throws
  value_type TempValue(std::forward<Args>(Arguments)...);

  if (Present_) {
    ValueStorage_.Destroy();
  }

  ValueStorage_.Create(std::move(TempValue));

  Present_ = true;

  return *this;

}

template <typename T> auto optional<T>::operator->() const -> const value_type * {

  return &ValueStorage_.Get();

}

template <typename T> auto optional<T>::operator->() -> value_type * {

  return &ValueStorage_.Get();

}

template <typename T> auto optional<T>::operator*() const -> const value_type & {

  return ValueStorage_.Get();

}

template <typename T> auto optional<T>::operator*() -> value_type & {

  return ValueStorage_.Get();

}

template <typename T> auto optional<T>::Get() const -> const value_type & {

  return ValueStorage_.Get();

}

template <typename T> auto optional<T>::Get() -> value_type & {

  return ValueStorage_.Get();

}

template <typename T> auto optional<T>::Release() -> value_type {

  value_type Value = std::move(ValueStorage_.Get());

  ValueStorage_.Destroy();

  Present_ = false;

  return Value;

}

template <typename T> void optional<T>::Reset() {

  if (Present_) {
    ValueStorage_.Destroy();
  }

  Present_ = false;

}

template <typename T> template <typename... Args> void optional<T>::value_storage::Create(Args &&...
  Arguments) {

  new (&Storage_) value_type(std::forward<Args>(Arguments)...);

}

template <typename T> void optional<T>::value_storage::Destroy() {

  value_type &Value = *reinterpret_cast<value_type *>(&Storage_);
  Value.~value_type();

}

template <typename T> auto optional<T>::value_storage::Get() const -> const value_type & {

  return *reinterpret_cast<const value_type *>(&Storage_);

}

template <typename T> auto optional<T>::value_storage::Get() -> value_type & {

  return *reinterpret_cast<value_type *>(&Storage_);

}

}
