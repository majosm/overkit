// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_POINTER_ITERATOR_HPP_LOADED
#define OVK_CORE_POINTER_ITERATOR_HPP_LOADED

#include <ovk/core/IteratorTraits.hpp>
#include <ovk/core/Requires.hpp>

#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utility>

namespace ovk {
namespace core {

// Iterator that wraps a pointer (sometimes can't use a pointer directly due to, e.g., function
// overload ambiguity issues)

template <typename ContainerType, typename PointerType> class pointer_iterator {

public:

  using difference_type = std::ptrdiff_t;
  using value_type = iterator_value_type<PointerType>;
  using pointer = PointerType;
  using reference = iterator_reference_type<PointerType>;
  using iterator_category = std::random_access_iterator_tag;

  constexpr pointer_iterator():
    Pointer_(pointer())
  {}

  explicit constexpr pointer_iterator(const pointer &Pointer):
    Pointer_(Pointer)
  {}

  template <typename OtherPointerType, OVK_FUNCTION_REQUIRES(!std::is_same<OtherPointerType,
    PointerType>::value && std::is_convertible<OtherPointerType, PointerType>::value)> constexpr
    pointer_iterator(const pointer_iterator<ContainerType, OtherPointerType> &Other):
    Pointer_(Other.Pointer())
  {}

  constexpr reference operator*() const { return *Pointer_; }

  constexpr pointer operator->() const { return Pointer_; }

  pointer_iterator &operator++() { ++Pointer_; return *this; }
  pointer_iterator operator++(int) { return pointer_iterator(Pointer_++); }

  pointer_iterator &operator--() { --Pointer_; return *this; }
  pointer_iterator operator--(int) { return pointer_iterator(Pointer_--); }

  constexpr reference operator[](difference_type i) const { return Pointer_[i]; }

  pointer_iterator &operator+=(difference_type i) { Pointer_ += i; return *this; }
  pointer_iterator &operator-=(difference_type i) { Pointer_ -= i; return *this; }

  constexpr pointer Pointer() const { return Pointer_; }

private:

  pointer Pointer_;

};

template <typename ContainerType, typename PointerType> constexpr pointer_iterator<ContainerType,
  PointerType> operator+(const pointer_iterator<ContainerType, PointerType> &Iter, std::ptrdiff_t i)
  {
  return pointer_iterator<ContainerType, PointerType>(Iter.Pointer() + i);
}

template <typename ContainerType, typename PointerType> constexpr pointer_iterator<ContainerType,
  PointerType> operator+(std::ptrdiff_t i, const pointer_iterator<ContainerType, PointerType> &Iter)
  {
  return pointer_iterator<ContainerType, PointerType>(Iter.Pointer() + i);
}

template <typename ContainerType, typename PointerType> constexpr pointer_iterator<ContainerType,
  PointerType> operator-(const pointer_iterator<ContainerType, PointerType> &Iter, std::ptrdiff_t i)
  {
  return pointer_iterator<ContainerType, PointerType>(Iter.Pointer() - i);
}

template <typename ContainerType, typename PointerType1, typename PointerType2,
  OVK_FUNCTION_REQUIRES(std::is_convertible<PointerType1, PointerType2>::value ||
  std::is_convertible<PointerType2, PointerType1>::value)> constexpr std::ptrdiff_t operator-(const
  pointer_iterator<ContainerType, PointerType1> &Left, const pointer_iterator<ContainerType,
  PointerType2> &Right) {
  return Left.Pointer() - Right.Pointer();
}

template <typename ContainerType, typename PointerType1, typename PointerType2,
  OVK_FUNCTION_REQUIRES(std::is_convertible<PointerType1, PointerType2>::value ||
  std::is_convertible<PointerType2, PointerType1>::value)> constexpr bool operator==(const
  pointer_iterator<ContainerType, PointerType1> &Left, const pointer_iterator<ContainerType,
  PointerType2> &Right) {
  return Left.Pointer() == Right.Pointer();
}

template <typename ContainerType, typename PointerType1, typename PointerType2,
  OVK_FUNCTION_REQUIRES(std::is_convertible<PointerType1, PointerType2>::value ||
  std::is_convertible<PointerType2, PointerType1>::value)> constexpr bool operator!=(const
  pointer_iterator<ContainerType, PointerType1> &Left, const pointer_iterator<ContainerType,
  PointerType2> &Right) {
  return Left.Pointer() != Right.Pointer();
}

template <typename ContainerType, typename PointerType1, typename PointerType2,
  OVK_FUNCTION_REQUIRES(std::is_convertible<PointerType1, PointerType2>::value ||
  std::is_convertible<PointerType2, PointerType1>::value)> constexpr bool operator<(const
  pointer_iterator<ContainerType, PointerType1> &Left, const pointer_iterator<ContainerType,
  PointerType2> &Right) {
  return Left.Pointer() < Right.Pointer();
}

template <typename ContainerType, typename PointerType1, typename PointerType2,
  OVK_FUNCTION_REQUIRES(std::is_convertible<PointerType1, PointerType2>::value ||
  std::is_convertible<PointerType2, PointerType1>::value)> constexpr bool operator>(const
  pointer_iterator<ContainerType, PointerType1> &Left, const pointer_iterator<ContainerType,
  PointerType2> &Right) {
  return Left.Pointer() > Right.Pointer();
}

template <typename ContainerType, typename PointerType1, typename PointerType2,
  OVK_FUNCTION_REQUIRES(std::is_convertible<PointerType1, PointerType2>::value ||
  std::is_convertible<PointerType2, PointerType1>::value)> constexpr bool operator<=(const
  pointer_iterator<ContainerType, PointerType1> &Left, const pointer_iterator<ContainerType,
  PointerType2> &Right) {
  return Left.Pointer() <= Right.Pointer();
}

template <typename ContainerType, typename PointerType1, typename PointerType2,
  OVK_FUNCTION_REQUIRES(std::is_convertible<PointerType1, PointerType2>::value ||
  std::is_convertible<PointerType2, PointerType1>::value)> constexpr bool operator>=(const
  pointer_iterator<ContainerType, PointerType1> &Left, const pointer_iterator<ContainerType,
  PointerType2> &Right) {
  return Left.Pointer() >= Right.Pointer();
}

}}

#endif
