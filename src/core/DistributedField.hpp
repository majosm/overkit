// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_DISTRIBUTED_FIELD_HPP_INCLUDED
#define OVK_CORE_DISTRIBUTED_FIELD_HPP_INCLUDED

#include <ovk/core/Cart.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Partition.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <mpi.h>

#include <memory>
#include <utility>

namespace ovk {

template <typename T> class distributed_field {

public:

  using value_type = T;
  using index_type = long long;
  using tuple_element_type = int;
  using tuple_type = tuple<int>;
  using interval_type = range;
  using indexer_type = range_indexer_c<long long>;
  using iterator = core::pointer_iterator<distributed_field, value_type *>;
  using const_iterator = core::pointer_iterator<distributed_field, const value_type *>;

  distributed_field() = default;

  explicit distributed_field(std::shared_ptr<const partition> Partition):
    Partition_(std::move(Partition)),
    Values_(Partition_->ExtendedRange())
  {}

  distributed_field(std::shared_ptr<const partition> Partition, const value_type &Value):
    Partition_(std::move(Partition)),
    Values_(Partition_->ExtendedRange(), Value)
  {}

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsInputIterator<IterType>())>
    distributed_field(std::shared_ptr<const partition> Partition, IterType First):
    Partition_(std::move(Partition)),
    Values_(Partition_->ExtendedRange(), First)
  {}

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U>::
    type, value_type>::value)> distributed_field(std::shared_ptr<const partition> Partition,
      const field_view<U> &View):
    Partition_(std::move(Partition)),
    Values_(Partition_->ExtendedRange(), View)
  {}

  template <typename FieldRefType, OVK_FUNCTION_REQUIRES(core::IsField<core::remove_cvref<
    FieldRefType>>() && !std::is_same<core::remove_cvref<FieldRefType>, distributed_field>::value &&
    !core::IsIterator<typename std::decay<FieldRefType>::type>() && std::is_convertible<
    core::array_access_type<FieldRefType &&>, value_type>::value)> distributed_field(
    std::shared_ptr<const partition> Partition, FieldRefType &&Field):
    Partition_(std::move(Partition)),
    Values_(Partition_->ExtendedRange(), std::forward<FieldRefType>(Field))
  {}

  distributed_field &Assign(std::shared_ptr<const partition> Partition) {
    Partition_ = std::move(Partition);
    Values_.Clear();
    Values_.Resize(Partition_->ExtendedRange());
    return *this;
  }

  distributed_field &Assign(std::shared_ptr<const partition> Partition, const value_type &Value) {
    Partition_ = std::move(Partition);
    Values_.Assign(Partition_->ExtendedRange(), Value);
    return *this;
  }

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsInputIterator<IterType>())>
    distributed_field &Assign(std::shared_ptr<const partition> Partition, IterType First) {
    Partition_ = std::move(Partition);
    Values_.Assign(Partition_->ExtendedRange(), First);
    return *this;
  }

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U>::
    type, value_type>::value)> distributed_field &Assign(std::shared_ptr<const partition> Partition,
    const field_view<U> &View) {
    Partition_ = std::move(Partition);
    Values_.Assign(Partition_->ExtendedRange(), View);
    return *this;
  }

  template <typename FieldRefType, OVK_FUNCTION_REQUIRES(core::IsField<core::remove_cvref<
    FieldRefType>>() && !std::is_same<core::remove_cvref<FieldRefType>, distributed_field>::value &&
    !core::IsIterator<typename std::decay<FieldRefType>::type>() && std::is_convertible<
    core::array_access_type<FieldRefType &&>, value_type>::value)> distributed_field &Assign(
    std::shared_ptr<const partition> Partition, FieldRefType &&Field) {
    Partition_ = std::move(Partition);
    Values_.Assign(Partition_->ExtendedRange(), std::forward<FieldRefType>(Field));
    return *this;
  }

  distributed_field &Assign(const distributed_field &Other) {
    Partition_ = Other.Partition_;
    Values_.Assign(Other.Values_);
    return *this;
  }

  distributed_field &Assign(distributed_field &&Other) noexcept {
    Partition_ = std::move(Other.Partition_);
    Values_.Assign(std::move(Other.Values_));
    return *this;
  }

  request Exchange() { return Partition_->Exchange(Values_); }

  distributed_field &Fill(const value_type &Value) {
    Values_.Fill(Value);
    return *this;
  }

  distributed_field &Fill(const range &Range, const value_type &Value) {
    const cart &Cart = Partition_->Cart();
    const range &LocalRange = Partition_->LocalRange();
    if (Cart.Range().Includes(Range)) {
      range IntersectRange = IntersectRanges(LocalRange, Range);
      Values_.Fill(IntersectRange, Value);
    } else {
      for (int k = Range.Begin(2); k < Range.End(2); ++k) {
        for (int j = Range.Begin(1); j < Range.End(1); ++j) {
          for (int i = Range.Begin(0); i < Range.End(0); ++i) {
            tuple<int> Point = {i,j,k};
            if (LocalRange.Contains(Point)) {
              Values_(Point) = Value;
            } else {
              Point = Cart.PeriodicAdjust(Point);
              if (LocalRange.Contains(Point)) {
                Values_(Point) = Value;
              }
            }
          }
        }
      }
    }
    Exchange();
    return *this;
  }

  template <typename IterType, OVK_FUNCTION_REQUIRES(core::IsInputIterator<IterType>() &&
    std::is_convertible<core::iterator_reference_type<IterType>, value_type>::value)>
    distributed_field &Fill(IterType First) {
    Values_.Fill(First);
    return *this;
  }

  template <typename U, OVK_FUNCTION_REQUIRES(std::is_convertible<typename std::remove_const<U>::
    type, value_type>::value)> distributed_field &Fill(const field_view<U> &View) {
    Values_.Fill(View);
    return *this;
  }

  template <typename FieldRefType, OVK_FUNCTION_REQUIRES(core::IsField<core::remove_cvref<
    FieldRefType>>() && !core::IsIterator<typename std::decay<FieldRefType>::type>() &&
    std::is_convertible<core::array_access_type<FieldRefType &&>, value_type>::value)>
    distributed_field &Fill(FieldRefType &&Field) {
    Values_.Fill(std::forward<FieldRefType>(Field));
    return *this;
  }

  explicit operator bool() const { return static_cast<bool>(Partition_); }

  const partition &Partition() const { return *Partition_; }
  const std::shared_ptr<const partition> &SharedPartition() const { return Partition_; }

  const cart &Cart() const { return Partition_->Cart(); }

  comm_view Comm() const { return Partition_->Comm(); }

  const range &GlobalRange() const { return Partition_->GlobalRange(); }
  const range &LocalRange() const { return Partition_->LocalRange(); }
  const range &ExtendedRange() const { return Partition_->ExtendedRange(); }

  const range &Extents() const { return Partition_->ExtendedRange(); }

  tuple<int> Size() const { return Values_.Size(); }
  int Size(int iDim) const { return Values_.Size(iDim); }

  long long Count() const { return Values_.Count(); }

  // Can't do this yet (field currently has different indexer type; may be able to fix this)
//   const indexer_type &Indexer() const { return Values_.Indexer(); }

  const field<value_type> &Values() const { return Values_; }
  field<value_type> &Values() { return Values_; }

  const value_type &operator()(const tuple<int> &Tuple) const { return Values_(Tuple); }
  value_type &operator()(const tuple<int> &Tuple) { return Values_(Tuple); }
  const value_type &operator()(int i, int j, int k) const { return Values_(i,j,k); }
  value_type &operator()(int i, int j, int k) { return Values_(i,j,k); }

  const value_type &operator[](long long iValue) const { return Values_[iValue]; }
  value_type &operator[](long long iValue) { return Values_[iValue]; }

  const value_type *Data() const { return Values_.Data(); }
  value_type *Data() { return Values_.Data(); }

  const value_type *Data(const tuple<int> &Tuple) const { return Values_.Data(Tuple); }
  value_type *Data(const tuple<int> &Tuple) { return Values_.Data(Tuple); }
  const value_type *Data(int i, int j, int k) const { return Values_.Data(i,j,k); }
  value_type *Data(int i, int j, int k) { return Values_.Data(i,j,k); }

  const_iterator Begin() const { return const_iterator(Values_.Data()); }
  iterator Begin() { return iterator(Values_.Data()); }
  const_iterator End() const { return const_iterator(Values_.Data()+Values_.Count()); }
  iterator End() { return iterator(Values_.Data()+Values_.Count()); }

  // Google Test doesn't use free begin/end functions and instead expects container to have
  // lowercase begin/end methods
  const_iterator begin() const { return Begin(); }
  iterator begin() { return Begin(); }
  const_iterator end() const { return End(); }
  iterator end() { return End(); }

private:

  std::shared_ptr<const partition> Partition_;
  field<value_type> Values_;

  friend class core::test_helper<distributed_field>;

};

template <typename T> typename distributed_field<T>::iterator begin(distributed_field<T> &Field) {
  return Field.Begin();
}

template <typename T> typename distributed_field<T>::const_iterator begin(const distributed_field<T>
  &Field) {
  return Field.Begin();
}

template <typename T> typename distributed_field<T>::iterator end(distributed_field<T> &Field) {
  return Field.End();
}

template <typename T> typename distributed_field<T>::const_iterator end(const distributed_field<T>
  &Field) {
  return Field.End();
}

template <typename T> struct array_traits<distributed_field<T>> {
  using value_type = T;
  static constexpr int Rank = MAX_DIMS;
  static constexpr array_layout Layout = array_layout::COLUMN_MAJOR;
  template <int iDim> static long long ExtentBegin(const distributed_field<T> &Field) {
    return Field.Extents().Begin(iDim);
  }
  template <int iDim> static long long ExtentEnd(const distributed_field<T> &Field) {
    return Field.Extents().End(iDim);
  }
  static const T *Data(const distributed_field<T> &Field) { return Field.Data(); }
  static T *Data(distributed_field<T> &Field) { return Field.Data(); }
};

}

#endif
