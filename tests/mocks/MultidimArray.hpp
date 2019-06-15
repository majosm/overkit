// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_TESTS_MOCKS_MULTIDIM_ARRAY_HPP_LOADED
#define OVK_TESTS_MOCKS_MULTIDIM_ARRAY_HPP_LOADED

#include <ovk/core/ArrayTraits.hpp>
#include <ovk/core/Interval.hpp>

#include <initializer_list>
#include <vector>

namespace tests {

// Minimal multidimensional array type
namespace multidim_array_internal {
template <typename T> class base {
public:
  base(const ovk::interval<long long,3> &Extents, T Value=T()):
    Values_(Extents.Count(), Value),
    Extents_(Extents)
  {}
  base(const ovk::interval<long long,3> &Extents, std::initializer_list<T> ValuesList):
    Values_(ValuesList),
    Extents_(Extents)
  {}
  long long Begin(int iDim) const { return Extents_.Begin(iDim); }
  long long End(int iDim) const { return Extents_.End(iDim); }
  const T &operator[](long long l) const { return Values_[l]; }
  T &operator[](long long l) { return Values_[l]; }
protected:
  std::vector<T> Values_;
  ovk::interval<long long,3> Extents_;
};
}

// Implementations are the same, just treated differently by array_traits
template <typename T> class multidim_array_row : public multidim_array_internal::base<T> {
protected:
  using parent_type = multidim_array_internal::base<T>;
  using parent_type::Values_;
  using parent_type::Extents_;
public:
  using parent_type::parent_type;
  const T &operator()(int i, int j, int k) const {
    return Values_[Extents_.Size(0)*(Extents_.Size(1)*(i-Extents_.Begin(0)) + (j-Extents_.Begin(1)))
      + (k-Extents_.Begin(2))];
  }
  T &operator()(int i, int j, int k) {
    return Values_[Extents_.Size(0)*(Extents_.Size(1)*(i-Extents_.Begin(0)) + (j-Extents_.Begin(1)))
      + (k-Extents_.Begin(2))];
  }
};
template <typename T> class multidim_array_col : public multidim_array_internal::base<T> {
protected:
  using parent_type = multidim_array_internal::base<T>;
  using parent_type::Values_;
  using parent_type::Extents_;
public:
  using parent_type::parent_type;
  const T &operator()(int i, int j, int k) const {
    return Values_[(i-Extents_.Begin(0)) + Extents_.Size(0)*((j-Extents_.Begin(1)) +
      Extents_.Size(1)*(k-Extents_.Begin(2)))];
  }
  T &operator()(int i, int j, int k) {
    return Values_[(i-Extents_.Begin(0)) + Extents_.Size(0)*((j-Extents_.Begin(1)) +
      Extents_.Size(1)*(k-Extents_.Begin(2)))];
  }
};

}

namespace ovk {

template <typename T> struct array_traits<tests::multidim_array_row<T>> {
private:
  using array_type = tests::multidim_array_row<T>;
public:
  using value_type = T;
  static constexpr int Rank = 3;
  static constexpr array_layout Layout = array_layout::ROW_MAJOR;
  template <int iDim> static long long ExtentBegin(const array_type &Array) {
    return Array.Begin(iDim);
  }
  template <int iDim> static long long ExtentEnd(const array_type &Array) {
    return Array.End(iDim);
  }
  static const value_type *Data(const array_type &Array) { return &Array[0]; }
  static value_type *Data(array_type &Array) { return &Array[0]; }
};

template <typename T> struct array_traits<tests::multidim_array_col<T>> {
private:
  using array_type = tests::multidim_array_col<T>;
public:
  using value_type = T;
  static constexpr int Rank = 3;
  static constexpr array_layout Layout = array_layout::COLUMN_MAJOR;
  template <int iDim> static long long ExtentBegin(const array_type &Array) {
    return Array.Begin(iDim);
  }
  template <int iDim> static long long ExtentEnd(const array_type &Array) {
    return Array.End(iDim);
  }
  static const value_type *Data(const array_type &Array) { return &Array[0]; }
  static value_type *Data(array_type &Array) { return &Array[0]; }
};

}

#endif
