// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CART_HPP_INCLUDED
#define OVK_CORE_CART_HPP_INCLUDED

#include <ovk/core/Cart.h>
#include <ovk/core/DataType.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Optional.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

namespace ovk {

enum class periodic_storage : typename std::underlying_type<ovk_periodic_storage>::type {
  UNIQUE = OVK_PERIODIC_STORAGE_UNIQUE,
  DUPLICATED = OVK_PERIODIC_STORAGE_DUPLICATED
};

inline bool ValidPeriodicStorage(periodic_storage PeriodicStorage) {
  return ovkValidPeriodicStorage(ovk_periodic_storage(PeriodicStorage));
}

namespace core {
template <> struct data_type_traits<periodic_storage> : data_type_traits<typename
  std::underlying_type<periodic_storage>::type> {};
}

class cart {

public:

  cart(int NumDims);
  cart(int NumDims, const range &Range, const tuple<bool> &Periodic, periodic_storage
    PeriodicStorage);

  int Dimension() const { return NumDims_; }

  const range &Range() const { return Range_; }
  range &Range() { return Range_; }

  const tuple<bool> &Periodic() const { return Periodic_; }
  tuple<bool> &Periodic() { return Periodic_; }
  const bool &Periodic(int iDim) const { return Periodic_(iDim); }
  bool &Periodic(int iDim) { return Periodic_(iDim); }

  const periodic_storage &PeriodicStorage() const { return PeriodicStorage_; }
  periodic_storage &PeriodicStorage() { return PeriodicStorage_; }

  tuple<int> GetPeriodSize() const;
  tuple<int> GetPeriod(const tuple<int> &Point) const;
  tuple<int> PeriodicAdjust(const tuple<int> &Point) const;

  optional<tuple<int>> MapToRange(const range &Range, const tuple<int> &Tuple) const;
  optional<range> MapToRange(const range &Range, const range &OtherRange) const;

private:

  int NumDims_;
  range Range_;
  tuple<bool> Periodic_;
  periodic_storage PeriodicStorage_;

  friend bool operator==(const cart &Left, const cart &Right);

};

inline bool operator==(const cart &Left, const cart &Right);
inline bool operator!=(const cart &Left, const cart &Right);

inline cart MakeEmptyCart(int NumDims);

}

#include <ovk/core/Cart.inl>

#endif
