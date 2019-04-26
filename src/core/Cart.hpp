// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CART_HPP_INCLUDED
#define OVK_CORE_CART_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

namespace ovk {

class cart {

public:

  // Remove this eventually?
  cart() = default;

  cart(int NumDims);
  cart(int NumDims, const range &Range, const tuple<bool> &Periodic, periodic_storage
    PeriodicStorage);

  int Dimension() const { return NumDims_; }

  const range &Range() const { return Range_; }
  range &Range() { return Range_; }

  const tuple<bool> &Periodic() const { return Periodic_; }
  tuple<bool> &Periodic() { return Periodic_; }
  const bool &Periodic(int iDim) const { return Periodic_[iDim]; }
  bool &Periodic(int iDim) { return Periodic_[iDim]; }

  const periodic_storage &PeriodicStorage() const { return PeriodicStorage_; }
  periodic_storage &PeriodicStorage() { return PeriodicStorage_; }

  tuple<int> GetPeriod(const tuple<int> &Point) const;
  tuple<int> PeriodicAdjust(const tuple<int> &Point) const;

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
inline cart CartPointToCell(const cart &PointCart);

}

#include <ovk/core/Cart.inl>

#endif
