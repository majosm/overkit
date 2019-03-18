// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_GRID_HPP_LOADED
#define OVK_SUPPORT_GRID_HPP_LOADED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/Range.hpp>

#include <mpi.h>

#include <functional>
#include <memory>
#include <string>

namespace support {

class grid {

public:

  class params {
  public:
    params():
      ID_(0),
      Name_("Grid"),
      Cart_(ovk::MakeEmptyCart(2)),
      PeriodicLength_(ovk::MakeUniformTuple<double>(0.)),
      LocalRange_(ovk::MakeEmptyRange(2)),
      ExtendAmount_(0),
      Layout_(ovk::array_layout::GRID)
    {}
    params &SetID(int ID) { ID_ = ID; return *this; }
    params &SetName(std::string Name) { Name_ = std::move(Name); return *this; }
    params &SetCart(const ovk::cart &Cart) { Cart_ = Cart; return *this; }
    params &SetComm(ovk::core::comm Comm) { Comm_ = std::move(Comm); return *this; }
    params &SetPeriodicLength(const ovk::tuple<double> &PeriodicLength) {
      PeriodicLength_ = PeriodicLength;
      return *this;
    }
    params &SetLocalRange(const ovk::range &LocalRange) { LocalRange_ = LocalRange; return *this; }
    params &SetExtendAmount(int ExtendAmount) { ExtendAmount_ = ExtendAmount; return *this; }
    params &SetLayout(ovk::array_layout Layout) { Layout_ = Layout; return *this; }
    template <typename F> params &AddManipulator(F Manipulator) {
      Manipulators_.Append(std::move(Manipulator));
      return *this;
    }
  private:
    int ID_;
    std::string Name_;
    ovk::cart Cart_;
    ovk::core::comm Comm_;
    ovk::tuple<double> PeriodicLength_;
    ovk::range LocalRange_;
    int ExtendAmount_;
    ovk::array_layout Layout_;
    ovk::array<std::function<void(grid &)>> Manipulators_;
    friend class grid;
  };

  explicit grid(params Params);

  int ID() const { return ID_; }

  const std::string &Name() const { return Name_; }

  const ovk::cart &Cart() const { return Cart_; }

  ovk::core::comm_view Comm() const { return Comm_; }

  const ovk::tuple<double> &PeriodicLength() const { return PeriodicLength_; }
  double PeriodicLength(int iDim) const { return PeriodicLength_[iDim]; }

  const ovk::range &LocalRange() const { return LocalRange_; }
  const ovk::range &ExtendedRange() const { return ExtendedRange_; }

  ovk::array_layout Layout() const { return Layout_; }

  template <ovk::array_layout Layout=ovk::array_layout::GRID> ovk::array_view<const ovk::array<
    double,OVK_MAX_DIMS,Layout>> XYZ() const {
    OVK_DEBUG_ASSERT(Layout == Layout_, "Wrong grid layout.");
    return static_cast<const fields<Layout> *>(Fields_.get())->XYZ;
  }
  template <ovk::array_layout Layout=ovk::array_layout::GRID> ovk::array_view<ovk::array<double,
    OVK_MAX_DIMS,Layout>> XYZ() {
    OVK_DEBUG_ASSERT(Layout == Layout_, "Wrong grid layout.");
    return static_cast<fields<Layout> *>(Fields_.get())->XYZ;
  }

  template <ovk::array_layout Layout=ovk::array_layout::GRID> ovk::array_view<const int,
    OVK_MAX_DIMS,Layout> IBlank() const {
    OVK_DEBUG_ASSERT(Layout == Layout_, "Wrong grid layout.");
    return static_cast<const fields<Layout> *>(Fields_.get())->IBlank;
  }
  template <ovk::array_layout Layout=ovk::array_layout::GRID> ovk::array_view<int,OVK_MAX_DIMS,
    Layout> IBlank() {
    OVK_DEBUG_ASSERT(Layout == Layout_, "Wrong grid layout.");
    return static_cast<fields<Layout> *>(Fields_.get())->IBlank;
  }

private:

  struct fields_base {
    virtual ~fields_base() {}
  };

  template <ovk::array_layout Layout> struct fields : fields_base {
    ovk::array<ovk::array<double,OVK_MAX_DIMS,Layout>> XYZ;
    ovk::array<int,OVK_MAX_DIMS,Layout> IBlank;
  };

  int ID_;
  std::string Name_;
  ovk::cart Cart_;
  ovk::core::comm Comm_;
  ovk::tuple<double> PeriodicLength_;
  ovk::range LocalRange_;
  ovk::range ExtendedRange_;
  ovk::array_layout Layout_;
  std::unique_ptr<fields_base> Fields_;

};

}

#endif
