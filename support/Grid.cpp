// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "support/Grid.hpp"

#include <ovk/core/Array.hpp>
#include <ovk/core/Cart.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Partition.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

#include <memory>
#include <utility>

namespace support {

grid::grid(params Params):
  ID_(Params.ID_),
  Name_(std::move(Params.Name_)),
  Cart_(Params.Cart_),
  Comm_(std::move(Params.Comm_)),
  PeriodicLength_(Params.PeriodicLength_),
  LocalRange_(Params.LocalRange_),
  ExtendedRange_(ovk::core::ExtendLocalRange(Cart_, LocalRange_, Params.ExtendAmount_)),
  Layout_(Params.Layout_)
{

  switch (Layout_) {
  case ovk::array_layout::ROW_MAJOR:
    {
      constexpr const ovk::array_layout RowMajor = ovk::array_layout::ROW_MAJOR;
      Fields_.reset(new fields<RowMajor>());
      auto &Fields = static_cast<fields<RowMajor> &>(*Fields_);
      Fields.XYZ.Resize({3});
      for (int iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
        Fields.XYZ(iDim).Resize(ExtendedRange_);
        for (int k = ExtendedRange_.Begin(2); k < ExtendedRange_.End(2); ++k) {
          for (int j = ExtendedRange_.Begin(1); j < ExtendedRange_.End(1); ++j) {
            for (int i = ExtendedRange_.Begin(0); i < ExtendedRange_.End(0); ++i) {
              ovk::tuple<int> Point = {i,j,k};
              Fields.XYZ(iDim)(Point) = double(Point[iDim]-1);
            }
          }
        }
      }
      Fields.IBlank.Resize(ExtendedRange_, 1);
    }
    break;
  case ovk::array_layout::COLUMN_MAJOR:
    {
      constexpr const ovk::array_layout ColumnMajor = ovk::array_layout::COLUMN_MAJOR;
      Fields_.reset(new fields<ColumnMajor>());
      auto &Fields = static_cast<fields<ColumnMajor> &>(*Fields_);
      Fields.XYZ.Resize({3});
      for (int iDim = 0; iDim < OVK_MAX_DIMS; ++iDim) {
        Fields.XYZ(iDim).Resize(ExtendedRange_);
        for (int k = ExtendedRange_.Begin(2); k < ExtendedRange_.End(2); ++k) {
          for (int j = ExtendedRange_.Begin(1); j < ExtendedRange_.End(1); ++j) {
            for (int i = ExtendedRange_.Begin(0); i < ExtendedRange_.End(0); ++i) {
              ovk::tuple<int> Point = {i,j,k};
              Fields.XYZ(iDim)(Point) = double(Point[iDim]-1);
            }
          }
        }
      }
      Fields.IBlank.Resize(ExtendedRange_, 1);
    }
    break;
  }

  for (auto &Manipulator : Params.Manipulators_) {
    std::move(Manipulator)(*this);
  }

}

}
