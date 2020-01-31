// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GEOMETRY_MANIPULATOR_HPP_INCLUDED
#define OVK_CORE_GEOMETRY_MANIPULATOR_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/GeometryOps.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/Tuple.hpp>

namespace ovk {
namespace core {

namespace geometry_manipulator_internal {

template <geometry_type Type_, int NumDims> struct geometry_manipulator_for_type_and_dim {
  constexpr geometry_type Type() const { return Type_; }
  bool OverlapsCell(const array_view<const field_view<const double>> &Coords, double Tolerance,
    const tuple<int> &Cell, const tuple<double> &PointCoords) const {
    return core::OverlapsCell<Type_, NumDims>(Coords, Tolerance, Cell, PointCoords);
  }
  optional<tuple<double>> CoordsInCell(const array_view<const field_view<const double>> &Coords,
    const tuple<int> &Cell, const tuple<double> &PointCoords) const {
    return core::CoordsInCell<Type_, NumDims>(Coords, Cell, PointCoords);
  }
  double CellVolume(const array_view<const field_view<const double>> &Coords, const tuple<int>
    &Cell) const {
    return core::CellVolume<Type_, NumDims>(Coords, Cell);
  }
  box CellBounds(const array_view<const field_view<const double>> &Coords, const tuple<int> &Cell)
    const {
    return core::CellBounds<Type_, NumDims>(Coords, Cell);
  }
};

}

class geometry_manipulator {

public:

  template <geometry_type Type, int NumDims> using manipulator_type =
    geometry_manipulator_internal::geometry_manipulator_for_type_and_dim<Type, NumDims>;

  geometry_manipulator(geometry_type Type, int NumDims):
    Type_(Type),
    NumDims_(NumDims)
  {}

  template <typename F, typename... Args> auto Apply(F &&Func, Args &&... Arguments) const ->
    decltype(std::forward<F>(Func)(manipulator_type<geometry_type::UNIFORM, 1>(),
    std::forward<Args>(Arguments)...)) {

    switch (Type_) {
    case geometry_type::UNIFORM:
      switch (NumDims_) {
      case 1:
        return std::forward<F>(Func)(manipulator_type<geometry_type::UNIFORM, 1>(),
          std::forward<Args>(Arguments)...);
      case 2:
        return std::forward<F>(Func)(manipulator_type<geometry_type::UNIFORM, 2>(),
          std::forward<Args>(Arguments)...);
      default:
        return std::forward<F>(Func)(manipulator_type<geometry_type::UNIFORM, 3>(),
          std::forward<Args>(Arguments)...);
      }
    case geometry_type::RECTILINEAR:
      switch (NumDims_) {
      case 1:
        return std::forward<F>(Func)(manipulator_type<geometry_type::RECTILINEAR, 1>(),
          std::forward<Args>(Arguments)...);
      case 2:
        return std::forward<F>(Func)(manipulator_type<geometry_type::RECTILINEAR, 2>(),
          std::forward<Args>(Arguments)...);
      default:
        return std::forward<F>(Func)(manipulator_type<geometry_type::RECTILINEAR, 3>(),
          std::forward<Args>(Arguments)...);
      }
    case geometry_type::ORIENTED_UNIFORM:
      switch (NumDims_) {
      case 1:
        return std::forward<F>(Func)(manipulator_type<geometry_type::ORIENTED_UNIFORM, 1>(),
          std::forward<Args>(Arguments)...);
      case 2:
        return std::forward<F>(Func)(manipulator_type<geometry_type::ORIENTED_UNIFORM, 2>(),
          std::forward<Args>(Arguments)...);
      default:
        return std::forward<F>(Func)(manipulator_type<geometry_type::ORIENTED_UNIFORM, 3>(),
          std::forward<Args>(Arguments)...);
      }
    case geometry_type::ORIENTED_RECTILINEAR:
      switch (NumDims_) {
      case 1:
        return std::forward<F>(Func)(manipulator_type<geometry_type::ORIENTED_RECTILINEAR, 1>(),
          std::forward<Args>(Arguments)...);
      case 2:
        return std::forward<F>(Func)(manipulator_type<geometry_type::ORIENTED_RECTILINEAR, 2>(),
          std::forward<Args>(Arguments)...);
      default:
        return std::forward<F>(Func)(manipulator_type<geometry_type::ORIENTED_RECTILINEAR, 3>(),
          std::forward<Args>(Arguments)...);
      }
    default: // geometry_type::CURVILINEAR
      OVK_DEBUG_ASSERT(Type_ == geometry_type::CURVILINEAR, "Unhandled enum value.");
      switch (NumDims_) {
      case 1:
        return std::forward<F>(Func)(manipulator_type<geometry_type::CURVILINEAR, 1>(),
          std::forward<Args>(Arguments)...);
      case 2:
        return std::forward<F>(Func)(manipulator_type<geometry_type::CURVILINEAR, 2>(),
          std::forward<Args>(Arguments)...);
      default:
        return std::forward<F>(Func)(manipulator_type<geometry_type::CURVILINEAR, 3>(),
          std::forward<Args>(Arguments)...);
      }
    }

  }

  geometry_type Type() const { return Type_; }
  int Dimension() const { return NumDims_; }

private:

  geometry_type Type_;
  int NumDims_;

};

}}

#endif
