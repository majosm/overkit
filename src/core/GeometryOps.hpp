// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GEOMETRY_OPS_HPP_INCLUDED
#define OVK_CORE_GEOMETRY_OPS_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Box.hpp>
#include <ovk/core/Debug.hpp>
#include <ovk/core/Field.hpp>
#include <ovk/core/GeometricPrimitiveOps.hpp>
#include <ovk/core/GeometryBase.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Range.hpp>
#include <ovk/core/ScalarOps.hpp>
#include <ovk/core/Tuple.hpp>

namespace ovk {
namespace core {

template <geometry_type Type, int NumDims> bool OverlapsCell(const array_view<const field_view<const
  double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double> &PointCoords);

template <geometry_type Type, int NumDims> optional<tuple<double>> CoordsInCell(const array_view<
  const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords);

template <geometry_type Type, int NumDims> double CellVolume(const array_view<const field_view<const
  double>> &Coords, const tuple<int> &Cell);

template <geometry_type Type, int NumDims> box CellBounds(const array_view<const field_view<const
  double>> &Coords, const tuple<int> &Cell);

}}

#include <ovk/core/GeometryOps.inl>

#endif
