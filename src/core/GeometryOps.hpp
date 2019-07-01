// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GEOMETRY_OPS_HPP_INCLUDED
#define OVK_CORE_GEOMETRY_OPS_HPP_INCLUDED

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Tuple.hpp>

#include <cmath>

namespace ovk {
namespace core {

template <typename IndexType=long long> IndexType CartesianGridCell1D(double Origin, double
  CellSize, double Coord);
template <typename IndexType=long long> tuple<IndexType> CartesianGridCell2D(const tuple<double>
  &Origin, const tuple<double> &CellSize, const tuple<double> &Coords);
template <typename IndexType=long long> tuple<IndexType> CartesianGridCell3D(const tuple<double>
  &Origin, const tuple<double> &CellSize, const tuple<double> &Coords);

double ColumnDeterminant2D(const tuple<double> &AI, const tuple<double> &AJ);
double ColumnDeterminant3D(const tuple<double> &AI, const tuple<double> &AJ, const tuple<double>
  &AK);

tuple<double> ColumnSolve2D(const tuple<double> &AI, const tuple<double> &AJ, const tuple<double>
  &B);
tuple<double> ColumnSolve3D(const tuple<double> &AI, const tuple<double> &AJ, const tuple<double>
  &AK, const tuple<double> &B);

elem<double,2> LagrangeInterpLinear(double U);
elem<double,2> LagrangeInterpLinearDeriv(double U);
elem<double,4> LagrangeInterpCubic(double U);
elem<double,4> LagrangeInterpCubicDeriv(double U);

double IsoLine2Node(double LowerNodeCoord, double UpperNodeCoord, double LocalCoord);
double IsoLine2NodeInverse(double LowerNodeCoord, double UpperNodeCoord, double Coord);

double IsoLine4Node(const array_view<const double> &NodeCoords, double LocalCoord);
double IsoLine4NodeInverse(const array_view<const double> &NodeCoords, double Coord, bool
  *MaybeSuccess=nullptr, double Tolerance=1.e-12, int MaxSteps=100);

tuple<double> IsoQuad4NodeUniform(const tuple<double> &LowerNodeCoords, const tuple<double>
  &UpperNodeCoords, const tuple<double> &LocalCoords);
tuple<double> IsoQuad4NodeUniformInverse(const tuple<double> &LowerNodeCoords, const tuple<double>
  &UpperNodeCoords, const tuple<double> &Coords);

tuple<double> IsoQuad4NodeOrientedUniform(const array_view<const tuple<double>> &NodeCoords,
  const tuple<double> &LocalCoords);
tuple<double> IsoQuad4NodeOrientedUniformInverse(const array_view<const tuple<double>> &NodeCoords,
  const tuple<double> &Coords);

tuple<double> IsoQuad4NodeNonUniform(const array_view<const tuple<double>> &NodeCoords, const
  tuple<double> &LocalCoords);
tuple<double> IsoQuad4NodeNonUniform(const array_view<const tuple<double>> &NodeCoords, const
  elem<double,2> &ShapeI, const elem<double,2> &ShapeJ);
tuple<double> IsoQuad4NodeNonUniformInverse(const array_view<const tuple<double>> &NodeCoords,
  const tuple<double> &Coords, bool *MaybeSuccess=nullptr, double Tolerance=1.e-12, int
  MaxSteps=100);

tuple<double> IsoQuad16Node(const array_view<const tuple<double>> &NodeCoords, const tuple<double>
  &LocalCoords);
tuple<double> IsoQuad16Node(const array_view<const tuple<double>> &NodeCoords, const elem<double,4>
  &ShapeI, const elem<double,4> &ShapeJ);
tuple<double> IsoQuad16NodeInverse(const array_view<const tuple<double>> &NodeCoords, const
  tuple<double> &Coords, bool *MaybeSuccess=nullptr, double Tolerance=1.e-12, int MaxSteps=100);

tuple<double> IsoHex8NodeUniform(const tuple<double> &LowerNodeCoords, const tuple<double>
  &UpperNodeCoords, const tuple<double> &LocalCoords);
tuple<double> IsoHex8NodeUniformInverse(const tuple<double> &LowerNodeCoords, const tuple<double>
  &UpperNodeCoords, const tuple<double> &Coords);

tuple<double> IsoHex8NodeOrientedUniform(const array_view<const tuple<double>> &NodeCoords, const
  tuple<double> &LocalCoords);
tuple<double> IsoHex8NodeOrientedUniformInverse(const array_view<const tuple<double>> &NodeCoords,
  const tuple<double> &Coords);

tuple<double> IsoHex8NodeNonUniform(const array_view<const tuple<double>> &NodeCoords, const
  tuple<double> &LocalCoords);
tuple<double> IsoHex8NodeNonUniform(const array_view<const tuple<double>> &NodeCoords, const
  elem<double,2> &ShapeI, const elem<double,2> &ShapeJ, const elem<double,2> &ShapeK);
tuple<double> IsoHex8NodeNonUniformInverse(const array_view<const tuple<double>> &NodeCoords,
  const tuple<double> &Coords, bool *MaybeSuccess=nullptr, double Tolerance=1.e-12, int
  MaxSteps=100);

tuple<double> IsoHex64Node(const array_view<const tuple<double>> &NodeCoords, const tuple<double>
  &LocalCoords);
tuple<double> IsoHex64Node(const array_view<const tuple<double>> &NodeCoords, const elem<double,4>
  &ShapeI, const elem<double,4> &ShapeJ, const elem<double,4> &ShapeK);
tuple<double> IsoHex64NodeInverse(const array_view<const tuple<double>> &NodeCoords, const
  tuple<double> &Coords, bool *MaybeSuccess=nullptr, double Tolerance=1.e-12, int MaxSteps=100);

bool OverlapsLine(double LowerNodeCoord, double UpperNodeCoord, double Coords, double
  Tolerance=1.e-12);
bool OverlapsQuadUniform(const tuple<double> &LowerNodeCoords, const tuple<double> &UpperNodeCoords,
  const tuple<double> &Coords, double Tolerance=1.e-12);
bool OverlapsQuadOrientedUniform(const array_view<const tuple<double>> &NodeCoords, const
  tuple<double> &Coords, double Tolerance=1.e-12);
bool OverlapsQuadNonUniform(const array_view<const tuple<double>> &NodeCoords, const tuple<double>
  &Coords, double Tolerance=1.e-12);
bool OverlapsHexUniform(const tuple<double> &LowerNodeCoords, const tuple<double> &UpperNodeCoords,
  const tuple<double> &Coords, double Tolerance=1.e-12);
bool OverlapsHexOrientedUniform(const array_view<const tuple<double>> &NodeCoords, const
  tuple<double> &Coords, double Tolerance=1.e-12);
bool OverlapsHexNonUniform(const array_view<const tuple<double>> &NodeCoords, const tuple<double>
  &Coords, double Tolerance=1.e-12);

double VolumeLine(double LowerNodeCoord, double UpperNodeCoord);
double VolumeQuadUniform(const tuple<double> &LowerNodeCoords, const tuple<double> &UpperNodeCoords);
double VolumeQuadOrientedUniform(const array_view<const tuple<double>> &NodeCoords);
double VolumeQuadNonUniform(const array_view<const tuple<double>> &NodeCoords);
double VolumeHexUniform(const tuple<double> &LowerNodeCoords, const tuple<double> &UpperNodeCoords);
double VolumeHexOrientedUniform(const array_view<const tuple<double>> &NodeCoords);
double VolumeHexNonUniform(const array_view<const tuple<double>> &NodeCoords);

}}

#include <ovk/core/GeometryOps.inl>
#include <ovk/core/GeometryOpsLine.inl>
#include <ovk/core/GeometryOpsQuad.inl>
#include <ovk/core/GeometryOpsHex.inl>

#endif
