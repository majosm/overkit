// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_GEOMETRIC_PRIMITIVE_OPS_HPP_INCLUDED
#define OVK_CORE_GEOMETRIC_PRIMITIVE_OPS_HPP_INCLUDED

#include <ovk/core/ArrayView.hpp>
#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Math.hpp>

#include <cmath>
#include <cstdlib>

namespace ovk {
namespace core {

double IsoLine2Node(double LowerNodeCoord, double UpperNodeCoord, double LocalCoord);
double IsoLine2NodeInverse(double LowerNodeCoord, double UpperNodeCoord, double Coord);

double IsoLine4Node(const array_view<const double> &NodeCoords, double LocalCoord);
double IsoLine4NodeInverse(const array_view<const double> &NodeCoords, double Coord, bool
  *MaybeSuccess=nullptr, double Tolerance=1.e-12, int MaxSteps=100);

elem<double,2> IsoQuad4NodeUniform(const elem<double,2> &LowerNodeCoords, const elem<double,2>
  &UpperNodeCoords, const elem<double,2> &LocalCoords);
elem<double,2> IsoQuad4NodeUniformInverse(const elem<double,2> &LowerNodeCoords, const
  elem<double,2> &UpperNodeCoords, const elem<double,2> &Coords);

elem<double,2> IsoQuad4NodeOrientedUniform(const array_view<const elem<double,2>> &NodeCoords,
  const elem<double,2> &LocalCoords);
elem<double,2> IsoQuad4NodeOrientedUniformInverse(const array_view<const elem<double,2>>
  &NodeCoords, const elem<double,2> &Coords);

elem<double,2> IsoQuad4NodeNonUniform(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,2> &LocalCoords);
elem<double,2> IsoQuad4NodeNonUniform(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,2> &ShapeI, const elem<double,2> &ShapeJ);
elem<double,2> IsoQuad4NodeNonUniformInverse(const array_view<const elem<double,2>> &NodeCoords,
  const elem<double,2> &Coords, bool *MaybeSuccess=nullptr, double Tolerance=1.e-12, int
  MaxSteps=100);

elem<double,2> IsoQuad16Node(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,2> &LocalCoords);
elem<double,2> IsoQuad16Node(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,4> &ShapeI, const elem<double,4> &ShapeJ);
elem<double,2> IsoQuad16NodeInverse(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,2> &Coords, bool *MaybeSuccess=nullptr, double Tolerance=1.e-12, int MaxSteps=100);

elem<double,3> IsoHex8NodeUniform(const elem<double,3> &LowerNodeCoords, const elem<double,3>
  &UpperNodeCoords, const elem<double,3> &LocalCoords);
elem<double,3> IsoHex8NodeUniformInverse(const elem<double,3> &LowerNodeCoords, const elem<double,3>
  &UpperNodeCoords, const elem<double,3> &Coords);

elem<double,3> IsoHex8NodeOrientedUniform(const array_view<const elem<double,3>> &NodeCoords, const
  elem<double,3> &LocalCoords);
elem<double,3> IsoHex8NodeOrientedUniformInverse(const array_view<const elem<double,3>> &NodeCoords,
  const elem<double,3> &Coords);

elem<double,3> IsoHex8NodeNonUniform(const array_view<const elem<double,3>> &NodeCoords, const
  elem<double,3> &LocalCoords);
elem<double,3> IsoHex8NodeNonUniform(const array_view<const elem<double,3>> &NodeCoords, const
  elem<double,2> &ShapeI, const elem<double,2> &ShapeJ, const elem<double,2> &ShapeK);
elem<double,3> IsoHex8NodeNonUniformInverse(const array_view<const elem<double,3>> &NodeCoords,
  const elem<double,3> &Coords, bool *MaybeSuccess=nullptr, double Tolerance=1.e-12, int
  MaxSteps=100);

elem<double,3> IsoHex64Node(const array_view<const elem<double,3>> &NodeCoords, const elem<double,3>
  &LocalCoords);
elem<double,3> IsoHex64Node(const array_view<const elem<double,3>> &NodeCoords, const elem<double,4>
  &ShapeI, const elem<double,4> &ShapeJ, const elem<double,4> &ShapeK);
elem<double,3> IsoHex64NodeInverse(const array_view<const elem<double,3>> &NodeCoords, const
  elem<double,3> &Coords, bool *MaybeSuccess=nullptr, double Tolerance=1.e-12, int MaxSteps=100);

bool OverlapsLine(double LowerNodeCoord, double UpperNodeCoord, double Coords, double
  Tolerance=1.e-12);
bool OverlapsQuadUniform(const elem<double,2> &LowerNodeCoords, const elem<double,2>
  &UpperNodeCoords, const elem<double,2> &Coords, double Tolerance=1.e-12);
bool OverlapsQuadOrientedUniform(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,2> &Coords, double Tolerance=1.e-12);
bool OverlapsQuadNonUniform(const array_view<const elem<double,2>> &NodeCoords, const elem<double,2>
  &Coords, double Tolerance=1.e-12);
bool OverlapsHexUniform(const elem<double,3> &LowerNodeCoords, const elem<double,3>
  &UpperNodeCoords, const elem<double,3> &Coords, double Tolerance=1.e-12);
bool OverlapsHexOrientedUniform(const array_view<const elem<double,3>> &NodeCoords, const
  elem<double,3> &Coords, double Tolerance=1.e-12);
bool OverlapsHexNonUniform(const array_view<const elem<double,3>> &NodeCoords, const elem<double,3>
  &Coords, double Tolerance=1.e-12);

double VolumeLine(double LowerNodeCoord, double UpperNodeCoord);
double VolumeQuadUniform(const elem<double,2> &LowerNodeCoords, const elem<double,2>
  &UpperNodeCoords);
double VolumeQuadOrientedUniform(const array_view<const elem<double,2>> &NodeCoords);
double VolumeQuadNonUniform(const array_view<const elem<double,2>> &NodeCoords);
double VolumeHexUniform(const elem<double,3> &LowerNodeCoords, const elem<double,3>
  &UpperNodeCoords);
double VolumeHexOrientedUniform(const array_view<const elem<double,3>> &NodeCoords);
double VolumeHexNonUniform(const array_view<const elem<double,3>> &NodeCoords);

}}

#include <ovk/core/GeometricPrimitiveOpsLine.inl>
#include <ovk/core/GeometricPrimitiveOpsQuad.inl>
#include <ovk/core/GeometricPrimitiveOpsHex.inl>

#endif
