// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MATH_HPP_INCLUDED
#define OVK_CORE_MATH_HPP_INCLUDED

#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Tuple.hpp>

#include <cmath>

namespace ovk {
namespace core {

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

}}

#include <ovk/core/Math.inl>

#endif
