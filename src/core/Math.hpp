// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_MATH_HPP_INCLUDED
#define OVK_CORE_MATH_HPP_INCLUDED

#include <ovk/core/Elem.hpp>
#include <ovk/core/Global.hpp>

#include <cmath>

namespace ovk {
namespace core {

double ColumnDeterminant2D(const elem<double,2> &AI, const elem<double,2> &AJ);
double ColumnDeterminant3D(const elem<double,3> &AI, const elem<double,3> &AJ, const elem<double,3>
  &AK);

elem<double,2> ColumnSolve2D(const elem<double,2> &AI, const elem<double,2> &AJ, const
  elem<double,2> &B);
elem<double,3> ColumnSolve3D(const elem<double,3> &AI, const elem<double,3> &AJ, const
  elem<double,3> &AK, const elem<double,3> &B);

elem<double,2> LagrangeInterpLinear(double U);
elem<double,2> LagrangeInterpLinearDeriv(double U);
elem<double,4> LagrangeInterpCubic(double U);
elem<double,4> LagrangeInterpCubicDeriv(double U);

}}

#include <ovk/core/Math.inl>

#endif
