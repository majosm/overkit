// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

inline double IsoLine2Node(double LowerNodeCoord, double UpperNodeCoord, double LocalCoord) {

  double U = LocalCoord;

  return (1.-U) * LowerNodeCoord + U * UpperNodeCoord;

}

inline double IsoLine2NodeInverse(double LowerNodeCoord, double UpperNodeCoord, double Coord) {

  double U = (Coord - LowerNodeCoord)/(UpperNodeCoord - LowerNodeCoord);

  return U;

}

inline double IsoLine4Node(const array_view<const double> &NodeCoords, double LocalCoord) {

  elem<double,4> Interp = LagrangeInterpCubic(LocalCoord);

  double Coord = 0.;

  for (int i = 0; i < 4; ++i) {
    Coord += Interp(i) * NodeCoords(i);
  }

  return Coord;

}

inline optional<double> IsoLine4NodeInverse(const array_view<const double> &NodeCoords, double
  Coord, double Tolerance, int MaxSteps) {

  double LocalCoord = 0.5;

  for (int iStep = 0; iStep < MaxSteps; ++iStep) {
    double Error = Coord - IsoLine4Node(NodeCoords, LocalCoord);
    if (std::abs(Error) <= Tolerance) break;
    elem<double,4> InterpDeriv = LagrangeInterpCubicDeriv(LocalCoord);
    double Deriv = 0.;
    for (int i = 0; i < 4; ++i) {
      Deriv += InterpDeriv(i) * NodeCoords(i);
    }
    LocalCoord = LocalCoord + Error/Deriv;
  }

  double Error = Coord - IsoLine4Node(NodeCoords, LocalCoord);
  if (Error <= Tolerance && !IsNaN(Error)) {
    return LocalCoord;
  } else {
    return {};
  }

}

inline bool OverlapsLine(double LowerNodeCoord, double UpperNodeCoord, double Coord, double
  Tolerance) {

  double LocalCoord = IsoLine2NodeInverse(LowerNodeCoord, UpperNodeCoord, Coord);

  if (LocalCoord < -Tolerance) {
    return false;
  } else if (LocalCoord > 1.+Tolerance) {
    return false;
  }

  return true;

}

inline double VolumeLine(double LowerNodeCoord, double UpperNodeCoord) {

  return UpperNodeCoord - LowerNodeCoord;

}

}}
