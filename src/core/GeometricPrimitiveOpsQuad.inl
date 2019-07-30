// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

inline elem<double,2> IsoQuad4NodeUniform(const elem<double,2> &LowerNodeCoords, const
  elem<double,2> &UpperNodeCoords, const elem<double,2> &LocalCoords) {

  double U = LocalCoords(0);
  double V = LocalCoords(1);

  return {
    (1.-U) * LowerNodeCoords(0) + U * UpperNodeCoords(0),
    (1.-V) * LowerNodeCoords(1) + V * UpperNodeCoords(1)
  };

}

inline elem<double,2> IsoQuad4NodeUniformInverse(const elem<double,2> &LowerNodeCoords, const
  elem<double,2> &UpperNodeCoords, const elem<double,2> &Coords) {

  double U = (Coords(0) - LowerNodeCoords(0))/(UpperNodeCoords(0) - LowerNodeCoords(0));
  double V = (Coords(1) - LowerNodeCoords(1))/(UpperNodeCoords(1) - LowerNodeCoords(1));

  return {U, V};

}

inline elem<double,2> IsoQuad4NodeOrientedUniform(const array_view<const elem<double,2>>
  &NodeCoords, const elem<double,2> &LocalCoords) {

  double U = LocalCoords(0);
  double V = LocalCoords(1);

  return {
    (1.-U-V)*NodeCoords(0)(0) + U*NodeCoords(1)(0) + V*NodeCoords(2)(0),
    (1.-U-V)*NodeCoords(0)(1) + U*NodeCoords(1)(1) + V*NodeCoords(2)(1),
  };

}

inline elem<double,2> IsoQuad4NodeOrientedUniformInverse(const array_view<const elem<double,2>>
  &NodeCoords, const elem<double,2> &Coords) {

  elem<double,2> I = {
    NodeCoords(1)(0) - NodeCoords(0)(0),
    NodeCoords(1)(1) - NodeCoords(0)(1),
  };

  elem<double,2> J = {
    NodeCoords(2)(0) - NodeCoords(0)(0),
    NodeCoords(2)(1) - NodeCoords(0)(1),
  };

  elem<double,2> RelativeCoords = {
    Coords(0) - NodeCoords(0)(0),
    Coords(1) - NodeCoords(0)(1),
  };

  auto Dot = [](const elem<double,2> &Vec1, const elem<double,2> &Vec2) -> double {
    return Vec1(0)*Vec2(0) + Vec1(1)*Vec2(1);
  };
  auto Project = [Dot](const elem<double,2> &Vec1, const elem<double,2> &Vec2) -> double {
    return Dot(Vec1, Vec2)/Dot(Vec2, Vec2);
  };

  double U = Project(RelativeCoords, I);
  double V = Project(RelativeCoords, J);

  return {U, V};

}

inline elem<double,2> IsoQuad4NodeNonUniform(const array_view<const elem<double,2>> &NodeCoords,
  const elem<double,2> &LocalCoords) {

  elem<double,2> ShapeI = LagrangeInterpLinear(LocalCoords(0));
  elem<double,2> ShapeJ = LagrangeInterpLinear(LocalCoords(1));

  return IsoQuad4NodeNonUniform(NodeCoords, ShapeI, ShapeJ);

}

inline elem<double,2> IsoQuad4NodeNonUniform(const array_view<const elem<double,2>> &NodeCoords,
  const elem<double,2> &ShapeI, const elem<double,2> &ShapeJ) {

  double Coefs[] = {
    ShapeI(0)*ShapeJ(0),
    ShapeI(1)*ShapeJ(0),
    ShapeI(0)*ShapeJ(1),
    ShapeI(1)*ShapeJ(1)
  };

  elem<double,2> Coords;

  for (int iDim = 0; iDim < 2; ++iDim) {
    Coords(iDim) =
      Coefs[0] * NodeCoords(0)(iDim) + Coefs[1] * NodeCoords(1)(iDim) +
      Coefs[2] * NodeCoords(2)(iDim) + Coefs[3] * NodeCoords(3)(iDim);
  }

  return Coords;

}

inline elem<double,2> IsoQuad4NodeNonUniformInverse(const array_view<const elem<double,2>>
  &NodeCoords, const elem<double,2> &Coords, bool *MaybeSuccess, double Tolerance, int MaxSteps) {

  auto SmallEnough = [Tolerance](const elem<double,2> &Tuple) -> bool {
    return
      std::abs(Tuple(0)) <= Tolerance &&
      std::abs(Tuple(1)) <= Tolerance;
  };

  auto IsNaN = [](const elem<double,2> &Tuple) -> bool {
    return
      std::isnan(Tuple(0)) ||
      std::isnan(Tuple(1));
  };

  elem<double,2> LocalCoords = {0.5, 0.5};

  for (int iStep = 0; iStep < MaxSteps; ++iStep) {
    elem<double,2> ShapeI = LagrangeInterpLinear(LocalCoords(0));
    elem<double,2> ShapeJ = LagrangeInterpLinear(LocalCoords(1));
    elem<double,2> Error = Coords - IsoQuad4NodeNonUniform(NodeCoords, ShapeI, ShapeJ);
    if (SmallEnough(Error)) break;
    elem<double,2> InterpDerivI = LagrangeInterpLinearDeriv(LocalCoords(0));
    elem<double,2> InterpDerivJ = LagrangeInterpLinearDeriv(LocalCoords(1));
    elem<double,2> JacobianI = {0., 0.};
    elem<double,2> JacobianJ = {0., 0.};
    int iNode = 0;
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        for (int iDim = 0; iDim < 2; ++iDim) {
          double Coord = NodeCoords(iNode)(iDim);
          JacobianI(iDim) += InterpDerivI(i) * ShapeJ(j) * Coord;
          JacobianJ(iDim) += ShapeI(i) * InterpDerivJ(j) * Coord;
        }
        ++iNode;
      }
    }
    LocalCoords = LocalCoords + ColumnSolve2D(JacobianI, JacobianJ, Error);
  }

  if (MaybeSuccess) {
    bool &Success = *MaybeSuccess;
    elem<double,2> ShapeI = LagrangeInterpLinear(LocalCoords(0));
    elem<double,2> ShapeJ = LagrangeInterpLinear(LocalCoords(1));
    elem<double,2> Error = Coords - IsoQuad4NodeNonUniform(NodeCoords, ShapeI, ShapeJ);
    Success = SmallEnough(Error) && !IsNaN(Error);
  }

  return LocalCoords;

}

inline elem<double,2> IsoQuad16Node(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,2> &LocalCoords) {

  elem<double,4> ShapeI = LagrangeInterpCubic(LocalCoords(0));
  elem<double,4> ShapeJ = LagrangeInterpCubic(LocalCoords(1));

  return IsoQuad16Node(NodeCoords, ShapeI, ShapeJ);

}

inline elem<double,2> IsoQuad16Node(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,4> &ShapeI, const elem<double,4> &ShapeJ) {

  elem<double,2> Coords = {0., 0.};

  int iNode = 0;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      for (int iDim = 0; iDim < 2; ++iDim) {
        Coords(iDim) += ShapeI(i) * ShapeJ(j) * NodeCoords(iNode)(iDim);
      }
      ++iNode;
    }
  }

  return Coords;

}

inline elem<double,2> IsoQuad16NodeInverse(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,2> &Coords, bool *MaybeSuccess, double Tolerance, int MaxSteps) {

  auto SmallEnough = [Tolerance](const elem<double,2> &Tuple) -> bool {
    return
      std::abs(Tuple(0)) <= Tolerance &&
      std::abs(Tuple(1)) <= Tolerance;
  };

  auto IsNaN = [](const elem<double,2> &Tuple) -> bool {
    return
      std::isnan(Tuple(0)) ||
      std::isnan(Tuple(1));
  };

  elem<double,2> LocalCoords = {0.5, 0.5};

  for (int iStep = 0; iStep < MaxSteps; ++iStep) {
    elem<double,4> ShapeI = LagrangeInterpCubic(LocalCoords(0));
    elem<double,4> ShapeJ = LagrangeInterpCubic(LocalCoords(1));
    elem<double,2> Error = Coords - IsoQuad16Node(NodeCoords, ShapeI, ShapeJ);
    if (SmallEnough(Error)) break;
    elem<double,4> InterpDerivI = LagrangeInterpCubicDeriv(LocalCoords(0));
    elem<double,4> InterpDerivJ = LagrangeInterpCubicDeriv(LocalCoords(1));
    elem<double,2> JacobianI = {0., 0.};
    elem<double,2> JacobianJ = {0., 0.};
    int iNode = 0;
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 4; ++i) {
        for (int iDim = 0; iDim < 2; ++iDim) {
          double Coord = NodeCoords(iNode)(iDim);
          JacobianI(iDim) += InterpDerivI(i) * ShapeJ(j) * Coord;
          JacobianJ(iDim) += ShapeI(i) * InterpDerivJ(j) * Coord;
        }
        ++iNode;
      }
    }
    LocalCoords = LocalCoords + ColumnSolve2D(JacobianI, JacobianJ, Error);
  }

  if (MaybeSuccess) {
    bool &Success = *MaybeSuccess;
    elem<double,4> ShapeI = LagrangeInterpCubic(LocalCoords(0));
    elem<double,4> ShapeJ = LagrangeInterpCubic(LocalCoords(1));
    elem<double,2> Error = Coords - IsoQuad16Node(NodeCoords, ShapeI, ShapeJ);
    Success = SmallEnough(Error) && !IsNaN(Error);
  }

  return LocalCoords;

}

inline bool OverlapsQuadUniform(const elem<double,2> &LowerNodeCoords, const elem<double,2>
  &UpperNodeCoords, const elem<double,2> &Coords, double Tolerance) {

  elem<double,2> LocalCoords = IsoQuad4NodeUniformInverse(LowerNodeCoords, UpperNodeCoords, Coords);

  for (int iDim = 0; iDim < 2; ++iDim) {
    if (LocalCoords(iDim) < -Tolerance) {
      return false;
    } else if (LocalCoords(iDim) > 1.+Tolerance) {
      return false;
    }
  }

  return true;

}

inline bool OverlapsQuadOrientedUniform(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,2> &Coords, double Tolerance) {

  elem<double,2> LocalCoords = IsoQuad4NodeOrientedUniformInverse(NodeCoords, Coords);

  for (int iDim = 0; iDim < 2; ++iDim) {
    if (LocalCoords(iDim) < -Tolerance) {
      return false;
    } else if (LocalCoords(iDim) > 1.+Tolerance) {
      return false;
    }
  }

  return true;

}

inline bool OverlapsQuadNonUniform(const array_view<const elem<double,2>> &NodeCoords, const
  elem<double,2> &Coords, double Tolerance) {

  // Not safe to invert isoparametric mapping; instead decompose into 2 triangles
  constexpr int Triangles[2][3] = {
    {0,1,2},
    {3,2,1}
  };

  for (int iTriangle = 0; iTriangle < 2; ++iTriangle) {
    int iVertex1 = Triangles[iTriangle][0];
    int iVertex2 = Triangles[iTriangle][1];
    int iVertex3 = Triangles[iTriangle][2];
    elem<double,2> I = {
      NodeCoords(iVertex2)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex2)(1) - NodeCoords(iVertex1)(1)
    };
    elem<double,2> J = {
      NodeCoords(iVertex3)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex3)(1) - NodeCoords(iVertex1)(1)
    };
    elem<double,2> RelativeCoords = {
      Coords(0) - NodeCoords(iVertex1)(0),
      Coords(1) - NodeCoords(iVertex1)(1)
    };
    elem<double,2> LocalCoords = ColumnSolve2D(I, J, RelativeCoords);
    bool Inside = true;
    double CoordSum = 0.;
    for (int iDim = 0; iDim < 2; ++iDim) {
      CoordSum += LocalCoords(iDim);
      if (LocalCoords(iDim) < -Tolerance || CoordSum > 1.+Tolerance) {
        Inside = false;
        break;
      }
    }
    if (Inside) return true;
  }

  return false;

}

inline double VolumeQuadUniform(const elem<double,2> &LowerNodeCoords, const elem<double,2>
  &UpperNodeCoords) {

  return
    (UpperNodeCoords(0) - LowerNodeCoords(0))*
    (UpperNodeCoords(1) - LowerNodeCoords(1));

}

inline double VolumeQuadOrientedUniform(const array_view<const elem<double,2>> &NodeCoords) {

  elem<double,2> I = {
    NodeCoords(1)(0) - NodeCoords(0)(0),
    NodeCoords(1)(1) - NodeCoords(0)(1)
  };

  elem<double,2> J = {
    NodeCoords(2)(0) - NodeCoords(0)(0),
    NodeCoords(2)(1) - NodeCoords(0)(1)
  };

  return std::abs(ColumnDeterminant2D(I, J));

}

inline double VolumeQuadNonUniform(const array_view<const elem<double,2>> &NodeCoords) {

  // Decompose into 2 triangles
  constexpr int Triangles[2][3] = {
    {0,1,2},
    {3,2,1}
  };

  double Volume = 0.;

  for (int iTriangle = 0; iTriangle < 2; ++iTriangle) {
    int iVertex1 = Triangles[iTriangle][0];
    int iVertex2 = Triangles[iTriangle][1];
    int iVertex3 = Triangles[iTriangle][2];
    elem<double,2> I = {
      NodeCoords(iVertex2)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex2)(1) - NodeCoords(iVertex1)(1)
    };
    elem<double,2> J = {
      NodeCoords(iVertex3)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex3)(1) - NodeCoords(iVertex1)(1)
    };
    Volume += 0.5 * std::abs(ColumnDeterminant2D(I, J));
  }

  return Volume;

}

}}
