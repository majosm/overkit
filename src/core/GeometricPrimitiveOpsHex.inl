// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

inline elem<double,3> IsoHex8NodeUniform(const elem<double,3> &LowerNodeCoords, const elem<double,3>
  &UpperNodeCoords, const elem<double,3> &LocalCoords) {

  double U = LocalCoords(0);
  double V = LocalCoords(1);
  double W = LocalCoords(2);

  return {
    (1.-U) * LowerNodeCoords(0) + U * UpperNodeCoords(0),
    (1.-V) * LowerNodeCoords(1) + V * UpperNodeCoords(1),
    (1.-W) * LowerNodeCoords(2) + W * UpperNodeCoords(2)
  };

}

inline elem<double,3> IsoHex8NodeUniformInverse(const elem<double,3> &LowerNodeCoords, const
  elem<double,3> &UpperNodeCoords, const elem<double,3> &Coords) {

  double U = (Coords(0) - LowerNodeCoords(0))/(UpperNodeCoords(0) - LowerNodeCoords(0));
  double V = (Coords(1) - LowerNodeCoords(1))/(UpperNodeCoords(1) - LowerNodeCoords(1));
  double W = (Coords(2) - LowerNodeCoords(2))/(UpperNodeCoords(2) - LowerNodeCoords(2));

  return {U, V, W};

}

inline elem<double,3> IsoHex8NodeOrientedUniform(const array_view<const elem<double,3>> &NodeCoords,
  const elem<double,3> &LocalCoords) {

  double U = LocalCoords(0);
  double V = LocalCoords(1);
  double W = LocalCoords(2);

  return {
    (1.-U-V-W)*NodeCoords(0)(0) + U*NodeCoords(1)(0) + V*NodeCoords(2)(0) + W*NodeCoords(4)(0),
    (1.-U-V-W)*NodeCoords(0)(1) + U*NodeCoords(1)(1) + V*NodeCoords(2)(1) + W*NodeCoords(4)(1),
    (1.-U-V-W)*NodeCoords(0)(2) + U*NodeCoords(1)(2) + V*NodeCoords(2)(2) + W*NodeCoords(4)(2)
  };

}

inline elem<double,3> IsoHex8NodeOrientedUniformInverse(const array_view<const elem<double,3>>
  &NodeCoords, const elem<double,3> &Coords) {

  elem<double,3> I = {
    NodeCoords(1)(0) - NodeCoords(0)(0),
    NodeCoords(1)(1) - NodeCoords(0)(1),
    NodeCoords(1)(2) - NodeCoords(0)(2)
  };

  elem<double,3> J = {
    NodeCoords(2)(0) - NodeCoords(0)(0),
    NodeCoords(2)(1) - NodeCoords(0)(1),
    NodeCoords(2)(2) - NodeCoords(0)(2)
  };

  elem<double,3> K = {
    NodeCoords(4)(0) - NodeCoords(0)(0),
    NodeCoords(4)(1) - NodeCoords(0)(1),
    NodeCoords(4)(2) - NodeCoords(0)(2)
  };

  elem<double,3> RelativeCoords = {
    Coords(0) - NodeCoords(0)(0),
    Coords(1) - NodeCoords(0)(1),
    Coords(2) - NodeCoords(0)(2)
  };

  auto Dot = [](const elem<double,3> &Vec1, const elem<double,3> &Vec2) -> double {
    return Vec1(0)*Vec2(0) + Vec1(1)*Vec2(1) + Vec1(2)*Vec2(2);
  };
  auto Project = [Dot](const elem<double,3> &Vec1, const elem<double,3> &Vec2) -> double {
    return Dot(Vec1, Vec2)/Dot(Vec2, Vec2);
  };

  double U = Project(RelativeCoords, I);
  double V = Project(RelativeCoords, J);
  double W = Project(RelativeCoords, K);

  return {U, V, W};

}

inline elem<double,3> IsoHex8NodeNonUniform(const array_view<const elem<double,3>> &NodeCoords,
  const elem<double,3> &LocalCoords) {

  elem<double,2> ShapeI = LagrangeInterpLinear(LocalCoords(0));
  elem<double,2> ShapeJ = LagrangeInterpLinear(LocalCoords(1));
  elem<double,2> ShapeK = LagrangeInterpLinear(LocalCoords(2));

  return IsoHex8NodeNonUniform(NodeCoords, ShapeI, ShapeJ, ShapeK);

}

inline elem<double,3> IsoHex8NodeNonUniform(const array_view<const elem<double,3>> &NodeCoords,
  const elem<double,2> &ShapeI, const elem<double,2> &ShapeJ, const elem<double,2> &ShapeK) {

  elem<double,3> Coords = {0., 0., 0.};

  int iNode = 0;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        for (int iDim = 0; iDim < 3; ++iDim) {
          Coords(iDim) += ShapeI(i) * ShapeJ(j) * ShapeK(k) * NodeCoords(iNode)(iDim);
        }
        ++iNode;
      }
    }
  }

  return Coords;

}

inline optional<elem<double,3>> IsoHex8NodeNonUniformInverse(const array_view<const elem<double,3>>
  &NodeCoords, const elem<double,3> &Coords, double Tolerance, int MaxSteps) {

  auto SmallEnough = [Tolerance](const elem<double,3> &Tuple) -> bool {
    return
      std::abs(Tuple(0)) <= Tolerance &&
      std::abs(Tuple(1)) <= Tolerance &&
      std::abs(Tuple(2)) <= Tolerance;
  };

  elem<double,3> LocalCoords = {0.5, 0.5, 0.5};

  for (int iStep = 0; iStep < MaxSteps; ++iStep) {
    elem<double,2> ShapeI = LagrangeInterpLinear(LocalCoords(0));
    elem<double,2> ShapeJ = LagrangeInterpLinear(LocalCoords(1));
    elem<double,2> ShapeK = LagrangeInterpLinear(LocalCoords(2));
    elem<double,3> Error = Coords - IsoHex8NodeNonUniform(NodeCoords, ShapeI, ShapeJ, ShapeK);
    if (SmallEnough(Error)) break;
    elem<double,2> InterpDerivI = LagrangeInterpLinearDeriv(LocalCoords(0));
    elem<double,2> InterpDerivJ = LagrangeInterpLinearDeriv(LocalCoords(1));
    elem<double,2> InterpDerivK = LagrangeInterpLinearDeriv(LocalCoords(2));
    elem<double,3> JacobianI = {0., 0., 0.};
    elem<double,3> JacobianJ = {0., 0., 0.};
    elem<double,3> JacobianK = {0., 0., 0.};
    int iNode = 0;
    for (int k = 0; k < 2; ++k) {
      for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 2; ++i) {
          for (int iDim = 0; iDim < 3; ++iDim) {
            double Coord = NodeCoords(iNode)(iDim);
            JacobianI(iDim) += InterpDerivI(i) * ShapeJ(j) * ShapeK(k) * Coord;
            JacobianJ(iDim) += ShapeI(i) * InterpDerivJ(j) * ShapeK(k) * Coord;
            JacobianK(iDim) += ShapeI(i) * ShapeJ(j) * InterpDerivK(k) * Coord;
          }
          ++iNode;
        }
      }
    }
    LocalCoords = LocalCoords + ColumnSolve3D(JacobianI, JacobianJ, JacobianK, Error);
  }

  elem<double,2> ShapeI = LagrangeInterpLinear(LocalCoords(0));
  elem<double,2> ShapeJ = LagrangeInterpLinear(LocalCoords(1));
  elem<double,2> ShapeK = LagrangeInterpLinear(LocalCoords(2));
  elem<double,3> Error = Coords - IsoHex8NodeNonUniform(NodeCoords, ShapeI, ShapeJ, ShapeK);
  if (SmallEnough(Error) && !IsNaN(Error)) {
    return LocalCoords;
  } else {
    return {};
  }

}

inline elem<double,3> IsoHex64Node(const array_view<const elem<double,3>> &NodeCoords, const
  elem<double,3> &LocalCoords) {

  elem<double,4> ShapeI = LagrangeInterpCubic(LocalCoords(0));
  elem<double,4> ShapeJ = LagrangeInterpCubic(LocalCoords(1));
  elem<double,4> ShapeK = LagrangeInterpCubic(LocalCoords(2));

  return IsoHex64Node(NodeCoords, ShapeI, ShapeJ, ShapeK);

}

inline elem<double,3> IsoHex64Node(const array_view<const elem<double,3>> &NodeCoords, const
  elem<double,4> &ShapeI, const elem<double,4> &ShapeJ, const elem<double,4> &ShapeK) {

  elem<double,3> Coords = {0., 0., 0.};

  int iNode = 0;
  for (int k = 0; k < 4; ++k) {
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 4; ++i) {
        for (int iDim = 0; iDim < 3; ++iDim) {
          Coords(iDim) += ShapeI(i) * ShapeJ(j) * ShapeK(k) * NodeCoords(iNode)(iDim);
        }
        ++iNode;
      }
    }
  }

  return Coords;

}

inline optional<elem<double,3>> IsoHex64NodeInverse(const array_view<const elem<double,3>>
  &NodeCoords, const elem<double,3> &Coords, double Tolerance, int MaxSteps) {

  auto SmallEnough = [Tolerance](const elem<double,3> &Tuple) -> bool {
    return
      std::abs(Tuple(0)) <= Tolerance &&
      std::abs(Tuple(1)) <= Tolerance &&
      std::abs(Tuple(2)) <= Tolerance;
  };

  elem<double,3> LocalCoords = {0.5, 0.5, 0.5};

  for (int iStep = 0; iStep < MaxSteps; ++iStep) {
    elem<double,4> ShapeI = LagrangeInterpCubic(LocalCoords(0));
    elem<double,4> ShapeJ = LagrangeInterpCubic(LocalCoords(1));
    elem<double,4> ShapeK = LagrangeInterpCubic(LocalCoords(2));
    elem<double,3> Error = Coords - IsoHex64Node(NodeCoords, ShapeI, ShapeJ, ShapeK);
    if (SmallEnough(Error)) break;
    elem<double,4> InterpDerivI = LagrangeInterpCubicDeriv(LocalCoords(0));
    elem<double,4> InterpDerivJ = LagrangeInterpCubicDeriv(LocalCoords(1));
    elem<double,4> InterpDerivK = LagrangeInterpCubicDeriv(LocalCoords(2));
    elem<double,3> JacobianI = {0., 0., 0.};
    elem<double,3> JacobianJ = {0., 0., 0.};
    elem<double,3> JacobianK = {0., 0., 0.};
    int iNode = 0;
    for (int k = 0; k < 4; ++k) {
      for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
          for (int iDim = 0; iDim < 3; ++iDim) {
            double Coord = NodeCoords(iNode)(iDim);
            JacobianI(iDim) += InterpDerivI(i) * ShapeJ(j) * ShapeK(k) * Coord;
            JacobianJ(iDim) += ShapeI(i) * InterpDerivJ(j) * ShapeK(k) * Coord;
            JacobianK(iDim) += ShapeI(i) * ShapeJ(j) * InterpDerivK(k) * Coord;
          }
          ++iNode;
        }
      }
    }
    LocalCoords = LocalCoords + ColumnSolve3D(JacobianI, JacobianJ, JacobianK, Error);
  }

  elem<double,4> ShapeI = LagrangeInterpCubic(LocalCoords(0));
  elem<double,4> ShapeJ = LagrangeInterpCubic(LocalCoords(1));
  elem<double,4> ShapeK = LagrangeInterpCubic(LocalCoords(2));
  elem<double,3> Error = Coords - IsoHex64Node(NodeCoords, ShapeI, ShapeJ, ShapeK);
  if (SmallEnough(Error) && !IsNaN(Error)) {
    return LocalCoords;
  } else {
    return {};
  }

}

inline bool OverlapsHexUniform(const elem<double,3> &LowerNodeCoords, const elem<double,3>
  &UpperNodeCoords, const elem<double,3> &Coords, double Tolerance) {

  elem<double,3> LocalCoords = IsoHex8NodeUniformInverse(LowerNodeCoords, UpperNodeCoords, Coords);

  for (int iDim = 0; iDim < 3; ++iDim) {
    if (LocalCoords(iDim) < -Tolerance) {
      return false;
    } else if (LocalCoords(iDim) > 1.+Tolerance) {
      return false;
    }
  }

  return true;

}

inline bool OverlapsHexOrientedUniform(const array_view<const elem<double,3>> &NodeCoords, const
  elem<double,3> &Coords, double Tolerance) {

  elem<double,3> LocalCoords = IsoHex8NodeOrientedUniformInverse(NodeCoords, Coords);

  for (int iDim = 0; iDim < 3; ++iDim) {
    if (LocalCoords(iDim) < -Tolerance) {
      return false;
    } else if (LocalCoords(iDim) > 1.+Tolerance) {
      return false;
    }
  }

  return true;

}

inline bool OverlapsHexNonUniform(const array_view<const elem<double,3>> &NodeCoords, const
  elem<double,3> &Coords, double Tolerance) {

  // Not safe to invert isoparametric mapping; instead decompose into 6 tetrahedra
  constexpr int Tetrahedra[6][4] = {
    {0,1,2,4},
    {1,2,4,5},
    {2,4,5,6},
    {1,2,5,3},
    {2,5,3,6},
    {3,6,5,7}
  };

  for (int iTetrahedron = 0; iTetrahedron < 6; ++iTetrahedron) {
    int iVertex1 = Tetrahedra[iTetrahedron][0];
    int iVertex2 = Tetrahedra[iTetrahedron][1];
    int iVertex3 = Tetrahedra[iTetrahedron][2];
    int iVertex4 = Tetrahedra[iTetrahedron][3];
    elem<double,3> I = {
      NodeCoords(iVertex2)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex2)(1) - NodeCoords(iVertex1)(1),
      NodeCoords(iVertex2)(2) - NodeCoords(iVertex1)(2)
    };
    elem<double,3> J = {
      NodeCoords(iVertex3)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex3)(1) - NodeCoords(iVertex1)(1),
      NodeCoords(iVertex3)(2) - NodeCoords(iVertex1)(2)
    };
    elem<double,3> K = {
      NodeCoords(iVertex4)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex4)(1) - NodeCoords(iVertex1)(1),
      NodeCoords(iVertex4)(2) - NodeCoords(iVertex1)(2)
    };
    elem<double,3> RelativeCoords = {
      Coords(0) - NodeCoords(iVertex1)(0),
      Coords(1) - NodeCoords(iVertex1)(1),
      Coords(2) - NodeCoords(iVertex1)(2)
    };
    elem<double,3> LocalCoords = ColumnSolve3D(I, J, K, RelativeCoords);
    bool Inside = true;
    double CoordSum = 0.;
    for (int iDim = 0; iDim < 3; ++iDim) {
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

inline double VolumeHexUniform(const elem<double,3> &LowerNodeCoords, const elem<double,3>
  &UpperNodeCoords) {

  return
    (UpperNodeCoords(0) - LowerNodeCoords(0))*
    (UpperNodeCoords(1) - LowerNodeCoords(1))*
    (UpperNodeCoords(2) - LowerNodeCoords(2));

}

inline double VolumeHexOrientedUniform(const array_view<const elem<double,3>> &NodeCoords) {

  elem<double,3> I = {
    NodeCoords(1)(0) - NodeCoords(0)(0),
    NodeCoords(1)(1) - NodeCoords(0)(1),
    NodeCoords(1)(2) - NodeCoords(0)(2)
  };

  elem<double,3> J = {
    NodeCoords(2)(0) - NodeCoords(0)(0),
    NodeCoords(2)(1) - NodeCoords(0)(1),
    NodeCoords(2)(2) - NodeCoords(0)(2)
  };

  elem<double,3> K = {
    NodeCoords(4)(0) - NodeCoords(0)(0),
    NodeCoords(4)(1) - NodeCoords(0)(1),
    NodeCoords(4)(2) - NodeCoords(0)(2)
  };

  return std::abs(ColumnDeterminant3D(I, J, K));

}

inline double VolumeHexNonUniform(const array_view<const elem<double,3>> &NodeCoords) {

  // Decompose into 6 tetrahedra
  constexpr int Tetrahedra[6][4] = {
    {0,1,2,4},
    {1,2,4,5},
    {2,4,5,6},
    {1,2,5,3},
    {2,5,3,6},
    {3,6,5,7}
  };

  double Volume = 0.;

  for (int iTetrahedron = 0; iTetrahedron < 6; ++iTetrahedron) {
    int iVertex1 = Tetrahedra[iTetrahedron][0];
    int iVertex2 = Tetrahedra[iTetrahedron][1];
    int iVertex3 = Tetrahedra[iTetrahedron][2];
    int iVertex4 = Tetrahedra[iTetrahedron][3];
    elem<double,3> I = {
      NodeCoords(iVertex2)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex2)(1) - NodeCoords(iVertex1)(1),
      NodeCoords(iVertex2)(2) - NodeCoords(iVertex1)(2)
    };
    elem<double,3> J = {
      NodeCoords(iVertex3)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex3)(1) - NodeCoords(iVertex1)(1),
      NodeCoords(iVertex3)(2) - NodeCoords(iVertex1)(2)
    };
    elem<double,3> K = {
      NodeCoords(iVertex4)(0) - NodeCoords(iVertex1)(0),
      NodeCoords(iVertex4)(1) - NodeCoords(iVertex1)(1),
      NodeCoords(iVertex4)(2) - NodeCoords(iVertex1)(2)
    };
    Volume += std::abs(ColumnDeterminant3D(I, J, K))/6.;
  }

  return Volume;

}

}}
