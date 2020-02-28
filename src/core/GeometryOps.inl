// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

// Have to declare these here because they are called in other functions before they are defined
// (otherwise get "specialization of '...' after instantiation" errors)
template <> inline optional<tuple<double>> CoordsInCell<geometry_type::CURVILINEAR, 1>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords);
template <> inline optional<tuple<double>> CoordsInCell<geometry_type::CURVILINEAR, 2>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords);
template <> inline optional<tuple<double>> CoordsInCell<geometry_type::CURVILINEAR, 3>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords);
template <> inline box CellBounds<geometry_type::CURVILINEAR, 2>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell);
template <> inline box CellBounds<geometry_type::CURVILINEAR, 3>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell);

template <> inline bool OverlapsCell<geometry_type::UNIFORM, 1>(const array_view<const field_view<
  const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  double LowerCornerCoord = Coords(0)(Cell);
  double UpperCornerCoord = Coords(0)(Cell(0)+1,0,0);

  double PointCoord = PointCoords(0);

  return OverlapsLine(LowerCornerCoord, UpperCornerCoord, PointCoord, Tolerance);

}

template <> inline bool OverlapsCell<geometry_type::UNIFORM, 2>(const array_view<const field_view<
  const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  constexpr tuple<int> UpperCornerOffset = {1,1,0};

  long long iLowerCorner = Coords(0).Indexer().ToIndex(Cell);
  long long iUpperCorner = Coords(0).Indexer().ToIndex(Cell+UpperCornerOffset);

  elem<double,2> LowerCornerCoords = {
    Coords(0)[iLowerCorner],
    Coords(1)[iLowerCorner]
  };
  elem<double,2> UpperCornerCoords = {
    Coords(0)[iUpperCorner],
    Coords(1)[iUpperCorner]
  };

  elem<double,2> PointCoords_ = {PointCoords(0), PointCoords(1)};

  return OverlapsQuadUniform(LowerCornerCoords, UpperCornerCoords, PointCoords_, Tolerance);

}

template <> inline bool OverlapsCell<geometry_type::UNIFORM, 3>(const array_view<const field_view<
  const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  constexpr tuple<int> UpperCornerOffset = {1,1,1};

  long long iLowerCorner = Coords(0).Indexer().ToIndex(Cell);
  long long iUpperCorner = Coords(0).Indexer().ToIndex(Cell+UpperCornerOffset);

  tuple<double> LowerCornerCoords = {
    Coords(0)[iLowerCorner],
    Coords(1)[iLowerCorner],
    Coords(2)[iLowerCorner]
  };
  tuple<double> UpperCornerCoords = {
    Coords(0)[iUpperCorner],
    Coords(1)[iUpperCorner],
    Coords(2)[iUpperCorner]
  };

  return OverlapsHexUniform(LowerCornerCoords, UpperCornerCoords, PointCoords, Tolerance);

}

template <> inline bool OverlapsCell<geometry_type::RECTILINEAR, 1>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return OverlapsCell<geometry_type::UNIFORM, 1>(Coords, Tolerance, Cell, PointCoords);

}

template <> inline bool OverlapsCell<geometry_type::RECTILINEAR, 2>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return OverlapsCell<geometry_type::UNIFORM, 2>(Coords, Tolerance, Cell, PointCoords);

}

template <> inline bool OverlapsCell<geometry_type::RECTILINEAR, 3>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return OverlapsCell<geometry_type::UNIFORM, 3>(Coords, Tolerance, Cell, PointCoords);

}

template <> inline bool OverlapsCell<geometry_type::ORIENTED_UNIFORM, 1>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return OverlapsCell<geometry_type::UNIFORM, 1>(Coords, Tolerance, Cell, PointCoords);

}

template <> inline bool OverlapsCell<geometry_type::ORIENTED_UNIFORM, 2>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  elem<double,2> NodeCoords[4];
  int iNode = 0;
  for (int j = Cell(1); j <= Cell(1)+1; ++j) {
    for (int i = Cell(0); i <= Cell(0)+1; ++i) {
      long long iPoint = Coords(0).Indexer().ToIndex(i,j,0);
      NodeCoords[iNode] = {
        Coords(0)[iPoint],
        Coords(1)[iPoint]
      };
      ++iNode;
    }
  }

  elem<double,2> PointCoords_ = {PointCoords(0), PointCoords(1)};

  return OverlapsQuadOrientedUniform(NodeCoords, PointCoords_, Tolerance);

}

template <> inline bool OverlapsCell<geometry_type::ORIENTED_UNIFORM, 3>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  tuple<double> NodeCoords[8];
  int iNode = 0;
  for (int k = Cell(2); k <= Cell(2)+1; ++k) {
    for (int j = Cell(1); j <= Cell(1)+1; ++j) {
      for (int i = Cell(0); i <= Cell(0)+1; ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
        NodeCoords[iNode] = {
          Coords(0)[iPoint],
          Coords(1)[iPoint],
          Coords(2)[iPoint]
        };
        ++iNode;
      }
    }
  }

  return OverlapsHexOrientedUniform(NodeCoords, PointCoords, Tolerance);

}

template <> inline bool OverlapsCell<geometry_type::ORIENTED_RECTILINEAR, 1>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return OverlapsCell<geometry_type::UNIFORM, 1>(Coords, Tolerance, Cell, PointCoords);

}

template <> inline bool OverlapsCell<geometry_type::ORIENTED_RECTILINEAR, 2>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return OverlapsCell<geometry_type::ORIENTED_UNIFORM, 2>(Coords, Tolerance, Cell, PointCoords);

}

template <> inline bool OverlapsCell<geometry_type::ORIENTED_RECTILINEAR, 3>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return OverlapsCell<geometry_type::ORIENTED_UNIFORM, 3>(Coords, Tolerance, Cell, PointCoords);

}

template <> inline bool OverlapsCell<geometry_type::CURVILINEAR, 1>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return OverlapsCell<geometry_type::UNIFORM, 1>(Coords, Tolerance, Cell, PointCoords);

}

template <> inline bool OverlapsCell<geometry_type::CURVILINEAR, 2>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  elem<double,2> NodeCoords[4];
  int iNode = 0;
  for (int j = Cell(1); j <= Cell(1)+1; ++j) {
    for (int i = Cell(0); i <= Cell(0)+1; ++i) {
      long long iPoint = Coords(0).Indexer().ToIndex(i,j,0);
      NodeCoords[iNode] = {
        Coords(0)[iPoint],
        Coords(1)[iPoint]
      };
      ++iNode;
    }
  }

  elem<double,2> PointCoords_ = {PointCoords(0), PointCoords(1)};

  return OverlapsQuadNonUniform(NodeCoords, PointCoords_, Tolerance);

}

template <> inline bool OverlapsCell<geometry_type::CURVILINEAR, 3>(const array_view<const
  field_view<const double>> &Coords, double Tolerance, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  tuple<double> NodeCoords[8];
  int iNode = 0;
  for (int k = Cell(2); k <= Cell(2)+1; ++k) {
    for (int j = Cell(1); j <= Cell(1)+1; ++j) {
      for (int i = Cell(0); i <= Cell(0)+1; ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
        NodeCoords[iNode] = {
          Coords(0)[iPoint],
          Coords(1)[iPoint],
          Coords(2)[iPoint]
        };
        ++iNode;
      }
    }
  }

  return OverlapsHexNonUniform(NodeCoords, PointCoords, Tolerance);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::UNIFORM, 1>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  double LowerCornerCoord = Coords(0)(Cell);
  double UpperCornerCoord = Coords(0)(Cell(0)+1,0,0);

  double PointCoord = PointCoords(0);

  double LocalCoord = IsoLine2NodeInverse(LowerCornerCoord, UpperCornerCoord, PointCoord);

  optional<tuple<double>> MaybeLocalCoords;
  tuple<double> &LocalCoords = *(MaybeLocalCoords.Assign());
  LocalCoords(0) = LocalCoord;
  LocalCoords(1) = 0.;
  LocalCoords(2) = 0.;

  return MaybeLocalCoords;

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::UNIFORM, 2>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  constexpr tuple<int> UpperCornerOffset = {1,1,0};

  long long iLowerCorner = Coords(0).Indexer().ToIndex(Cell);
  long long iUpperCorner = Coords(0).Indexer().ToIndex(Cell+UpperCornerOffset);

  elem<double,2> LowerCornerCoords = {
    Coords(0)[iLowerCorner],
    Coords(1)[iLowerCorner]
  };
  elem<double,2> UpperCornerCoords = {
    Coords(0)[iUpperCorner],
    Coords(1)[iUpperCorner]
  };

  elem<double,2> PointCoords_ = {PointCoords(0), PointCoords(1)};

  elem<double,2> LocalCoords_ = IsoQuad4NodeUniformInverse(LowerCornerCoords, UpperCornerCoords,
    PointCoords_);

  optional<tuple<double>> MaybeLocalCoords;
  tuple<double> &LocalCoords = *(MaybeLocalCoords.Assign());
  LocalCoords(0) = LocalCoords_(0);
  LocalCoords(1) = LocalCoords_(1);
  LocalCoords(2) = 0.;

  return MaybeLocalCoords;

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::UNIFORM, 3>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  constexpr tuple<int> UpperCornerOffset = {1,1,1};

  long long iLowerCorner = Coords(0).Indexer().ToIndex(Cell);
  long long iUpperCorner = Coords(0).Indexer().ToIndex(Cell+UpperCornerOffset);

  tuple<double> LowerCornerCoords = {
    Coords(0)[iLowerCorner],
    Coords(1)[iLowerCorner],
    Coords(2)[iLowerCorner]
  };
  tuple<double> UpperCornerCoords = {
    Coords(0)[iUpperCorner],
    Coords(1)[iUpperCorner],
    Coords(2)[iUpperCorner]
  };

  return IsoHex8NodeUniformInverse(LowerCornerCoords, UpperCornerCoords, PointCoords);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::RECTILINEAR, 1>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return CoordsInCell<geometry_type::CURVILINEAR, 1>(Coords, Cell, PointCoords);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::RECTILINEAR, 2>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return CoordsInCell<geometry_type::CURVILINEAR, 2>(Coords, Cell, PointCoords);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::RECTILINEAR, 3>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return CoordsInCell<geometry_type::CURVILINEAR, 3>(Coords, Cell, PointCoords);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::ORIENTED_UNIFORM, 1>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  return CoordsInCell<geometry_type::UNIFORM, 1>(Coords, Cell, PointCoords);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::ORIENTED_UNIFORM, 2>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  elem<double,2> NodeCoords[4];
  int iNode = 0;
  for (int j = Cell(1); j <= Cell(1)+1; ++j) {
    for (int i = Cell(0); i <= Cell(0)+1; ++i) {
      long long iPoint = Coords(0).Indexer().ToIndex(i,j,0);
      NodeCoords[iNode] = {
        Coords(0)[iPoint],
        Coords(1)[iPoint]
      };
      ++iNode;
    }
  }

  elem<double,2> PointCoords_ = {PointCoords(0), PointCoords(1)};

  elem<double,2> LocalCoords_ = IsoQuad4NodeOrientedUniformInverse(NodeCoords, PointCoords_);

  optional<tuple<double>> MaybeLocalCoords;
  tuple<double> &LocalCoords = *(MaybeLocalCoords.Assign());
  LocalCoords(0) = LocalCoords_(0);
  LocalCoords(1) = LocalCoords_(1);
  LocalCoords(2) = 0.;

  return MaybeLocalCoords;

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::ORIENTED_UNIFORM, 3>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  tuple<double> NodeCoords[8];
  int iNode = 0;
  for (int k = Cell(2); k <= Cell(2)+1; ++k) {
    for (int j = Cell(1); j <= Cell(1)+1; ++j) {
      for (int i = Cell(0); i <= Cell(0)+1; ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
        NodeCoords[iNode] = {
          Coords(0)[iPoint],
          Coords(1)[iPoint],
          Coords(2)[iPoint]
        };
        ++iNode;
      }
    }
  }

  return IsoHex8NodeOrientedUniformInverse(NodeCoords, PointCoords);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::ORIENTED_RECTILINEAR, 1>(
  const array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const
  tuple<double> &PointCoords) {

  return CoordsInCell<geometry_type::CURVILINEAR, 1>(Coords, Cell, PointCoords);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::ORIENTED_RECTILINEAR, 2>(
  const array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const
  tuple<double> &PointCoords) {

  return CoordsInCell<geometry_type::CURVILINEAR, 2>(Coords, Cell, PointCoords);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::ORIENTED_RECTILINEAR, 3>(
  const array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const
  tuple<double> &PointCoords) {

  return CoordsInCell<geometry_type::CURVILINEAR, 3>(Coords, Cell, PointCoords);

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::CURVILINEAR, 1>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  optional<tuple<double>> MaybeLocalCoords;

  int ShiftedCell = Cell(0);
  ShiftedCell = Max<int>(ShiftedCell, Coords(0).Extents().Begin(0)+1);
  ShiftedCell = Min<int>(ShiftedCell, Coords(0).Extents().End(0)-3);

  bool StencilFits = ShiftedCell-1 < Coords(0).Extents().Begin(0) &&
    ShiftedCell+2 >= Coords(0).Extents().End(0);

  if (StencilFits) {

    double NodeCoords[4];
    int iNode = 0;
    for (int iPoint = ShiftedCell-1; iPoint <= ShiftedCell+2; ++iPoint) {
      NodeCoords[iNode] = Coords(0)[iPoint];
      ++iNode;
    }

    double PointCoord = PointCoords(0);

    auto MaybeLocalCoord = IsoLine4NodeInverse(NodeCoords, PointCoord);

    if (MaybeLocalCoord) {
      double LocalCoord = *MaybeLocalCoord;
      LocalCoord += double(ShiftedCell - Cell(0));
      tuple<double> &LocalCoords = *(MaybeLocalCoords.Assign());
      LocalCoords(0) = LocalCoord;
      LocalCoords(1) = 0.;
      LocalCoords(2) = 0.;
    }

  }

  return MaybeLocalCoords;

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::CURVILINEAR, 2>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  optional<tuple<double>> MaybeLocalCoords;

  elem<int,2> ShiftedCell = {Cell(0),Cell(1)};
  bool StencilFits = true;
  for (int iDim = 0; iDim < 2; ++iDim) {
    ShiftedCell(iDim) = Max<int>(ShiftedCell(iDim), Coords(0).Extents().Begin(iDim)+1);
    ShiftedCell(iDim) = Min<int>(ShiftedCell(iDim), Coords(0).Extents().End(iDim)-3);
    if (ShiftedCell(iDim)-1 < Coords(0).Extents().Begin(iDim) || ShiftedCell(iDim)+2 >=
      Coords(0).Extents().End(iDim)) {
      StencilFits = false;
      break;
    }
  }

  if (StencilFits) {

    elem<double,2> NodeCoords[16];
    int iNode = 0;
    for (int j = ShiftedCell(1)-1; j <= ShiftedCell(1)+2; ++j) {
      for (int i = ShiftedCell(0)-1; i <= ShiftedCell(0)+2; ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,0);
        NodeCoords[iNode] = {
          Coords(0)[iPoint],
          Coords(1)[iPoint]
        };
        ++iNode;
      }
    }

    elem<double,2> PointCoords_ = {PointCoords(0), PointCoords(1)};

    auto MaybeLocalCoords_ = IsoQuad16NodeInverse(NodeCoords, PointCoords_);

    if (MaybeLocalCoords_) {
      elem<double,2> &LocalCoords_ = *MaybeLocalCoords_;
      LocalCoords_(0) += ShiftedCell(0) - Cell(0);
      LocalCoords_(1) += ShiftedCell(1) - Cell(1);
      tuple<double> &LocalCoords  = *(MaybeLocalCoords.Assign());
      LocalCoords(0) = LocalCoords_(0);
      LocalCoords(1) = LocalCoords_(1);
      LocalCoords(2) = 0.;
    }

  }

  return MaybeLocalCoords;

}

template <> inline optional<tuple<double>> CoordsInCell<geometry_type::CURVILINEAR, 3>(const
  array_view<const field_view<const double>> &Coords, const tuple<int> &Cell, const tuple<double>
  &PointCoords) {

  optional<tuple<double>> MaybeLocalCoords;

  tuple<int> ShiftedCell = Cell;
  bool StencilFits = true;
  for (int iDim = 0; iDim < 3; ++iDim) {
    ShiftedCell(iDim) = Max<int>(ShiftedCell(iDim), Coords(0).Extents().Begin(iDim)+1);
    ShiftedCell(iDim) = Min<int>(ShiftedCell(iDim), Coords(0).Extents().End(iDim)-3);
    if (ShiftedCell(iDim)-1 < Coords(0).Extents().Begin(iDim) || ShiftedCell(iDim)+2 >=
      Coords(0).Extents().End(iDim)) {
      StencilFits = false;
      break;
    }
  }

  if (StencilFits) {

    tuple<double> NodeCoords[64];
    int iNode = 0;
    for (int k = ShiftedCell(2)-1; k <= ShiftedCell(2)+2; ++k) {
      for (int j = ShiftedCell(1)-1; j <= ShiftedCell(1)+2; ++j) {
        for (int i = ShiftedCell(0)-1; i <= ShiftedCell(0)+2; ++i) {
          long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
          NodeCoords[iNode] = {
            Coords(0)[iPoint],
            Coords(1)[iPoint],
            Coords(2)[iPoint]
          };
          ++iNode;
        }
      }
    }

    MaybeLocalCoords = IsoHex64NodeInverse(NodeCoords, PointCoords);

    if (MaybeLocalCoords) {
      tuple<double> &LocalCoords = *MaybeLocalCoords;
      LocalCoords += ShiftedCell - Cell;
    }

  }

  return MaybeLocalCoords;

}

template <> inline double CellVolume<geometry_type::UNIFORM, 1>(const array_view<const field_view<
  const double>> &Coords, const tuple<int> &Cell) {

  double LowerCornerCoord = Coords(0)(Cell(0),0,0);
  double UpperCornerCoord = Coords(0)(Cell(0)+1,0,0);

  return VolumeLine(LowerCornerCoord, UpperCornerCoord);

}

template <> inline double CellVolume<geometry_type::UNIFORM, 2>(const array_view<const field_view<
  const double>> &Coords, const tuple<int> &Cell) {

  constexpr tuple<int> UpperCornerOffset = {1,1,0};

  long long iLowerCorner = Coords(0).Indexer().ToIndex(Cell);
  long long iUpperCorner = Coords(0).Indexer().ToIndex(Cell+UpperCornerOffset);

  elem<double,2> LowerCornerCoords = {
    Coords(0)[iLowerCorner],
    Coords(1)[iLowerCorner]
  };
  elem<double,2> UpperCornerCoords = {
    Coords(0)[iUpperCorner],
    Coords(1)[iUpperCorner]
  };

  return VolumeQuadUniform(LowerCornerCoords, UpperCornerCoords);

}

template <> inline double CellVolume<geometry_type::UNIFORM, 3>(const array_view<const field_view<
  const double>> &Coords, const tuple<int> &Cell) {

  constexpr tuple<int> UpperCornerOffset = {1,1,1};

  long long iLowerCorner = Coords(0).Indexer().ToIndex(Cell);
  long long iUpperCorner = Coords(0).Indexer().ToIndex(Cell+UpperCornerOffset);

  tuple<double> LowerCornerCoords = {
    Coords(0)[iLowerCorner],
    Coords(1)[iLowerCorner],
    Coords(2)[iLowerCorner]
  };
  tuple<double> UpperCornerCoords = {
    Coords(0)[iUpperCorner],
    Coords(1)[iUpperCorner],
    Coords(2)[iUpperCorner]
  };

  return VolumeHexUniform(LowerCornerCoords, UpperCornerCoords);

}

template <> inline double CellVolume<geometry_type::RECTILINEAR, 1>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellVolume<geometry_type::UNIFORM, 1>(Coords, Cell);

}

template <> inline double CellVolume<geometry_type::RECTILINEAR, 2>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellVolume<geometry_type::UNIFORM, 2>(Coords, Cell);

}

template <> inline double CellVolume<geometry_type::RECTILINEAR, 3>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellVolume<geometry_type::UNIFORM, 3>(Coords, Cell);

}

template <> inline double CellVolume<geometry_type::ORIENTED_UNIFORM, 1>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellVolume<geometry_type::UNIFORM, 1>(Coords, Cell);

}

template <> inline double CellVolume<geometry_type::ORIENTED_UNIFORM, 2>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  elem<double,2> NodeCoords[4];
  int iNode = 0;
  for (int j = Cell(1); j <= Cell(1)+1; ++j) {
    for (int i = Cell(0); i <= Cell(0)+1; ++i) {
      long long iPoint = Coords(0).Indexer().ToIndex(i,j,0);
      NodeCoords[iNode] = {
        Coords(0)[iPoint],
        Coords(1)[iPoint]
      };
      ++iNode;
    }
  }

  return VolumeQuadOrientedUniform(NodeCoords);

}

template <> inline double CellVolume<geometry_type::ORIENTED_UNIFORM, 3>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  tuple<double> NodeCoords[8];
  int iNode = 0;
  for (int k = Cell(2); k <= Cell(2)+1; ++k) {
    for (int j = Cell(1); j <= Cell(1)+1; ++j) {
      for (int i = Cell(0); i <= Cell(0)+1; ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
        NodeCoords[iNode] = {
          Coords(0)[iPoint],
          Coords(1)[iPoint],
          Coords(2)[iPoint]
        };
        ++iNode;
      }
    }
  }

  return VolumeHexOrientedUniform(NodeCoords);

}

template <> inline double CellVolume<geometry_type::ORIENTED_RECTILINEAR, 1>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellVolume<geometry_type::UNIFORM, 1>(Coords, Cell);

}

template <> inline double CellVolume<geometry_type::ORIENTED_RECTILINEAR, 2>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellVolume<geometry_type::ORIENTED_UNIFORM, 2>(Coords, Cell);

}

template <> inline double CellVolume<geometry_type::ORIENTED_RECTILINEAR, 3>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellVolume<geometry_type::ORIENTED_UNIFORM, 3>(Coords, Cell);

}

template <> inline double CellVolume<geometry_type::CURVILINEAR, 1>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellVolume<geometry_type::UNIFORM, 1>(Coords, Cell);

}

template <> inline double CellVolume<geometry_type::CURVILINEAR, 2>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  elem<double,2> NodeCoords[4];
  int iNode = 0;
  for (int j = Cell(1); j <= Cell(1)+1; ++j) {
    for (int i = Cell(0); i <= Cell(0)+1; ++i) {
      long long iPoint = Coords(0).Indexer().ToIndex(i,j,0);
      NodeCoords[iNode] = {
        Coords(0)[iPoint],
        Coords(1)[iPoint]
      };
      ++iNode;
    }
  }

  return VolumeQuadNonUniform(NodeCoords);

}

template <> inline double CellVolume<geometry_type::CURVILINEAR, 3>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  tuple<double> NodeCoords[8];
  int iNode = 0;
  for (int k = Cell(2); k <= Cell(2)+1; ++k) {
    for (int j = Cell(1); j <= Cell(1)+1; ++j) {
      for (int i = Cell(0); i <= Cell(0)+1; ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
        NodeCoords[iNode] = {
          Coords(0)[iPoint],
          Coords(1)[iPoint],
          Coords(2)[iPoint]
        };
        ++iNode;
      }
    }
  }

  return VolumeHexNonUniform(NodeCoords);

}

template <> inline box CellBounds<geometry_type::UNIFORM, 1>(const array_view<const field_view<
  const double>> &Coords, const tuple<int> &Cell) {

  box Bounds;

  Bounds.Begin(0) = Coords(0)(Cell(0),0,0);
  Bounds.End(0) = Coords(0)(Cell(0)+1,0,0);

  return Bounds;

}

template <> inline box CellBounds<geometry_type::UNIFORM, 2>(const array_view<const field_view<
  const double>> &Coords, const tuple<int> &Cell) {

  constexpr tuple<int> UpperCornerOffset = {1,1,0};

  long long iLowerCorner = Coords(0).Indexer().ToIndex(Cell);
  long long iUpperCorner = Coords(0).Indexer().ToIndex(Cell+UpperCornerOffset);

  box Bounds;

  Bounds.Begin() = {
    Coords(0)[iLowerCorner],
    Coords(1)[iLowerCorner],
    0.
  };
  Bounds.End() = {
    Coords(0)[iUpperCorner],
    Coords(1)[iUpperCorner],
    0.
  };

  return Bounds;

}

template <> inline box CellBounds<geometry_type::UNIFORM, 3>(const array_view<const field_view<
  const double>> &Coords, const tuple<int> &Cell) {

  constexpr tuple<int> UpperCornerOffset = {1,1,1};

  long long iLowerCorner = Coords(0).Indexer().ToIndex(Cell);
  long long iUpperCorner = Coords(0).Indexer().ToIndex(Cell+UpperCornerOffset);

  box Bounds;

  Bounds.Begin() = {
    Coords(0)[iLowerCorner],
    Coords(1)[iLowerCorner],
    Coords(2)[iLowerCorner]
  };
  Bounds.End() = {
    Coords(0)[iUpperCorner],
    Coords(1)[iUpperCorner],
    Coords(2)[iUpperCorner]
  };

  return Bounds;

}

template <> inline box CellBounds<geometry_type::RECTILINEAR, 1>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::UNIFORM, 1>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::RECTILINEAR, 2>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::UNIFORM, 2>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::RECTILINEAR, 3>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::UNIFORM, 3>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::ORIENTED_UNIFORM, 1>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::UNIFORM, 1>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::ORIENTED_UNIFORM, 2>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::CURVILINEAR, 2>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::ORIENTED_UNIFORM, 3>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::CURVILINEAR, 3>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::ORIENTED_RECTILINEAR, 1>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::UNIFORM, 1>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::ORIENTED_RECTILINEAR, 2>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::CURVILINEAR, 2>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::ORIENTED_RECTILINEAR, 3>(const array_view<const
  field_view<const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::CURVILINEAR, 3>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::CURVILINEAR, 1>(const array_view<const field_view<
  const double>> &Coords, const tuple<int> &Cell) {

  return CellBounds<geometry_type::UNIFORM, 1>(Coords, Cell);

}

template <> inline box CellBounds<geometry_type::CURVILINEAR, 2>(const array_view<const field_view<
  const double>> &Coords, const tuple<int> &Cell) {

  box Bounds = MakeEmptyBox(2);

  for (int j = Cell(1); j <= Cell(1)+1; ++j) {
    for (int i = Cell(0); i <= Cell(0)+1; ++i) {
      long long iPoint = Coords(0).Indexer().ToIndex(i,j,0);
      tuple<double> NodeCoords = {
        Coords(0)[iPoint],
        Coords(1)[iPoint],
        0.
      };
      Bounds = ExtendBox(Bounds, NodeCoords);
    }
  }

  return Bounds;

}

template <> inline box CellBounds<geometry_type::CURVILINEAR, 3>(const array_view<const field_view<
  const double>> &Coords, const tuple<int> &Cell) {

  box Bounds = MakeEmptyBox(3);

  for (int k = Cell(2); k <= Cell(2)+1; ++k) {
    for (int j = Cell(1); j <= Cell(1)+1; ++j) {
      for (int i = Cell(0); i <= Cell(0)+1; ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
        tuple<double> NodeCoords = {
          Coords(0)[iPoint],
          Coords(1)[iPoint],
          Coords(2)[iPoint]
        };
        Bounds = ExtendBox(Bounds, NodeCoords);
      }
    }
  }

  return Bounds;

}

}}
