// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

namespace overlaps_cell_internal {

inline bool OverlapsCellUniform(int NumDims, const array<field<double>> &Coords, double Tolerance,
  const tuple<int> &Cell, const tuple<double> &PointCoords) {

  tuple<int> UpperCornerOffset = MakeUniformTuple<int>(NumDims, 1, 0);

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

  bool Overlaps;

  switch (NumDims) {
  case 1:
    Overlaps = core::OverlapsLine(LowerCornerCoords(0), UpperCornerCoords(0), PointCoords(0),
      Tolerance);
    break;
  case 2:
    Overlaps = core::OverlapsQuadUniform(LowerCornerCoords, UpperCornerCoords, PointCoords,
      Tolerance);
    break;
  case 3:
    Overlaps = core::OverlapsHexUniform(LowerCornerCoords, UpperCornerCoords, PointCoords,
      Tolerance);
    break;
  }

  return Overlaps;

}

inline bool OverlapsCellOrientedUniform(int NumDims, const array<field<double>> &Coords, double
  Tolerance, const tuple<int> &Cell, const tuple<double> &PointCoords) {

  tuple<double> NodeCoordData[1 << MAX_DIMS];
//   static_array<tuple<double>,1 << MAX_DIMS> NodeCoords;

  range CellExtents;
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    CellExtents.Begin(iDim) = Cell(iDim);
    CellExtents.End(iDim) = Cell(iDim)+2;
  }
  for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
    CellExtents.Begin(iDim) = 0;
    CellExtents.End(iDim) = 1;
  }

  int NumNodes = 0;
  for (int k = CellExtents.Begin(2); k < CellExtents.End(2); ++k) {
    for (int j = CellExtents.Begin(1); j < CellExtents.End(1); ++j) {
      for (int i = CellExtents.Begin(0); i < CellExtents.End(0); ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
        NodeCoordData[NumNodes] = {
          Coords(0)[iPoint],
          Coords(1)[iPoint],
          Coords(2)[iPoint]
        };
        ++NumNodes;
      }
    }
  }

  array_view<const tuple<double>> NodeCoords(NodeCoordData, {NumNodes});

  bool Overlaps;

  switch (NumDims) {
  case 2:
    Overlaps = core::OverlapsQuadOrientedUniform(NodeCoords, PointCoords, Tolerance);
    break;
  case 3:
    Overlaps = core::OverlapsHexOrientedUniform(NodeCoords, PointCoords, Tolerance);
    break;
  }

  return Overlaps;

}

inline bool OverlapsCellNonUniform(int NumDims, const array<field<double>> &Coords, double
  Tolerance, const tuple<int> &Cell, const tuple<double> &PointCoords) {

  tuple<double> NodeCoordData[1 << MAX_DIMS];
//   static_array<tuple<double>,1 << MAX_DIMS> NodeCoords;

  range CellExtents;
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    CellExtents.Begin(iDim) = Cell(iDim);
    CellExtents.End(iDim) = Cell(iDim)+2;
  }
  for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
    CellExtents.Begin(iDim) = 0;
    CellExtents.End(iDim) = 1;
  }

  int NumNodes = 0;
  for (int k = CellExtents.Begin(2); k < CellExtents.End(2); ++k) {
    for (int j = CellExtents.Begin(1); j < CellExtents.End(1); ++j) {
      for (int i = CellExtents.Begin(0); i < CellExtents.End(0); ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
        NodeCoordData[NumNodes] = {
          Coords(0)[iPoint],
          Coords(1)[iPoint],
          Coords(2)[iPoint]
        };
        ++NumNodes;
      }
    }
  }

  array_view<const tuple<double>> NodeCoords(NodeCoordData, {NumNodes});

  bool Overlaps;

  switch (NumDims) {
  case 2:
    Overlaps = core::OverlapsQuadNonUniform(NodeCoords, PointCoords, Tolerance);
    break;
  case 3:
    Overlaps = core::OverlapsHexNonUniform(NodeCoords, PointCoords, Tolerance);
    break;
  }

  return Overlaps;

}

}

inline bool OverlapsCell(int NumDims, const array<field<double>> &Coords, geometry_type
  GeometryType, double Tolerance, const tuple<int> &Cell, const tuple<double> &PointCoords) {

  using overlaps_cell_internal::OverlapsCellUniform;
  using overlaps_cell_internal::OverlapsCellOrientedUniform;
  using overlaps_cell_internal::OverlapsCellNonUniform;

  bool Overlaps;

  switch (GeometryType) {
  case geometry_type::UNIFORM:
  case geometry_type::RECTILINEAR:
    Overlaps = OverlapsCellUniform(NumDims, Coords, Tolerance, Cell, PointCoords);
    break;
  case geometry_type::ORIENTED_UNIFORM:
  case geometry_type::ORIENTED_RECTILINEAR:
    Overlaps = OverlapsCellOrientedUniform(NumDims, Coords, Tolerance, Cell, PointCoords);
    break;
  case geometry_type::CURVILINEAR:
    Overlaps = OverlapsCellNonUniform(NumDims, Coords, Tolerance, Cell, PointCoords);
    break;
  }

  return Overlaps;

}

namespace coords_in_cell_internal {

inline tuple<double> CoordsInCellUniform(int NumDims, const array<field<double>> &Coords, const
  tuple<int> &Cell, const tuple<double> &PointCoords) {

  tuple<int> UpperCornerOffset = MakeUniformTuple<int>(NumDims, 1, 0);

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

  tuple<double> LocalCoords;

  switch (NumDims) {
  case 1:
    LocalCoords(0) = core::IsoLine2NodeInverse(LowerCornerCoords(0), UpperCornerCoords(0),
      PointCoords(0));
    break;
  case 2:
    LocalCoords = core::IsoQuad4NodeUniformInverse(LowerCornerCoords, UpperCornerCoords,
      PointCoords);
    break;
  case 3:
    LocalCoords = core::IsoHex8NodeUniformInverse(LowerCornerCoords, UpperCornerCoords,
      PointCoords);
    break;
  }

  return LocalCoords;

}

inline tuple<double> CoordsInCellOrientedUniform(int NumDims, const array<field<double>> &Coords,
  const tuple<int> &Cell, const tuple<double> &PointCoords) {

  tuple<double> NodeCoordData[1 << MAX_DIMS];
//   static_array<tuple<double>,1 << MAX_DIMS> NodeCoords;

  range CellExtents;
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    CellExtents.Begin(iDim) = Cell(iDim);
    CellExtents.End(iDim) = Cell(iDim)+2;
  }
  for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
    CellExtents.Begin(iDim) = 0;
    CellExtents.End(iDim) = 1;
  }

  int NumNodes = 0;
  for (int k = CellExtents.Begin(2); k < CellExtents.End(2); ++k) {
    for (int j = CellExtents.Begin(1); j < CellExtents.End(1); ++j) {
      for (int i = CellExtents.Begin(0); i < CellExtents.End(0); ++i) {
        long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
        NodeCoordData[NumNodes] = {
          Coords(0)[iPoint],
          Coords(1)[iPoint],
          Coords(2)[iPoint]
        };
        ++NumNodes;
      }
    }
  }

  array_view<const tuple<double>> NodeCoords(NodeCoordData, {NumNodes});

  tuple<double> LocalCoords;

  switch (NumDims) {
  case 2:
    LocalCoords = core::IsoQuad4NodeOrientedUniformInverse(NodeCoords, PointCoords);
    break;
  case 3:
    LocalCoords = core::IsoHex8NodeOrientedUniformInverse(NodeCoords, PointCoords);
    break;
  }

  return LocalCoords;

}

inline tuple<double> CoordsInCellNonUniform(int NumDims, const array<field<double>> &Coords, const
  tuple<int> &Cell, const tuple<double> &PointCoords) {

  tuple<double> LocalCoords;

  if (NumDims == 1) {

    double NodeCoordData[4];

    int ShiftedCell = Cell(0);
    ShiftedCell = Max<int>(ShiftedCell, Coords(0).Extents().Begin(0)+1);
    ShiftedCell = Min<int>(ShiftedCell, Coords(0).Extents().End(0)-3);

    int iNode = 0;
    for (int iPoint = ShiftedCell-1; iPoint < ShiftedCell+3; ++iPoint) {
      NodeCoordData[iNode] = Coords(0)[iPoint];
      ++iNode;
    }

    array_view<const double> NodeCoords(NodeCoordData, {4});

    LocalCoords(0) = core::IsoLine4NodeInverse(NodeCoords, PointCoords(0));
    LocalCoords(1) = 0;
    LocalCoords(2) = 0;

  } else {

    tuple<double> NodeCoordData[1 << 2*MAX_DIMS];
  //   static_array<tuple<double>,1 << 2*MAX_DIMS> NodeCoords;

    tuple<int> ShiftedCell = Cell;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      ShiftedCell(iDim) = Max<int>(ShiftedCell(iDim), Coords(0).Extents().Begin(iDim)+1);
      ShiftedCell(iDim) = Min<int>(ShiftedCell(iDim), Coords(0).Extents().End(iDim)-2);
    }

    range CellExtents;
    for (int iDim = 0; iDim < NumDims; ++iDim) {
      CellExtents.Begin(iDim) = ShiftedCell(iDim)-1;
      CellExtents.End(iDim) = ShiftedCell(iDim)+2;
    }
    for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
      CellExtents.Begin(iDim) = 0;
      CellExtents.End(iDim) = 1;
    }

    int NumNodes = 0;
    for (int k = CellExtents.Begin(2); k < CellExtents.End(2); ++k) {
      for (int j = CellExtents.Begin(1); j < CellExtents.End(1); ++j) {
        for (int i = CellExtents.Begin(0); i < CellExtents.End(0); ++i) {
          long long iPoint = Coords(0).Indexer().ToIndex(i,j,k);
          NodeCoordData[NumNodes] = {
            Coords(0)[iPoint],
            Coords(1)[iPoint],
            Coords(2)[iPoint]
          };
          ++NumNodes;
        }
      }
    }

    array_view<const tuple<double>> NodeCoords(NodeCoordData, {NumNodes});

    switch (NumDims) {
    case 2:
      LocalCoords = core::IsoQuad16NodeInverse(NodeCoords, PointCoords);
      break;
    case 3:
      LocalCoords = core::IsoHex64NodeInverse(NodeCoords, PointCoords);
      break;
    }

  }

  return LocalCoords;

}

}

inline tuple<double> CoordsInCell(int NumDims, const array<field<double>> &Coords, geometry_type
  GeometryType, const tuple<int> &Cell, const tuple<double> &PointCoords) {

  using coords_in_cell_internal::CoordsInCellUniform;
  using coords_in_cell_internal::CoordsInCellOrientedUniform;
  using coords_in_cell_internal::CoordsInCellNonUniform;

  tuple<double> LocalCoords;

  switch (GeometryType) {
  case geometry_type::UNIFORM:
    LocalCoords = CoordsInCellUniform(NumDims, Coords, Cell, PointCoords);
    break;
  case geometry_type::ORIENTED_UNIFORM:
    LocalCoords = CoordsInCellOrientedUniform(NumDims, Coords, Cell, PointCoords);
    break;
  case geometry_type::RECTILINEAR:
  case geometry_type::ORIENTED_RECTILINEAR:
  case geometry_type::CURVILINEAR:
    LocalCoords = CoordsInCellNonUniform(NumDims, Coords, Cell, PointCoords);
    break;
  }

  return LocalCoords;

}

}}
