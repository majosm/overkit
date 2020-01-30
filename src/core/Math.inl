// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

namespace ovk {
namespace core {

inline bool IsNaN(double Value) {

// std::isnan is unreliable if using GCC's -ffast-math, so just disable NaN handling
#ifndef __FAST_MATH__
  return std::isnan(Value);
#else
  return false;
#endif

}

template <int N> bool IsNaN(const elem<double,N> &Value) {

  bool HasNaN = false;

  for (int iDim = 0; iDim < N; ++iDim) {
    HasNaN = HasNaN || IsNaN(Value(iDim));
  }

  return HasNaN;

}

inline double ColumnDeterminant2D(const elem<double,2> &AI, const elem<double,2> &AJ) {

  return AI(0)*AJ(1) - AI(1)*AJ(0);

}

inline double ColumnDeterminant3D(const elem<double,3> &AI, const elem<double,3> &AJ, const
  elem<double,3> &AK) {

  return
    AI(0) * (AJ(1)*AK(2) - AJ(2)*AK(1)) +
    AJ(0) * (AK(1)*AI(2) - AK(2)*AI(1)) +
    AK(0) * (AI(1)*AJ(2) - AI(2)*AJ(1));

}

inline elem<double,2> ColumnSolve2D(const elem<double,2> &AI, const elem<double,2> &AJ, const
  elem<double,2> &B) {

  elem<double,2> X;

  double Det = ColumnDeterminant2D(AI, AJ);

  X(0) = ColumnDeterminant2D(B, AJ)/Det;
  X(1) = ColumnDeterminant2D(AI, B)/Det;

  return X;

}

inline elem<double,3> ColumnSolve3D(const elem<double,3> &AI, const elem<double,3> &AJ, const
  elem<double,3> &AK, const elem<double,3> &B) {

  elem<double,3> X;

  double Det = ColumnDeterminant3D(AI, AJ, AK);

  X(0) = ColumnDeterminant3D(B, AJ, AK)/Det;
  X(1) = ColumnDeterminant3D(AI, B, AK)/Det;
  X(2) = ColumnDeterminant3D(AI, AJ, B)/Det;

  return X;

}

inline elem<double,2> LagrangeInterpLinear(double U) {

  return {1. - U, U};

}

inline elem<double,2> LagrangeInterpLinearDeriv(double U) {

  return {-1., 1.};

}

inline elem<double,4> LagrangeInterpCubic(double U) {

  return {
    -(       U * (U - 1.) * (U - 2.)) / 6.,
     ((U + 1.) * (U - 1.) * (U - 2.)) / 2.,
    -((U + 1.) *        U * (U - 2.)) / 2.,
     ((U + 1.) *        U * (U - 1.)) / 6.
  };

}

inline elem<double,4> LagrangeInterpCubicDeriv(double U) {

  constexpr double Sqrt3 = 1.732050807568877293527446341505872366942805253810380628055;
  constexpr double Sqrt7 = 2.645751311064590590501615753639260425710259183082450180368;

  constexpr double Roots[] = {
      1. + 1./Sqrt3,
      1. - 1./Sqrt3,
    (2. + Sqrt7)/3.,
    (2. - Sqrt7)/3.,
    (1. + Sqrt7)/3.,
    (1. - Sqrt7)/3.,
           1./Sqrt3,
          -1./Sqrt3
  };

  return {
         -((U - Roots[0]) * (U - Roots[1]))/2.,
     (3. * (U - Roots[2]) * (U - Roots[3]))/2.,
    -(3. * (U - Roots[4]) * (U - Roots[5]))/2.,
          ((U - Roots[6]) * (U - Roots[7]))/2.
  };

}

inline tuple<int> MapToUniformGridCell(int NumDims, const tuple<int> &Origin, const tuple<int>
  &CellSize, const tuple<int> &Point) {

  tuple<int> Cell = MakeUniformTuple<int>(NumDims, 0);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    int Offset = Point(iDim) - Origin(iDim);
    // Division rounding down
    Cell(iDim) = Offset/CellSize(iDim) - (Offset % CellSize(iDim) < 0);
  }

  return Cell;

}

inline tuple<int> MapToUniformGridCell(int NumDims, const tuple<double> &Origin, const tuple<double>
  &CellSize, const tuple<double> &Point) {

  tuple<int> Cell = MakeUniformTuple<int>(NumDims, 0);

  for (int iDim = 0; iDim < NumDims; ++iDim) {
    double Offset = Point(iDim) - Origin(iDim);
    Cell(iDim) = int(std::floor(Offset/CellSize(iDim)));
  }

  return Cell;

}

}}
