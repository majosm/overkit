// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "tests/fixtures/WavyInWavy.hpp"

#include "support/Constants.hpp"
#include "support/Decomp.hpp"

#include <overkit.hpp>

#include <mpi.h>

#include <cmath>

using support::PI;

namespace tests {

ovk::domain WavyInWavy(int NumDims, ovk::comm_view Comm, int Size, bool PreCutHole) {

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(Comm)
    .SetStatusLoggingThreshold(0)
  ));

  ovk::domain Domain = ovk::CreateDomain(std::move(Context), ovk::domain::params()
    .SetDimension(NumDims)
    .SetComm(Comm)
  );

  OVK_DEBUG_ASSERT(Size % 2 == 0, "Size must be even.");
  OVK_DEBUG_ASSERT(!PreCutHole || Size >= 8, "Size is too small; no hole will be cut.");

  ovk::tuple<int> BackgroundSize = {1,1,1};
  ovk::tuple<int> ForegroundSize = {1,1,1};
  for (int iDim = 0; iDim < NumDims; ++iDim) {
    BackgroundSize(iDim) = Size;
    ForegroundSize(iDim) = Size;
  }

  bool BackgroundIsLocal = Comm.Rank() < ovk::Max(Comm.Size()/2, 1);
  bool ForegroundIsLocal = Comm.Rank() >= Comm.Size()/2;

  ovk::comm BackgroundComm = ovk::CreateSubsetComm(Comm, BackgroundIsLocal);
  if (BackgroundIsLocal) {
    ovk::tuple<int> CartDims = ovk::MakeUniformTuple<int>(NumDims, 0, 1);
    BackgroundComm = ovk::CreateCartComm(BackgroundComm, NumDims, CartDims, {false,false,false});
  }

  ovk::comm ForegroundComm = ovk::CreateSubsetComm(Comm, ForegroundIsLocal);
  if (ForegroundIsLocal) {
    ovk::tuple<int> CartDims = ovk::MakeUniformTuple<int>(NumDims, 0, 1);
    ForegroundComm = ovk::CreateCartComm(ForegroundComm, NumDims, CartDims, {false,false,false});
  }

  ovk::array<int> GridIDs({2}, {1, 2});
  ovk::array<ovk::optional<ovk::grid::params>> MaybeGridParams({2});

  if (BackgroundIsLocal) {
    ovk::range LocalRange = support::CartesianDecomp(NumDims, {BackgroundSize}, BackgroundComm);
    MaybeGridParams(0) = Domain.MakeGridParams()
      .SetName("Background")
      .SetComm(BackgroundComm)
      .SetGlobalRange({BackgroundSize})
      .SetLocalRange(LocalRange);
  }

  if (ForegroundIsLocal) {
    ovk::range LocalRange = support::CartesianDecomp(NumDims, {ForegroundSize}, ForegroundComm);
    MaybeGridParams(1) = Domain.MakeGridParams()
      .SetName("Foreground")
      .SetComm(ForegroundComm)
      .SetGlobalRange({ForegroundSize})
      .SetLocalRange(LocalRange);
  }

  Domain.CreateGrids(GridIDs, MaybeGridParams);

  Domain.CreateComponent<ovk::geometry_component>(1);

  {

    auto GeometryComponentEditHandle = Domain.EditComponent<ovk::geometry_component>(1);
    ovk::geometry_component &GeometryComponent = *GeometryComponentEditHandle;

    GeometryComponent.CreateGeometries(GridIDs);

    double WaveAmplitude = 0.125;

    auto CoordFunc = [&](const ovk::tuple<double> &Xi) -> ovk::tuple<double> {
      ovk::tuple<double> Coords = {0.,0.,0.};
      switch (NumDims) {
      case 1: {
        double U = Xi(0)/double(Size-1);
        Coords(0) = (1.-U)*(-1.) + U + WaveAmplitude*std::sin(2.*PI*U);
        break;
      }
      case 2: {
        double U = Xi(0)/double(Size-1);
        double V = Xi(1)/double(Size-1);
        Coords(0) = (1.-U)*(-1.) + U + WaveAmplitude*std::sin(2.*PI*V);
        Coords(1) = (1.-V)*(-1.) + V + WaveAmplitude*std::sin(2.*PI*U);
        break;
      }
      default: {
        double U = Xi(0)/double(Size-1);
        double V = Xi(1)/double(Size-1);
        double W = Xi(2)/double(Size-1);
        Coords(0) = (1.-U)*(-1.) + U + WaveAmplitude*std::sin(2.*PI*V)*std::sin(2.*PI*W);
        Coords(1) = (1.-V)*(-1.) + V + WaveAmplitude*std::sin(2.*PI*W)*std::sin(2.*PI*U);
        Coords(2) = (1.-W)*(-1.) + W + WaveAmplitude*std::sin(2.*PI*U)*std::sin(2.*PI*V);
        break;
      }}
      return Coords;
    };

    if (BackgroundIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      auto GeometryEditHandle = GeometryComponent.EditGeometry(1);
      ovk::geometry &Geometry = *GeometryEditHandle;
      auto CoordsEditHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsEditHandle;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ovk::tuple<double> Xi = {double(i),double(j),double(k)};
            ovk::tuple<double> CoordTuple = CoordFunc(Xi);
            for (int iDim = 0; iDim < NumDims; ++iDim) {
              Coords(iDim)(i,j,k) = CoordTuple(iDim);
            }
          }
        }
      }
    }

    if (ForegroundIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      auto GeometryEditHandle = GeometryComponent.EditGeometry(2);
      ovk::geometry &Geometry = *GeometryEditHandle;
      auto CoordsEditHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsEditHandle;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            ovk::tuple<int> Point = {i,j,k};
            ovk::tuple<double> Xi = {0.,0.,0.};
            for (int iDim = 0; iDim < NumDims; ++iDim) {
              Xi(iDim) = 0.25*double(Size-1) + 0.5*double(Point(iDim));
            }
            ovk::tuple<double> CoordTuple = CoordFunc(Xi);
            for (int iDim = 0; iDim < NumDims; ++iDim) {
              Coords(iDim)(i,j,k) = CoordTuple(iDim);
            }
          }
        }
      }
    }

  }

  Domain.CreateComponent<ovk::state_component>(2);

  {

    auto StateComponentEditHandle = Domain.EditComponent<ovk::state_component>(2);
    ovk::state_component &StateComponent = *StateComponentEditHandle;

    StateComponent.CreateStates(GridIDs);

    if (BackgroundIsLocal && PreCutHole) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &ExtendedRange = Grid.ExtendedRange();
      auto StateEditHandle = StateComponent.EditState(1);
      ovk::state &State = *StateEditHandle;
      auto FlagsEditHandle = State.EditFlags();
      ovk::distributed_field<ovk::state_flags> &Flags = *FlagsEditHandle;
      ovk::range HoleRange = ovk::MakeEmptyRange(NumDims);
      int HoleSize = 2*(Size/4-1);
      for (int iDim = 0; iDim < NumDims; ++iDim) {
        HoleRange.Begin(iDim) = (Size-HoleSize)/2;
        HoleRange.End(iDim) = HoleRange.Begin(iDim) + HoleSize;
      }
      for (int k = ExtendedRange.Begin(2); k < ExtendedRange.End(2); ++k) {
        for (int j = ExtendedRange.Begin(1); j < ExtendedRange.End(1); ++j) {
          for (int i = ExtendedRange.Begin(0); i < ExtendedRange.End(0); ++i) {
            if (HoleRange.Contains({i,j,k})) {
              Flags(i,j,k) &= ~ovk::state_flags::ACTIVE;
            }
          }
        }
      }
    }

  }

  return Domain;

}

}
