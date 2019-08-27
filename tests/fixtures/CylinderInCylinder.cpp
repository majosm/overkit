// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "tests/fixtures/CylinderInCylinder.hpp"

#include "support/Constants.hpp"
#include "support/Decomp.hpp"

#include <overkit.hpp>

#include <mpi.h>

#include <cmath>

using support::PI;

namespace tests {

ovk::domain CylinderInCylinder(ovk::comm_view Comm, int Size, int OverlapAmount,
  bool Stagger, const ovk::tuple<bool> &DecompDirs, ovk::periodic_storage PeriodicStorage) {

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(Comm)
//     .SetLogLevel(ovk::log_level::ERRORS | ovk::log_level::WARNINGS | ovk::log_level::STATUS |
//       ovk::log_level::DEBUG)
  ));

  ovk::domain Domain = ovk::CreateDomain(std::move(Context), ovk::domain::params()
    .SetDimension(2)
    .SetComm(Comm)
  );

  ovk::tuple<int> InnerSize = {Size/2+1 + OverlapAmount-OverlapAmount/2, 2*Size, 1};
  ovk::tuple<int> OuterSize = {Size-Size/2+1 + OverlapAmount/2, 2*Size, 1};
  if (PeriodicStorage == ovk::periodic_storage::DUPLICATED) {
    ++InnerSize(1);
    ++OuterSize(1);
  }

  ovk::tuple<bool> Periodic = {false, true, false};

  bool InnerIsLocal = Comm.Rank() < ovk::Max(Comm.Size()/2, 1);
  bool OuterIsLocal = Comm.Rank() >= Comm.Size()/2;

  ovk::comm InnerComm = ovk::CreateSubsetComm(Comm, InnerIsLocal);
  if (InnerIsLocal) {
    ovk::tuple<int> CartDims = {int(!DecompDirs(0)),int(!DecompDirs(1)),1};
    InnerComm = ovk::CreateCartComm(InnerComm, 2, CartDims, Periodic);
  }

  ovk::comm OuterComm = ovk::CreateSubsetComm(Comm, OuterIsLocal);
  if (OuterIsLocal) {
    ovk::tuple<int> CartDims = {int(!DecompDirs(0)),int(!DecompDirs(1)),1};
    OuterComm = ovk::CreateCartComm(OuterComm, 2, CartDims, Periodic);
  }

  ovk::array<int> GridIDs({2}, {1, 2});
  ovk::array<ovk::optional<ovk::grid::params>> MaybeGridParams({2});

  if (InnerIsLocal) {
    ovk::range LocalRange = support::CartesianDecomp(2, {InnerSize}, InnerComm);
    MaybeGridParams(0) = Domain.MakeGridParams()
      .SetName("Inner")
      .SetComm(InnerComm)
      .SetGlobalRange({InnerSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  if (OuterIsLocal) {
    ovk::range LocalRange = support::CartesianDecomp(2, {OuterSize}, OuterComm);
    MaybeGridParams(1) = Domain.MakeGridParams()
      .SetName("Outer")
      .SetComm(OuterComm)
      .SetGlobalRange({OuterSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  Domain.CreateGrids(GridIDs, MaybeGridParams);

  Domain.CreateComponent<ovk::geometry_component>(1);

  {

    auto GeometryComponentEditHandle = Domain.EditComponent<ovk::geometry_component>(1);
    ovk::geometry_component &GeometryComponent = *GeometryComponentEditHandle;

    GeometryComponent.CreateGeometries(GridIDs);

    double RMin = 0.1;
    double RMax = 1.;

    double StaggerOffset = Stagger ? 0.5 : 0.;

    auto RadiusFunc = [&](double Xi) -> double {
      double XiMax = double(Size) + StaggerOffset;
      double U = Xi/XiMax;
      return RMin*std::pow(RMax/RMin, U);
    };

    auto ThetaFunc = [&](double Eta) -> double {
      double EtaMax = double(2*Size);
      double V = Eta/EtaMax;
      return 2.*PI*V;
    };

    if (InnerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      auto GeometryEditHandle = GeometryComponent.EditGeometry(1);
      ovk::geometry &Geometry = *GeometryEditHandle;
      auto CoordsEditHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsEditHandle;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double Radius = RadiusFunc(double(i));
          double Theta = ThetaFunc(double(j));
          Coords(0)(i,j,0) = Radius*std::cos(Theta);
          Coords(1)(i,j,0) = Radius*std::sin(Theta);
        }
      }
    }

    if (OuterIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      auto GeometryEditHandle = GeometryComponent.EditGeometry(2);
      ovk::geometry &Geometry = *GeometryEditHandle;
      auto CoordsEditHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsEditHandle;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double Radius = RadiusFunc(double(InnerSize(0))-StaggerOffset-OverlapAmount+double(i));
          double Theta = ThetaFunc(StaggerOffset+double(j));
          Coords(0)(i,j,0) = Radius*std::cos(Theta);
          Coords(1)(i,j,0) = Radius*std::sin(Theta);
        }
      }
    }

  }

  Domain.CreateComponent<ovk::state_component>(2);

  {
    auto StateComponentEditHandle = Domain.EditComponent<ovk::state_component>(2);
    ovk::state_component &StateComponent = *StateComponentEditHandle;
    StateComponent.CreateStates(GridIDs);
  }

  return Domain;

}

}
