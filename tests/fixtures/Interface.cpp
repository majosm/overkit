// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "tests/fixtures/Interface.hpp"

#include "support/Decomp.hpp"
#include "support/XDMF.hpp"

#include <overkit.hpp>

namespace tests {

ovk::domain Interface2D(ovk::comm_view Comm, const ovk::box &Bounds, const ovk::tuple<int> &Size,
  const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(!Periodic(1), "Can't be periodic in interface-normal direction.");

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(Comm)
    .SetStatusLoggingThreshold(0)
  ));

  ovk::domain Domain = ovk::CreateDomain(std::move(Context), ovk::domain::params()
    .SetDimension(2)
    .SetComm(Comm)
  );

  bool LowerIsLocal = Comm.Rank() < ovk::Max(Comm.Size()/2, 1);
  bool UpperIsLocal = Comm.Rank() >= Comm.Size()/2;

  ovk::comm LowerComm = ovk::CreateSubsetComm(Comm, LowerIsLocal);
  if (LowerIsLocal) {
    LowerComm = ovk::CreateCartComm(LowerComm, 2, {0,0,1}, Periodic);
  }

  ovk::comm UpperComm = ovk::CreateSubsetComm(Comm, UpperIsLocal);
  if (UpperIsLocal) {
    UpperComm = ovk::CreateCartComm(UpperComm, 2, {0,0,1}, Periodic);
  }

  ovk::array<int> GridIDs({2}, {1, 2});
  ovk::array<ovk::optional<ovk::grid::params>> MaybeGridParams({2});

  ovk::tuple<int> LowerSize = {Size(0), (Size(1)+2)/2, Size(2)};
  ovk::tuple<int> UpperSize = {Size(0), Size(1)+2-(Size(1)+2)/2, Size(2)};

  if (LowerIsLocal) {
    ovk::range LocalRange = support::CartesianDecomp(2, {LowerSize}, LowerComm);
    MaybeGridParams(0) = Domain.MakeGridParams()
      .SetName("Lower")
      .SetComm(LowerComm)
      .SetGlobalRange({LowerSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  if (UpperIsLocal) {
    ovk::range LocalRange = support::CartesianDecomp(2, {UpperSize}, UpperComm);
    MaybeGridParams(1) = Domain.MakeGridParams()
      .SetName("Upper")
      .SetComm(UpperComm)
      .SetGlobalRange({UpperSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  Domain.CreateGrids(GridIDs, MaybeGridParams);

  Domain.CreateComponent<ovk::geometry_component>(1);

  {

    auto GeometryComponentHandle = Domain.EditComponent<ovk::geometry_component>(1);
    ovk::geometry_component &GeometryComponent = *GeometryComponentHandle;

    ovk::array<ovk::optional<ovk::geometry::params>> MaybeGeometryParams({2});

    if (LowerIsLocal) {
      MaybeGeometryParams(0) = ovk::geometry::params()
        .SetType(ovk::geometry_type::UNIFORM);
    }

    if (UpperIsLocal) {
      MaybeGeometryParams(1) = ovk::geometry::params()
        .SetType(ovk::geometry_type::UNIFORM);
    }

    GeometryComponent.CreateGeometries(GridIDs, MaybeGeometryParams);

    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      auto GeometryHandle = GeometryComponent.EditGeometry(1);
      ovk::geometry &Geometry = *GeometryHandle;
      auto CoordsHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsHandle;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i)/double(LowerSize(0)-1);
          double V = double(j)/double(Size(1)-1);
          Coords(0)(i,j,0) = 2.*(U-0.5);
          Coords(1)(i,j,0) = 2.*(V-0.5);
        }
      }
    }

    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      auto GeometryHandle = GeometryComponent.EditGeometry(2);
      ovk::geometry &Geometry = *GeometryHandle;
      auto CoordsHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsHandle;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          double U = double(i)/double(UpperSize(0)-1);
          double V = double(j+LowerSize(1)-2)/double(Size(1)-1);
          Coords(0)(i,j,0) = 2.*(U-0.5);
          Coords(1)(i,j,0) = 2.*(V-0.5);
        }
      }
    }

  }

// #ifdef OVK_HAVE_XDMF
//   auto &GeometryComponent = Domain.Component<ovk::geometry_component>(1);

//   ovk::elem<support::xdmf_grid_meta,2> XDMFGrids = {
//     {"Lower", LowerSize},
//     {"Upper", UpperSize}
//   };

//   support::CreateXDMF("Interface2D.xmf", 2, Comm, std::move(XDMFGrids), {});

//   if (LowerIsLocal) {
//     const ovk::grid &Grid = Domain.Grid(1);
//     auto &Coords = GeometryComponent.Geometry(1).Coords();
//     support::xdmf XDMF = support::OpenXDMF("Interface2D.xmf", Grid.Comm());
//     for (int iDim = 0; iDim < 2; ++iDim) {
//       XDMF.WriteGeometry("Lower", iDim, Coords(iDim), Grid.LocalRange());
//     }
//   }

//   if (UpperIsLocal) {
//     const ovk::grid &Grid = Domain.Grid(2);
//     auto &Coords = GeometryComponent.Geometry(2).Coords();
//     support::xdmf XDMF = support::OpenXDMF("Interface2D.xmf", Grid.Comm());
//     for (int iDim = 0; iDim < 2; ++iDim) {
//       XDMF.WriteGeometry("Upper", iDim, Coords(iDim), Grid.LocalRange());
//     }
//   }
// #endif

  return Domain;

}

ovk::domain Interface2DManualConnectivity(ovk::comm_view Comm, const ovk::box &Bounds, const
  ovk::tuple<int> &Size, const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage) {

  ovk::domain Domain = Interface2D(Comm, Bounds, Size, Periodic, PeriodicStorage);

  bool LowerIsLocal = Domain.GridIsLocal(1);
  bool UpperIsLocal = Domain.GridIsLocal(2);

  ovk::tuple<int> LowerSize = Domain.GridInfo(1).GlobalRange().Size();

  Domain.CreateComponent<ovk::connectivity_component>(4);

  auto ConnectivityComponentEditHandle = Domain.EditComponent<ovk::connectivity_component>(4);
  ovk::connectivity_component &ConnectivityComponent = *ConnectivityComponentEditHandle;

  ovk::array<ovk::elem<int,2>> ConnectivityIDs({2}, {{1,2}, {2,1}});

  ConnectivityComponent.CreateConnectivities(ConnectivityIDs);

  if (LowerIsLocal) {

    const ovk::grid &Grid = Domain.Grid(1);
    const ovk::range &LocalRange = Grid.LocalRange();

    bool HasInterface = LocalRange.End(1) == LowerSize(1);

    auto ConnectivityMEditHandle = ConnectivityComponent.EditConnectivityM({1,2});
    ovk::connectivity_m &ConnectivityM = *ConnectivityMEditHandle;

    long long NumDonors = HasInterface ? LocalRange.Size(0) : 0;
    ConnectivityM.Resize(NumDonors, 1);

    auto ExtentsEditHandle = ConnectivityM.EditExtents();
    auto CoordsEditHandle = ConnectivityM.EditCoords();
    auto InterpCoefsEditHandle = ConnectivityM.EditInterpCoefs();
    auto DestinationsEditHandle = ConnectivityM.EditDestinations();

    if (HasInterface) {
      ovk::array<int,3> &Extents = *ExtentsEditHandle;
      ovk::array<double,2> &Coords = *CoordsEditHandle;
      ovk::array<double,3> &InterpCoefs = *InterpCoefsEditHandle;
      ovk::array<int,2> &Destinations = *DestinationsEditHandle;
      long long iDonor = 0;
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        Extents(0,0,iDonor) = i;
        Extents(0,1,iDonor) = LowerSize(1)-2;
        Extents(1,0,iDonor) = Extents(0,0,iDonor)+1;
        Extents(1,1,iDonor) = Extents(0,1,iDonor)+1;
        Coords(0,iDonor) = 0.;
        Coords(1,iDonor) = 0.;
        InterpCoefs(0,0,iDonor) = 1.;
        InterpCoefs(1,0,iDonor) = 1.;
        Destinations(0,iDonor) = i;
        Destinations(1,iDonor) = 0;
        ++iDonor;
      }
    }

    auto ConnectivityNEditHandle = ConnectivityComponent.EditConnectivityN({2,1});
    ovk::connectivity_n &ConnectivityN = *ConnectivityNEditHandle;

    long long NumReceivers = HasInterface ? LocalRange.Size(0) : 0;
    ConnectivityN.Resize(NumReceivers);

    auto PointsEditHandle = ConnectivityN.EditPoints();
    auto SourcesEditHandle = ConnectivityN.EditSources();

    if (HasInterface) {
      ovk::array<int,2> &Points = *PointsEditHandle;
      ovk::array<int,2> &Sources = *SourcesEditHandle;
      long long iReceiver = 0;
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        Points(0,iReceiver) = i; 
        Points(1,iReceiver) = LowerSize(1)-1;
        Sources(0,iReceiver) = i;
        Sources(1,iReceiver) = 1;
        ++iReceiver;
      }
    }

  }

  if (UpperIsLocal) {

    const ovk::grid &Grid = Domain.Grid(2);
    const ovk::range &LocalRange = Grid.LocalRange();

    bool HasInterface = LocalRange.Begin(1) == 0;

    auto ConnectivityMEditHandle = ConnectivityComponent.EditConnectivityM({2,1});
    ovk::connectivity_m &ConnectivityM = *ConnectivityMEditHandle;

    long long NumDonors = HasInterface ? LocalRange.Size(0) : 0;
    ConnectivityM.Resize(NumDonors, 1);

    auto ExtentsEditHandle = ConnectivityM.EditExtents();
    auto CoordsEditHandle = ConnectivityM.EditCoords();
    auto InterpCoefsEditHandle = ConnectivityM.EditInterpCoefs();
    auto DestinationsEditHandle = ConnectivityM.EditDestinations();

    if (HasInterface) {
      ovk::array<int,3> &Extents = *ExtentsEditHandle;
      ovk::array<double,2> &Coords = *CoordsEditHandle;
      ovk::array<double,3> &InterpCoefs = *InterpCoefsEditHandle;
      ovk::array<int,2> &Destinations = *DestinationsEditHandle;
      long long iDonor = 0;
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        Extents(0,0,iDonor) = i;
        Extents(0,1,iDonor) = 1;
        Extents(1,0,iDonor) = Extents(0,0,iDonor)+1;
        Extents(1,1,iDonor) = Extents(0,1,iDonor)+1;
        Coords(0,iDonor) = 0.;
        Coords(1,iDonor) = 0.;
        InterpCoefs(0,0,iDonor) = 1.;
        InterpCoefs(1,0,iDonor) = 1.;
        Destinations(0,iDonor) = i;
        Destinations(1,iDonor) = LowerSize(1)-1;
        ++iDonor;
      }
    }

    auto ConnectivityNEditHandle = ConnectivityComponent.EditConnectivityN({1,2});
    ovk::connectivity_n &ConnectivityN = *ConnectivityNEditHandle;

    long long NumReceivers = HasInterface ? LocalRange.Size(0) : 0;
    ConnectivityN.Resize(NumReceivers);

    auto PointsEditHandle = ConnectivityN.EditPoints();
    auto SourcesEditHandle = ConnectivityN.EditSources();

    if (HasInterface) {
      ovk::array<int,2> &Points = *PointsEditHandle;
      ovk::array<int,2> &Sources = *SourcesEditHandle;
      long long iReceiver = 0;
      for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
        Points(0,iReceiver) = i;
        Points(1,iReceiver) = 0;
        Sources(0,iReceiver) = i;
        Sources(1,iReceiver) = LowerSize(1)-2;
        ++iReceiver;
      }
    }

  }

  return Domain;

}

ovk::domain Interface3D(ovk::comm_view Comm, const ovk::box &Bounds, const ovk::tuple<int> &Size,
  const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(!Periodic(2), "Can't be periodic in interface-normal direction.");

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(Comm)
    .SetStatusLoggingThreshold(0)
  ));

  ovk::domain Domain = ovk::CreateDomain(std::move(Context), ovk::domain::params()
    .SetDimension(3)
    .SetComm(Comm)
  );

  bool LowerIsLocal = Comm.Rank() < ovk::Max(Comm.Size()/2, 1);
  bool UpperIsLocal = Comm.Rank() >= Comm.Size()/2;

  ovk::comm LowerComm = ovk::CreateSubsetComm(Comm, LowerIsLocal);
  if (LowerIsLocal) {
    LowerComm = ovk::CreateCartComm(LowerComm, 3, {0,0,0}, Periodic);
  }

  ovk::comm UpperComm = ovk::CreateSubsetComm(Comm, UpperIsLocal);
  if (UpperIsLocal) {
    UpperComm = ovk::CreateCartComm(UpperComm, 3, {0,0,0}, Periodic);
  }

  ovk::array<int> GridIDs({2}, {1, 2});
  ovk::array<ovk::optional<ovk::grid::params>> MaybeGridParams({2});

  ovk::tuple<int> LowerSize = {Size(0), Size(1), (Size(2)+2)/2};
  ovk::tuple<int> UpperSize = {Size(0), Size(1), Size(2)+2-(Size(2)+2)/2};

  if (LowerIsLocal) {
    ovk::range LocalRange = support::CartesianDecomp(3, {LowerSize}, LowerComm);
    MaybeGridParams(0) = Domain.MakeGridParams()
      .SetName("Lower")
      .SetComm(LowerComm)
      .SetGlobalRange({LowerSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  if (UpperIsLocal) {
    ovk::range LocalRange = support::CartesianDecomp(3, {UpperSize}, UpperComm);
    MaybeGridParams(1) = Domain.MakeGridParams()
      .SetName("Upper")
      .SetComm(UpperComm)
      .SetGlobalRange({UpperSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  Domain.CreateGrids(GridIDs, MaybeGridParams);

  Domain.CreateComponent<ovk::geometry_component>(1);

  {

    auto GeometryComponentHandle = Domain.EditComponent<ovk::geometry_component>(1);
    ovk::geometry_component &GeometryComponent = *GeometryComponentHandle;

    ovk::array<ovk::optional<ovk::geometry::params>> MaybeGeometryParams({2});

    if (LowerIsLocal) {
      MaybeGeometryParams(0) = ovk::geometry::params()
        .SetType(ovk::geometry_type::UNIFORM);
    }

    if (UpperIsLocal) {
      MaybeGeometryParams(1) = ovk::geometry::params()
        .SetType(ovk::geometry_type::UNIFORM);
    }

    GeometryComponent.CreateGeometries(GridIDs, MaybeGeometryParams);

    if (LowerIsLocal) {
      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &LocalRange = Grid.LocalRange();
      auto GeometryHandle = GeometryComponent.EditGeometry(1);
      ovk::geometry &Geometry = *GeometryHandle;
      auto CoordsHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsHandle;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i)/double(LowerSize(0)-1);
            double V = double(j)/double(LowerSize(1)-1);
            double W = double(k)/double(Size(2)-1);
            Coords(0)(i,j,k) = 2.*(U-0.5);
            Coords(1)(i,j,k) = 2.*(V-0.5);
            Coords(2)(i,j,k) = 2.*(W-0.5);
          }
        }
      }
    }

    if (UpperIsLocal) {
      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &LocalRange = Grid.LocalRange();
      auto GeometryHandle = GeometryComponent.EditGeometry(2);
      ovk::geometry &Geometry = *GeometryHandle;
      auto CoordsHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsHandle;
      for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
            double U = double(i)/double(UpperSize(0)-1);
            double V = double(j)/double(UpperSize(1)-1);
            double W = double(k+LowerSize(2)-2)/double(Size(2)-1);
            Coords(0)(i,j,k) = 2.*(U-0.5);
            Coords(1)(i,j,k) = 2.*(V-0.5);
            Coords(2)(i,j,k) = 2.*(W-0.5);
          }
        }
      }
    }

  }

// #ifdef OVK_HAVE_XDMF
//   auto &GeometryComponent = Domain.Component<ovk::geometry_component>(1);

//   ovk::elem<support::xdmf_grid_meta,2> XDMFGrids = {
//     {"Lower", LowerSize},
//     {"Upper", UpperSize}
//   };

//   support::CreateXDMF("Interface3D.xmf", 3, Comm, std::move(XDMFGrids), {});

//   if (LowerIsLocal) {
//     const ovk::grid &Grid = Domain.Grid(1);
//     auto &Coords = GeometryComponent.Geometry(1).Coords();
//     support::xdmf XDMF = support::OpenXDMF("Interface3D.xmf", Grid.Comm());
//     for (int iDim = 0; iDim < 3; ++iDim) {
//       XDMF.WriteGeometry("Lower", iDim, Coords(iDim), Grid.LocalRange());
//     }
//   }

//   if (UpperIsLocal) {
//     const ovk::grid &Grid = Domain.Grid(2);
//     auto &Coords = GeometryComponent.Geometry(2).Coords();
//     support::xdmf XDMF = support::OpenXDMF("Interface3D.xmf", Grid.Comm());
//     for (int iDim = 0; iDim < 3; ++iDim) {
//       XDMF.WriteGeometry("Upper", iDim, Coords(iDim), Grid.LocalRange());
//     }
//   }
// #endif

  return Domain;

}

ovk::domain Interface3DManualConnectivity(ovk::comm_view Comm, const ovk::box &Bounds, const
  ovk::tuple<int> &Size, const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage) {

  ovk::domain Domain = Interface3D(Comm, Bounds, Size, Periodic, PeriodicStorage);

  bool LowerIsLocal = Domain.GridIsLocal(1);
  bool UpperIsLocal = Domain.GridIsLocal(2);

  ovk::tuple<int> LowerSize = Domain.GridInfo(1).GlobalRange().Size();

  Domain.CreateComponent<ovk::connectivity_component>(4);

  auto ConnectivityComponentEditHandle = Domain.EditComponent<ovk::connectivity_component>(4);
  ovk::connectivity_component &ConnectivityComponent = *ConnectivityComponentEditHandle;

  ovk::array<ovk::elem<int,2>> ConnectivityIDs({2}, {{1,2}, {2,1}});

  ConnectivityComponent.CreateConnectivities(ConnectivityIDs);

  if (LowerIsLocal) {

    const ovk::grid &Grid = Domain.Grid(1);
    const ovk::range &LocalRange = Grid.LocalRange();

    bool HasInterface = LocalRange.End(2) == LowerSize(2);

    auto ConnectivityMEditHandle = ConnectivityComponent.EditConnectivityM({1,2});
    ovk::connectivity_m &ConnectivityM = *ConnectivityMEditHandle;

    long long NumDonors = HasInterface ? LocalRange.Size(0)*LocalRange.Size(1) : 0;
    ConnectivityM.Resize(NumDonors, 1);

    auto ExtentsEditHandle = ConnectivityM.EditExtents();
    auto CoordsEditHandle = ConnectivityM.EditCoords();
    auto InterpCoefsEditHandle = ConnectivityM.EditInterpCoefs();
    auto DestinationsEditHandle = ConnectivityM.EditDestinations();

    if (HasInterface) {
      ovk::array<int,3> &Extents = *ExtentsEditHandle;
      ovk::array<double,2> &Coords = *CoordsEditHandle;
      ovk::array<double,3> &InterpCoefs = *InterpCoefsEditHandle;
      ovk::array<int,2> &Destinations = *DestinationsEditHandle;
      long long iDonor = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          Extents(0,0,iDonor) = i;
          Extents(0,1,iDonor) = j;
          Extents(0,2,iDonor) = LowerSize(2)-2;
          Extents(1,0,iDonor) = Extents(0,0,iDonor)+1;
          Extents(1,1,iDonor) = Extents(0,1,iDonor)+1;
          Extents(1,2,iDonor) = Extents(0,2,iDonor)+1;
          Coords(0,iDonor) = 0.;
          Coords(1,iDonor) = 0.;
          Coords(2,iDonor) = 0.;
          InterpCoefs(0,0,iDonor) = 1.;
          InterpCoefs(1,0,iDonor) = 1.;
          InterpCoefs(2,0,iDonor) = 1.;
          Destinations(0,iDonor) = i;
          Destinations(1,iDonor) = j;
          Destinations(2,iDonor) = 0;
          ++iDonor;
        }
      }
    }

    auto ConnectivityNEditHandle = ConnectivityComponent.EditConnectivityN({2,1});
    ovk::connectivity_n &ConnectivityN = *ConnectivityNEditHandle;

    long long NumReceivers = HasInterface ? LocalRange.Size(0)*LocalRange.Size(1) : 0;
    ConnectivityN.Resize(NumReceivers);

    auto PointsEditHandle = ConnectivityN.EditPoints();
    auto SourcesEditHandle = ConnectivityN.EditSources();

    if (HasInterface) {
      ovk::array<int,2> &Points = *PointsEditHandle;
      ovk::array<int,2> &Sources = *SourcesEditHandle;
      long long iReceiver = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          Points(0,iReceiver) = i;
          Points(1,iReceiver) = j;
          Points(2,iReceiver) = LowerSize(2)-1;
          Sources(0,iReceiver) = i;
          Sources(1,iReceiver) = j;
          Sources(2,iReceiver) = 1;
          ++iReceiver;
        }
      }
    }

  }

  if (UpperIsLocal) {

    const ovk::grid &Grid = Domain.Grid(2);
    const ovk::range &LocalRange = Grid.LocalRange();

    bool HasInterface = LocalRange.Begin(2) == 0;

    auto ConnectivityMEditHandle = ConnectivityComponent.EditConnectivityM({2,1});
    ovk::connectivity_m &ConnectivityM = *ConnectivityMEditHandle;

    long long NumDonors = HasInterface ? LocalRange.Size(0)*LocalRange.Size(1) : 0;
    ConnectivityM.Resize(NumDonors, 1);

    auto ExtentsEditHandle = ConnectivityM.EditExtents();
    auto CoordsEditHandle = ConnectivityM.EditCoords();
    auto InterpCoefsEditHandle = ConnectivityM.EditInterpCoefs();
    auto DestinationsEditHandle = ConnectivityM.EditDestinations();

    if (HasInterface) {
      ovk::array<int,3> &Extents = *ExtentsEditHandle;
      ovk::array<double,2> &Coords = *CoordsEditHandle;
      ovk::array<double,3> &InterpCoefs = *InterpCoefsEditHandle;
      ovk::array<int,2> &Destinations = *DestinationsEditHandle;
      long long iDonor = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          Extents(0,0,iDonor) = i;
          Extents(0,1,iDonor) = j;
          Extents(0,2,iDonor) = 1;
          Extents(1,0,iDonor) = Extents(0,0,iDonor)+1;
          Extents(1,1,iDonor) = Extents(0,1,iDonor)+1;
          Extents(1,2,iDonor) = Extents(0,2,iDonor)+1;
          Coords(0,iDonor) = 0.;
          Coords(1,iDonor) = 0.;
          Coords(2,iDonor) = 0.;
          InterpCoefs(0,0,iDonor) = 1.;
          InterpCoefs(1,0,iDonor) = 1.;
          InterpCoefs(2,0,iDonor) = 1.;
          Destinations(0,iDonor) = i;
          Destinations(1,iDonor) = j;
          Destinations(2,iDonor) = LowerSize(2)-1;
          ++iDonor;
        }
      }
    }

    auto ConnectivityNEditHandle = ConnectivityComponent.EditConnectivityN({1,2});
    ovk::connectivity_n &ConnectivityN = *ConnectivityNEditHandle;

    long long NumReceivers = HasInterface ? LocalRange.Size(0)*LocalRange.Size(1) : 0;
    ConnectivityN.Resize(NumReceivers);

    auto PointsEditHandle = ConnectivityN.EditPoints();
    auto SourcesEditHandle = ConnectivityN.EditSources();

    if (HasInterface) {
      ovk::array<int,2> &Points = *PointsEditHandle;
      ovk::array<int,2> &Sources = *SourcesEditHandle;
      long long iReceiver = 0;
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          Points(0,iReceiver) = i;
          Points(1,iReceiver) = j;
          Points(2,iReceiver) = 0;
          Sources(0,iReceiver) = i;
          Sources(1,iReceiver) = j;
          Sources(2,iReceiver) = LowerSize(2)-2;
          ++iReceiver;
        }
      }
    }

  }

  return Domain;

}

}
