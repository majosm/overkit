// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "tests/fixtures/Interface.hpp"

#include "support/Decomp.hpp"

#include <overkit.hpp>

namespace tests {

ovk::domain Interface2D(ovk::comm_view Comm, const ovk::box &Bounds, const ovk::tuple<int> &Size,
  const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage) {

  OVK_DEBUG_ASSERT(!Periodic(1), "Can't be periodic in interface-normal direction.");

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(Comm)
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

  if (LowerIsLocal) {
    ovk::tuple<int> GridSize = {Size(0), (Size(1)+2)/2, Size(2)};
    ovk::range LocalRange = support::CartesianDecomp(2, {GridSize}, LowerComm);
    MaybeGridParams(0) = Domain.MakeGridParams()
      .SetName("Lower")
      .SetComm(LowerComm)
      .SetGlobalRange({GridSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  if (UpperIsLocal) {
    ovk::tuple<int> GridSize = {Size(0), Size(1)+2-(Size(1)+2)/2, Size(2)};
    ovk::range LocalRange = support::CartesianDecomp(2, {GridSize}, UpperComm);
    MaybeGridParams(1) = Domain.MakeGridParams()
      .SetName("Upper")
      .SetComm(UpperComm)
      .SetGlobalRange({GridSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  Domain.CreateGrids(GridIDs, MaybeGridParams);

  return Domain;

}

ovk::domain Interface2DManualConnectivity(ovk::comm_view Comm, const ovk::box &Bounds, const
  ovk::tuple<int> &Size, const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage) {

  ovk::domain Domain = Interface2D(Comm, Bounds, Size, Periodic, PeriodicStorage);

  bool LowerIsLocal = Domain.GridIsLocal(1);
  bool UpperIsLocal = Domain.GridIsLocal(2);

  ovk::tuple<int> LowerSize = Domain.GridInfo(1).GlobalRange().Size();

  Domain.CreateComponent<ovk::connectivity_component>(1);

  auto ConnectivityComponentEditHandle = Domain.EditComponent<ovk::connectivity_component>(1);
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

  if (LowerIsLocal) {
    ovk::tuple<int> GridSize = {Size(0), Size(1), (Size(2)+2)/2};
    ovk::range LocalRange = support::CartesianDecomp(3, {GridSize}, LowerComm);
    MaybeGridParams(0) = Domain.MakeGridParams()
      .SetName("Lower")
      .SetComm(LowerComm)
      .SetGlobalRange({GridSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  if (UpperIsLocal) {
    ovk::tuple<int> GridSize = {Size(0), Size(1), Size(2)+2-(Size(2)+2)/2};
    ovk::range LocalRange = support::CartesianDecomp(3, {GridSize}, UpperComm);
    MaybeGridParams(1) = Domain.MakeGridParams()
      .SetName("Upper")
      .SetComm(UpperComm)
      .SetGlobalRange({GridSize})
      .SetLocalRange(LocalRange)
      .SetPeriodic(Periodic)
      .SetPeriodicStorage(PeriodicStorage);
  }

  Domain.CreateGrids(GridIDs, MaybeGridParams);

  return Domain;

}

ovk::domain Interface3DManualConnectivity(ovk::comm_view Comm, const ovk::box &Bounds, const
  ovk::tuple<int> &Size, const ovk::tuple<bool> &Periodic, ovk::periodic_storage PeriodicStorage) {

  ovk::domain Domain = Interface3D(Comm, Bounds, Size, Periodic, PeriodicStorage);

  bool LowerIsLocal = Domain.GridIsLocal(1);
  bool UpperIsLocal = Domain.GridIsLocal(2);

  ovk::tuple<int> LowerSize = Domain.GridInfo(1).GlobalRange().Size();

  Domain.CreateComponent<ovk::connectivity_component>(1);

  auto ConnectivityComponentEditHandle = Domain.EditComponent<ovk::connectivity_component>(1);
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
