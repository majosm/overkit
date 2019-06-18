// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <overkit.hpp>

#include "examples/Common.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <exception>
#include <memory>
#include <utility>
#include <vector>

using examples::CartesianDecomp;

#define Print(...) printf(__VA_ARGS__); fflush(stdout)

namespace {
void Interface();
}

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  try {
    Interface();
  } catch (const std::exception &Exception) {
    Print("Encountered error:\n%s\n", Exception.what());
  } catch (...) {
    Print("Unknown error occurred.\n");
  }

  MPI_Finalize();

  return 0;

}

namespace {

struct grid_data {
  MPI_Comm Comm = MPI_COMM_NULL;
  std::array<int,3> Size = {{0,0,1}};
  std::array<int,6> LocalRange = {{0,0,0,0,0,1}};
  long long NumLocalPoints = 0;
  std::array<int,6> ExtendedRange = {{0,0,0,0,0,1}};
  long long NumExtendedPoints = 0;
  ~grid_data() noexcept {
    if (Comm != MPI_COMM_NULL) {
      MPI_Comm_free(&Comm);
    }
  }
};

void Interface() {

  int NumWorldProcs, WorldRank;
  MPI_Comm_size(MPI_COMM_WORLD, &NumWorldProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

  auto Context = std::make_shared<ovk::context>(ovk::CreateContext(ovk::context::params()
    .SetComm(MPI_COMM_WORLD)
    .SetLogLevel(ovk::log_level::ERRORS | ovk::log_level::WARNINGS | ovk::log_level::STATUS)
  ));

  ovk::domain Domain = ovk::CreateDomain(Context, ovk::domain::params()
    .SetDimension(2)
    .SetComm(MPI_COMM_WORLD)
  );

  std::array<int,3> Size = {{64,64,1}};

  std::array<int,2> GridIDs = {{1, 2}};

  bool Grid1IsLocal = WorldRank < std::max(NumWorldProcs/2, 1);
  bool Grid2IsLocal = WorldRank >= NumWorldProcs/2;

  std::array<int,3> Grid1Size = {{(Size[0]+2)/2, Size[1], Size[2]}};
  std::array<int,3> Grid2Size = {{Size[0]+2-(Size[0]+2)/2, Size[1], Size[2]}};

  grid_data Grid1Data, Grid2Data;

  if (Grid1IsLocal) {
    grid_data &Data = Grid1Data;
    MPI_Comm TempComm;
    MPI_Comm_split(MPI_COMM_WORLD, 0, WorldRank, &TempComm);
    int NumGridProcs;
    MPI_Comm_size(TempComm, &NumGridProcs);
    std::array<int,3> CartDims = {{0,0,1}};
    std::array<int,3> CartPeriods = {{0,0,0}};
    MPI_Dims_create(NumGridProcs, 2, CartDims.data());
    MPI_Cart_create(TempComm, 2, CartDims.data(), CartPeriods.data(), 1, &Data.Comm);
    MPI_Comm_free(&TempComm);
    Data.Size = Grid1Size;
    Data.LocalRange = CartesianDecomp(2, Data.Size, Data.Comm);
    Data.NumLocalPoints =
      (long long)(Data.LocalRange[3] - Data.LocalRange[0]) *
      (long long)(Data.LocalRange[4] - Data.LocalRange[1]) *
      (long long)(Data.LocalRange[5] - Data.LocalRange[2]);
    // Pretend we have a halo
    Data.ExtendedRange = Data.LocalRange;
    for (int iDim = 0; iDim < 2; ++iDim) {
      if (Data.LocalRange[iDim] > 0) --Data.ExtendedRange[iDim];
      if (Data.LocalRange[3+iDim] < Data.Size[iDim]) ++Data.ExtendedRange[3+iDim];
    }
    Data.NumExtendedPoints =
      (long long)(Data.ExtendedRange[3] - Data.ExtendedRange[0]) *
      (long long)(Data.ExtendedRange[4] - Data.ExtendedRange[1]) *
      (long long)(Data.ExtendedRange[5] - Data.ExtendedRange[2]);
  } else {
    MPI_Comm DummyComm;
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, WorldRank, &DummyComm);
  }

  if (Grid2IsLocal) {
    grid_data &Data = Grid2Data;
    MPI_Comm TempComm;
    MPI_Comm_split(MPI_COMM_WORLD, 0, WorldRank, &TempComm);
    int NumGridProcs;
    MPI_Comm_size(TempComm, &NumGridProcs);
    std::array<int,3> CartDims = {{0,0,1}};
    std::array<int,3> CartPeriods = {{0,0,0}};
    MPI_Dims_create(NumGridProcs, 2, CartDims.data());
    MPI_Cart_create(TempComm, 2, CartDims.data(), CartPeriods.data(), 1, &Data.Comm);
    MPI_Comm_free(&TempComm);
    Data.Size = Grid2Size;
    Data.LocalRange = CartesianDecomp(2, Data.Size, Data.Comm);
    Data.NumLocalPoints =
      (long long)(Data.LocalRange[3] - Data.LocalRange[0]) *
      (long long)(Data.LocalRange[4] - Data.LocalRange[1]) *
      (long long)(Data.LocalRange[5] - Data.LocalRange[2]);
    // Pretend we have a halo
    Data.ExtendedRange = Data.LocalRange;
    for (int iDim = 0; iDim < 2; ++iDim) {
      if (Data.LocalRange[iDim] > 0) --Data.ExtendedRange[iDim];
      if (Data.LocalRange[3+iDim] < Data.Size[iDim]) ++Data.ExtendedRange[3+iDim];
    }
    Data.NumExtendedPoints =
      (long long)(Data.ExtendedRange[3] - Data.ExtendedRange[0]) *
      (long long)(Data.ExtendedRange[4] - Data.ExtendedRange[1]) *
      (long long)(Data.ExtendedRange[5] - Data.ExtendedRange[2]);
  } else {
    MPI_Comm DummyComm;
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, WorldRank, &DummyComm);
  }

  std::array<ovk::optional<ovk::grid::params>,2> MaybeGridParams;

  if (Grid1IsLocal) {
    MaybeGridParams[0] = Domain.MakeGridParams()
      .SetName("Left")
      .SetComm(Grid1Data.Comm)
      .SetGlobalRange({Grid1Data.Size})
      .SetLocalRange({&Grid1Data.LocalRange[0], &Grid1Data.LocalRange[3]})
      .SetGeometryType(ovk::geometry_type::UNIFORM);
  }

  if (Grid2IsLocal) {
    MaybeGridParams[1] = Domain.MakeGridParams()
      .SetName("Right")
      .SetComm(Grid2Data.Comm)
      .SetGlobalRange({Grid2Data.Size})
      .SetLocalRange({&Grid2Data.LocalRange[0], &Grid2Data.LocalRange[3]})
      .SetGeometryType(ovk::geometry_type::UNIFORM);
  }

  Domain.CreateGrids(GridIDs, MaybeGridParams);

  constexpr int CONNECTIVITY_ID = 1;
  Domain.CreateComponent<ovk::connectivity_component>(CONNECTIVITY_ID);

  {

    auto ConnectivityComponentEditHandle = Domain.EditComponent<ovk::connectivity_component>(
      CONNECTIVITY_ID);
    ovk::connectivity_component &ConnectivityComponent = *ConnectivityComponentEditHandle;

    std::array<int,2> MGridIDs = {{1, 2}};
    std::array<int,2> NGridIDs = {{2, 1}};

    ConnectivityComponent.CreateConnectivities(MGridIDs, NGridIDs);

    if (Grid1IsLocal) {

      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.End(0) == GlobalRange.End(0);

      auto ConnectivityMEditHandle = ConnectivityComponent.EditConnectivityM(1, 2);
      ovk::connectivity_m &ConnectivityM = *ConnectivityMEditHandle;

      long long NumDonors = HasInterface ? LocalRange.Size(1) : 0;
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
          Extents(0,0,iDonor) = GlobalRange.End(0)-2;
          Extents(0,1,iDonor) = j;
          Extents(1,0,iDonor) = Extents(0,0,iDonor)+1;
          Extents(1,1,iDonor) = Extents(0,1,iDonor)+1;
          Coords(0,iDonor) = 0.;
          Coords(1,iDonor) = 0.;
          InterpCoefs(0,0,iDonor) = 1.;
          InterpCoefs(1,0,iDonor) = 1.;
          Destinations(0,iDonor) = 0;
          Destinations(1,iDonor) = j;
          ++iDonor;
        }
      }

      auto ConnectivityNEditHandle = ConnectivityComponent.EditConnectivityN(2, 1);
      ovk::connectivity_n &ConnectivityN = *ConnectivityNEditHandle;

      long long NumReceivers = HasInterface ? LocalRange.Size(1) : 0;
      ConnectivityN.Resize(NumReceivers);

      auto PointsEditHandle = ConnectivityN.EditPoints();
      auto SourcesEditHandle = ConnectivityN.EditSources();

      if (HasInterface) {
        ovk::array<int,2> &Points = *PointsEditHandle;
        ovk::array<int,2> &Sources = *SourcesEditHandle;
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          Points(0,iReceiver) = GlobalRange.End(0)-1;
          Points(1,iReceiver) = j;
          Sources(0,iReceiver) = 1;
          Sources(1,iReceiver) = j;
          ++iReceiver;
        }
      }

    }

    if (Grid2IsLocal) {

      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.Begin(0) == GlobalRange.Begin(0);

      auto ConnectivityMEditHandle = ConnectivityComponent.EditConnectivityM(2, 1);
      ovk::connectivity_m &ConnectivityM = *ConnectivityMEditHandle;

      long long NumDonors = HasInterface ? LocalRange.Size(1) : 0;
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
          Extents(0,0,iDonor) = 1;
          Extents(0,1,iDonor) = j;
          Extents(1,0,iDonor) = Extents(0,0,iDonor)+1;
          Extents(1,1,iDonor) = Extents(0,1,iDonor)+1;
          Coords(0,iDonor) = 0.;
          Coords(1,iDonor) = 0.;
          InterpCoefs(0,0,iDonor) = 1.;
          InterpCoefs(1,0,iDonor) = 1.;
          Destinations(0,iDonor) = Grid2Size[0]-1;
          Destinations(1,iDonor) = j;
          ++iDonor;
        }
      }

      auto ConnectivityNEditHandle = ConnectivityComponent.EditConnectivityN(1, 2);
      ovk::connectivity_n &ConnectivityN = *ConnectivityNEditHandle;

      long long NumReceivers = HasInterface ? LocalRange.Size(1) : 0;
      ConnectivityN.Resize(NumReceivers);

      auto PointsEditHandle = ConnectivityN.EditPoints();
      auto SourcesEditHandle = ConnectivityN.EditSources();

      if (HasInterface) {
        ovk::array<int,2> &Points = *PointsEditHandle;
        ovk::array<int,2> &Sources = *SourcesEditHandle;
        long long iReceiver = 0;
        for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
          Points(0,iReceiver) = 0;
          Points(1,iReceiver) = j;
          Sources(0,iReceiver) = Grid1Size[0]-2;
          Sources(1,iReceiver) = j;
          ++iReceiver;
        }
      }

    }

  }

  ovk::exchanger Exchanger = ovk::CreateExchanger(Context);

  Exchanger.Bind(Domain, ovk::exchanger::bindings()
    .SetConnectivityComponentID(CONNECTIVITY_ID)
  );

  auto &ConnectivityComponent = Domain.Component<ovk::connectivity_component>(CONNECTIVITY_ID);

  std::vector<double> Grid1DonorValues, Grid1ReceiverValues;
  if (Grid1IsLocal) {
    ovk::range ExtendedRange(&Grid1Data.ExtendedRange[0], &Grid1Data.ExtendedRange[3]);
    const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(1,2);
    Exchanger.CreateCollect(1, 2, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
      ExtendedRange, ovk::array_layout::ROW_MAJOR);
    Exchanger.CreateSend(1, 2, 1, ovk::data_type::DOUBLE, 1, 1);
    Grid1DonorValues.resize(ConnectivityM.Count());
    const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(2,1);
    Exchanger.CreateReceive(2, 1, 1, ovk::data_type::DOUBLE, 1, 1);
    Exchanger.CreateDisperse(2, 1, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
      ExtendedRange, ovk::array_layout::ROW_MAJOR);
    Grid1ReceiverValues.resize(ConnectivityN.Count());
  }

  std::vector<double> Grid2DonorValues, Grid2ReceiverValues;
  if (Grid2IsLocal) {
    ovk::range ExtendedRange(&Grid2Data.ExtendedRange[0], &Grid2Data.ExtendedRange[3]);
    const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM(2,1);
    Exchanger.CreateCollect(2, 1, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
      ExtendedRange, ovk::array_layout::ROW_MAJOR);
    Exchanger.CreateSend(2, 1, 1, ovk::data_type::DOUBLE, 1, 1);
    Grid2DonorValues.resize(ConnectivityM.Count());
    const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(1,2);
    Exchanger.CreateReceive(1, 2, 1, ovk::data_type::DOUBLE, 1, 1);
    Exchanger.CreateDisperse(1, 2, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
      ExtendedRange, ovk::array_layout::ROW_MAJOR);
    Grid2ReceiverValues.resize(ConnectivityN.Count());
  }

  std::vector<double> Grid1FieldValues;
  if (Grid1IsLocal) {
    Grid1FieldValues.resize(Grid1Data.NumExtendedPoints, -1.);
  }

  std::vector<double> Grid2FieldValues;
  if (Grid2IsLocal) {
    Grid2FieldValues.resize(Grid2Data.NumExtendedPoints, 1.);
  }

  std::vector<ovk::request> Requests;

  if (Grid1IsLocal) {
    double *ReceiverValues = Grid1ReceiverValues.data();
    ovk::request Request = Exchanger.Receive(2, 1, 1, &ReceiverValues);
    Requests.push_back(std::move(Request));
  }

  if (Grid2IsLocal) {
    double *ReceiverValues = Grid2ReceiverValues.data();
    ovk::request Request = Exchanger.Receive(1, 2, 1, &ReceiverValues);
    Requests.push_back(std::move(Request));
  }

  if (Grid1IsLocal) {
    const double *FieldValues = Grid1FieldValues.data();
    double *DonorValues = Grid1DonorValues.data();
    Exchanger.Collect(1, 2, 1, &FieldValues, &DonorValues);
    ovk::request Request = Exchanger.Send(1, 2, 1, &DonorValues);
    Requests.push_back(std::move(Request));
  }

  if (Grid2IsLocal) {
    const double *FieldValues = Grid2FieldValues.data();
    double *DonorValues = Grid2DonorValues.data();
    Exchanger.Collect(2, 1, 1, &FieldValues, &DonorValues);
    ovk::request Request = Exchanger.Send(2, 1, 1, &DonorValues);
    Requests.push_back(std::move(Request));
  }

  ovk::WaitAll(Requests);

  if (Grid1IsLocal) {
    const double *ReceiverValues = Grid1ReceiverValues.data();
    double *FieldValues = Grid1FieldValues.data();
    Exchanger.Disperse(2, 1, 1, &ReceiverValues, &FieldValues);
  }

  if (Grid2IsLocal) {
    const double *ReceiverValues = Grid2ReceiverValues.data();
    double *FieldValues = Grid2FieldValues.data();
    Exchanger.Disperse(1, 2, 1, &ReceiverValues, &FieldValues);
  }

}

}
