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

using examples::CreateCartesianDecompDims;
using examples::CartesianDecomp;

namespace {
void Interface();
}

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  try {
    Interface();
  } catch (const std::exception &Exception) {
    std::printf("Encountered error:\n%s\n", Exception.what()); std::fflush(stdout);
  } catch (...) {
    std::printf("Unknown error occurred.\n"); std::fflush(stdout);
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
    .SetStatusLoggingThreshold(4)
  ));

  ovk::domain Domain = ovk::CreateDomain(Context, ovk::domain::params()
    .SetDimension(2)
    .SetComm(MPI_COMM_WORLD)
  );

  std::array<int,3> Size = {{64,64,1}};

  std::array<int,2> GridIDs = {{1, 2}};

  bool LeftIsLocal = WorldRank < std::max(NumWorldProcs/2, 1);
  bool RightIsLocal = WorldRank >= NumWorldProcs/2;

  std::array<int,3> LeftSize = {{(Size[0]+2)/2, Size[1], Size[2]}};
  std::array<int,3> RightSize = {{Size[0]+2-(Size[0]+2)/2, Size[1], Size[2]}};

  grid_data LeftData, RightData;

  if (LeftIsLocal) {
    grid_data &Data = LeftData;
    MPI_Comm TempComm;
    MPI_Comm_split(MPI_COMM_WORLD, 0, WorldRank, &TempComm);
    int NumGridProcs;
    MPI_Comm_size(TempComm, &NumGridProcs);
    std::array<int,3> CartDims = CreateCartesianDecompDims(NumGridProcs, 2, {{0,0,1}});
    std::array<int,3> CartPeriods = {{0,0,0}};
    MPI_Cart_create(TempComm, 2, CartDims.data(), CartPeriods.data(), 1, &Data.Comm);
    MPI_Comm_free(&TempComm);
    Data.Size = LeftSize;
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

  if (RightIsLocal) {
    grid_data &Data = RightData;
    MPI_Comm TempComm;
    MPI_Comm_split(MPI_COMM_WORLD, 0, WorldRank, &TempComm);
    int NumGridProcs;
    MPI_Comm_size(TempComm, &NumGridProcs);
    std::array<int,3> CartDims = CreateCartesianDecompDims(NumGridProcs, 2, {{0,0,1}});
    std::array<int,3> CartPeriods = {{0,0,0}};
    MPI_Cart_create(TempComm, 2, CartDims.data(), CartPeriods.data(), 1, &Data.Comm);
    MPI_Comm_free(&TempComm);
    Data.Size = RightSize;
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

  if (LeftIsLocal) {
    MaybeGridParams[0] = Domain.MakeGridParams()
      .SetName("Left")
      .SetComm(LeftData.Comm)
      .SetGlobalRange({LeftData.Size})
      .SetLocalRange({&LeftData.LocalRange[0], &LeftData.LocalRange[3]});
  }

  if (RightIsLocal) {
    MaybeGridParams[1] = Domain.MakeGridParams()
      .SetName("Right")
      .SetComm(RightData.Comm)
      .SetGlobalRange({RightData.Size})
      .SetLocalRange({&RightData.LocalRange[0], &RightData.LocalRange[3]});
  }

  Domain.CreateGrids(GridIDs, MaybeGridParams);

  constexpr int CONNECTIVITY_ID = 1;
  Domain.CreateComponent<ovk::connectivity_component>(CONNECTIVITY_ID);

  {

    auto ConnectivityComponentEditHandle = Domain.EditComponent<ovk::connectivity_component>(
      CONNECTIVITY_ID);
    ovk::connectivity_component &ConnectivityComponent = *ConnectivityComponentEditHandle;

    std::vector<ovk::elem<int,2>> ConnectivityIDs = {{1,2}, {2,1}};

    ConnectivityComponent.CreateConnectivities(ConnectivityIDs);

    if (LeftIsLocal) {

      const ovk::grid &Grid = Domain.Grid(1);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.End(0) == GlobalRange.End(0);

      auto ConnectivityMEditHandle = ConnectivityComponent.EditConnectivityM({1,2});
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

      auto ConnectivityNEditHandle = ConnectivityComponent.EditConnectivityN({2,1});
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

    if (RightIsLocal) {

      const ovk::grid &Grid = Domain.Grid(2);
      const ovk::range &GlobalRange = Grid.GlobalRange();
      const ovk::range &LocalRange = Grid.LocalRange();

      bool HasInterface = LocalRange.Begin(0) == GlobalRange.Begin(0);

      auto ConnectivityMEditHandle = ConnectivityComponent.EditConnectivityM({2,1});
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
          Destinations(0,iDonor) = RightSize[0]-1;
          Destinations(1,iDonor) = j;
          ++iDonor;
        }
      }

      auto ConnectivityNEditHandle = ConnectivityComponent.EditConnectivityN({1,2});
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
          Sources(0,iReceiver) = LeftSize[0]-2;
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

  std::vector<double> LeftDonorValues, LeftReceiverValues;
  if (LeftIsLocal) {
    ovk::range ExtendedRange(&LeftData.ExtendedRange[0], &LeftData.ExtendedRange[3]);
    const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM({1,2});
    Exchanger.CreateCollect({1,2}, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
      ExtendedRange, ovk::array_layout::ROW_MAJOR);
    Exchanger.CreateSend({1,2}, 1, ovk::data_type::DOUBLE, 1, 1);
    LeftDonorValues.resize(ConnectivityM.Size());
    const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN({2,1});
    Exchanger.CreateReceive({2,1}, 1, ovk::data_type::DOUBLE, 1, 1);
    Exchanger.CreateDisperse({2,1}, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
      ExtendedRange, ovk::array_layout::ROW_MAJOR);
    LeftReceiverValues.resize(ConnectivityN.Size());
  }

  std::vector<double> RightDonorValues, RightReceiverValues;
  if (RightIsLocal) {
    ovk::range ExtendedRange(&RightData.ExtendedRange[0], &RightData.ExtendedRange[3]);
    const ovk::connectivity_m &ConnectivityM = ConnectivityComponent.ConnectivityM({2,1});
    Exchanger.CreateCollect({2,1}, 1, ovk::collect_op::INTERPOLATE, ovk::data_type::DOUBLE, 1,
      ExtendedRange, ovk::array_layout::ROW_MAJOR);
    Exchanger.CreateSend({2,1}, 1, ovk::data_type::DOUBLE, 1, 1);
    RightDonorValues.resize(ConnectivityM.Size());
    const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN({1,2});
    Exchanger.CreateReceive({1,2}, 1, ovk::data_type::DOUBLE, 1, 1);
    Exchanger.CreateDisperse({1,2}, 1, ovk::disperse_op::OVERWRITE, ovk::data_type::DOUBLE, 1,
      ExtendedRange, ovk::array_layout::ROW_MAJOR);
    RightReceiverValues.resize(ConnectivityN.Size());
  }

  std::vector<double> LeftFieldValues;
  if (LeftIsLocal) {
    LeftFieldValues.resize(LeftData.NumExtendedPoints, -1.);
  }

  std::vector<double> RightFieldValues;
  if (RightIsLocal) {
    RightFieldValues.resize(RightData.NumExtendedPoints, 1.);
  }

  std::vector<ovk::request> Requests;

  if (LeftIsLocal) {
    double *ReceiverValues = LeftReceiverValues.data();
    ovk::request Request = Exchanger.Receive({2,1}, 1, &ReceiverValues);
    Requests.push_back(std::move(Request));
  }

  if (RightIsLocal) {
    double *ReceiverValues = RightReceiverValues.data();
    ovk::request Request = Exchanger.Receive({1,2}, 1, &ReceiverValues);
    Requests.push_back(std::move(Request));
  }

  if (LeftIsLocal) {
    const double *FieldValues = LeftFieldValues.data();
    double *DonorValues = LeftDonorValues.data();
    Exchanger.Collect({1,2}, 1, &FieldValues, &DonorValues);
    ovk::request Request = Exchanger.Send({1,2}, 1, &DonorValues);
    Requests.push_back(std::move(Request));
  }

  if (RightIsLocal) {
    const double *FieldValues = RightFieldValues.data();
    double *DonorValues = RightDonorValues.data();
    Exchanger.Collect({2,1}, 1, &FieldValues, &DonorValues);
    ovk::request Request = Exchanger.Send({2,1}, 1, &DonorValues);
    Requests.push_back(std::move(Request));
  }

  ovk::WaitAll(Requests);

  if (LeftIsLocal) {
    const double *ReceiverValues = LeftReceiverValues.data();
    double *FieldValues = LeftFieldValues.data();
    Exchanger.Disperse({2,1}, 1, &ReceiverValues, &FieldValues);
  }

  if (RightIsLocal) {
    const double *ReceiverValues = RightReceiverValues.data();
    double *FieldValues = RightFieldValues.data();
    Exchanger.Disperse({1,2}, 1, &ReceiverValues, &FieldValues);
  }

}

}
