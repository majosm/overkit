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
using examples::command_args;
using examples::command_args_parser;
#ifdef OVK_HAVE_XDMF
using examples::xdmf;
using examples::xdmf_grid_meta;
using examples::xdmf_attribute_meta;
using examples::xdmf_attribute_type;
using examples::CreateXDMF;
using examples::OpenXDMF;
#endif

namespace {
void GetCommandLineArguments(int argc, char **argv, bool &Help, int &N);
void Interface(int N);
}

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);

  int WorldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

  try {
    bool Help;
    int N;
    GetCommandLineArguments(argc, argv, Help, N);
    if (!Help) {
      Interface(N);
    }
  } catch (const std::exception &Exception) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (WorldRank == 0) {
      std::fprintf(stderr, "Encountered error:\n%s\n", Exception.what()); std::fflush(stderr);
    }
  } catch (...) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (WorldRank == 0) {
      std::fprintf(stderr, "Unknown error occurred.\n"); std::fflush(stderr);
    }
  }

  MPI_Finalize();

  return 0;

}

namespace {

void GetCommandLineArguments(int argc, char **argv, bool &Help, int &N) {

  int WorldRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &WorldRank);

  command_args_parser CommandArgsParser(WorldRank == 0);
  CommandArgsParser.SetHelpUsage("Interface [<options> ...]");
  CommandArgsParser.SetHelpDescription("Generates an overset mesh consisting of two grids "
    "overlapping along an interface.");
  CommandArgsParser.AddOption<int>("size", 'N', "Characteristic size of grids [ Default: 64 ]");

  command_args CommandArgs = CommandArgsParser.Parse({{argc}, argv});

  Help = CommandArgs.GetOptionValue<bool>("help", false);
  N = CommandArgs.GetOptionValue<int>("size", 64);

}

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

void Interface(int N) {

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

  std::array<int,3> Size = {{N,N,1}};

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
          Destinations(0,iDonor) = LeftSize[0]-1;
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

#ifdef OVK_HAVE_XDMF
  std::array<xdmf_grid_meta,2> XDMFGrids = {{
    {"Left", LeftSize},
    {"Right", RightSize}
  }};

  std::array<xdmf_attribute_meta,3> XDMFAttributes = {{
    {"State", xdmf_attribute_type::INT},
    {"BeforeExchange", xdmf_attribute_type::DOUBLE},
    {"AfterExchange", xdmf_attribute_type::DOUBLE}
  }};

  CreateXDMF("Interface.xmf", 2, MPI_COMM_WORLD, std::move(XDMFGrids), std::move(XDMFAttributes));

  auto CreateOutputState = [&Domain, CONNECTIVITY_ID](int GridID) -> ovk::field<int> {
    auto &ConnectivityComponent = Domain.Component<ovk::connectivity_component>(CONNECTIVITY_ID);
    const ovk::grid &Grid = Domain.Grid(GridID);
    ovk::field<int> OutputState(Grid.LocalRange(), 1);
    for (auto &ConnectivityID : ConnectivityComponent.LocalConnectivityNIDs()) {
      int MGridID = ConnectivityID(0);
      int NGridID = ConnectivityID(1);
      if (NGridID != GridID) continue;
      const ovk::connectivity_n &ConnectivityN = ConnectivityComponent.ConnectivityN(ConnectivityID);
      const ovk::array<int,2> &Points = ConnectivityN.Points();
      for (long long iReceiver = 0; iReceiver < ConnectivityN.Size(); ++iReceiver) {
        ovk::tuple<int> Point = {
          Points(0,iReceiver),
          Points(1,iReceiver),
          Points(2,iReceiver)
        };
        OutputState(Point) = -MGridID;
      }
    }
    return OutputState;
  };

  if (LeftIsLocal) {
    const grid_data &Data = LeftData;
    const std::array<int,6> &LocalRange = Data.LocalRange;
    xdmf XDMF = OpenXDMF("Interface.xmf", Data.Comm);
    ovk::array<ovk::field<double>> Coords({2});
    Coords(0).Resize({&LocalRange[0], &LocalRange[3]});
    Coords(1).Resize({&LocalRange[0], &LocalRange[3]});
    for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
      for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
        double U = double(i)/double(Size[0]-1);
        double V = double(j)/double(LeftSize[1]-1);
        Coords(0)(i,j,0) = 2.*(U-0.5);
        Coords(1)(i,j,0) = 2.*(V-0.5);
      }
    }
    ovk::field<int> OutputState = CreateOutputState(1);
    for (int iDim = 0; iDim < 2; ++iDim) {
      XDMF.WriteGeometry("Left", iDim, Coords(iDim));
    }
    XDMF.WriteAttribute("Left", "State", OutputState);
  }

  if (RightIsLocal) {
    const grid_data &Data = RightData;
    const std::array<int,6> &LocalRange = Data.LocalRange;
    xdmf XDMF = OpenXDMF("Interface.xmf", Data.Comm);
    ovk::array<ovk::field<double>> Coords({2});
    Coords(0).Resize({&LocalRange[0], &LocalRange[3]});
    Coords(1).Resize({&LocalRange[0], &LocalRange[3]});
    for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
      for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
        double U = double(i+LeftSize[0]-2)/double(Size[0]-1);
        double V = double(j)/double(RightSize[1]-1);
        Coords(0)(i,j,0) = 2.*(U-0.5);
        Coords(1)(i,j,0) = 2.*(V-0.5);
      }
    }
    ovk::field<int> OutputState = CreateOutputState(2);
    for (int iDim = 0; iDim < 2; ++iDim) {
      XDMF.WriteGeometry("Right", iDim, Coords(iDim));
    }
    XDMF.WriteAttribute("Right", "State", OutputState);
  }
#endif

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
    const grid_data &Data = LeftData;
    const std::array<int,6> &ExtendedRange = Data.ExtendedRange;
    LeftFieldValues.resize(Data.NumExtendedPoints);
    for (int j = ExtendedRange[1]; j < ExtendedRange[4]; ++j) {
      for (int i = ExtendedRange[0]; i < ExtendedRange[3]; ++i) {
        long long l = (ExtendedRange[4]-ExtendedRange[1])*(i-ExtendedRange[0]) +
          (j-ExtendedRange[1]);
        double U = double(i)/double(LeftSize[0]-1);
        double V = double(j)/double(LeftSize[1]-1);
        LeftFieldValues[l] = U*V;
      }
    }
  }

  std::vector<double> RightFieldValues;
  if (RightIsLocal) {
    const grid_data &Data = RightData;
    const std::array<int,6> &ExtendedRange = Data.ExtendedRange;
    RightFieldValues.resize(Data.NumExtendedPoints);
    for (int j = ExtendedRange[1]; j < ExtendedRange[4]; ++j) {
      for (int i = ExtendedRange[0]; i < ExtendedRange[3]; ++i) {
        long long l = (ExtendedRange[4]-ExtendedRange[1])*(i-ExtendedRange[0]) +
          (j-ExtendedRange[1]);
        double U = double(i)/double(RightSize[0]-1);
        double V = double(j)/double(RightSize[1]-1);
        RightFieldValues[l] = (1.-U)*(1.-V);
      }
    }
  }

#ifdef OVK_HAVE_XDMF
  if (LeftIsLocal) {
    const grid_data &Data = LeftData;
    const std::array<int,6> &LocalRange = Data.LocalRange;
    const std::array<int,6> &ExtendedRange = Data.ExtendedRange;
    xdmf XDMF = OpenXDMF("Interface.xmf", Data.Comm);
    ovk::field<double> ValuesTransposed({&LocalRange[0], &LocalRange[3]});
    for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
      for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
        long long l = (ExtendedRange[4]-ExtendedRange[1])*(i-ExtendedRange[0]) +
          (j-ExtendedRange[1]);
        ValuesTransposed(i,j,0) = LeftFieldValues[l];
      }
    }
    XDMF.WriteAttribute("Left", "BeforeExchange", ValuesTransposed);
  }

  if (RightIsLocal) {
    const grid_data &Data = RightData;
    const std::array<int,6> &LocalRange = Data.LocalRange;
    const std::array<int,6> &ExtendedRange = Data.ExtendedRange;
    xdmf XDMF = OpenXDMF("Interface.xmf", Data.Comm);
    ovk::field<double> ValuesTransposed({&LocalRange[0], &LocalRange[3]});
    for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
      for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
        long long l = (ExtendedRange[4]-ExtendedRange[1])*(i-ExtendedRange[0]) +
          (j-ExtendedRange[1]);
        ValuesTransposed(i,j,0) = RightFieldValues[l];
      }
    }
    XDMF.WriteAttribute("Right", "BeforeExchange", ValuesTransposed);
  }
#endif

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

#ifdef OVK_HAVE_XDMF
  if (LeftIsLocal) {
    const grid_data &Data = LeftData;
    const std::array<int,6> &LocalRange = Data.LocalRange;
    const std::array<int,6> &ExtendedRange = Data.ExtendedRange;
    xdmf XDMF = OpenXDMF("Interface.xmf", Data.Comm);
    ovk::field<double> ValuesTransposed({&LocalRange[0], &LocalRange[3]});
    for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
      for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
        long long l = (ExtendedRange[4]-ExtendedRange[1])*(i-ExtendedRange[0]) +
          (j-ExtendedRange[1]);
        ValuesTransposed(i,j,0) = LeftFieldValues[l];
      }
    }
    XDMF.WriteAttribute("Left", "AfterExchange", ValuesTransposed);
  }

  if (RightIsLocal) {
    const grid_data &Data = RightData;
    const std::array<int,6> &LocalRange = Data.LocalRange;
    const std::array<int,6> &ExtendedRange = Data.ExtendedRange;
    xdmf XDMF = OpenXDMF("Interface.xmf", Data.Comm);
    ovk::field<double> ValuesTransposed({&LocalRange[0], &LocalRange[3]});
    for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
      for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
        long long l = (ExtendedRange[4]-ExtendedRange[1])*(i-ExtendedRange[0]) +
          (j-ExtendedRange[1]);
        ValuesTransposed(i,j,0) = RightFieldValues[l];
      }
    }
    XDMF.WriteAttribute("Right", "AfterExchange", ValuesTransposed);
  }
#endif

}

}
