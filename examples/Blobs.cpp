// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include <overkit.hpp>

#include "examples/Common.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <memory>
#include <utility>
#include <vector>

using examples::DecomposeDomain;
using examples::CreateCartesianDecompDims;
using examples::CartesianDecomp;
using examples::PI;
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
void Blobs(int N);
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
      Blobs(N);
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
  CommandArgsParser.SetHelpUsage("Blobs [<options> ...]");
  CommandArgsParser.SetHelpDescription("Generates an overset mesh representing a rectangular "
    "domain with several blob-like holes in it.");
  CommandArgsParser.AddOption<int>("size", 'N', "Characteristic size of grids [ Default: 81 ]");

  command_args CommandArgs = CommandArgsParser.Parse({{argc}, argv});

  Help = CommandArgs.GetOptionValue<bool>("help", false);
  N = CommandArgs.GetOptionValue<int>("size", 81);

}

struct grid_data {
  MPI_Comm Comm = MPI_COMM_NULL;
  std::array<int,3> Size = {{0,0,1}};
  std::array<bool,3> Periodic = {{false,false,false}};
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

void Blobs(int N) {

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

  std::array<int,3> BackgroundSize = {{N,N,1}};
  std::array<int,3> BlobSize = {{N/2,N,1}};

  long long NumBackgroundPoints = BackgroundSize[0]*BackgroundSize[1];
  long long NumBlobPoints = BlobSize[0]*BlobSize[1];

  std::array<long long,4> NumPointsPerGrid = {{
    NumBackgroundPoints,
    NumBlobPoints,
    NumBlobPoints,
    NumBlobPoints
  }};

  std::array<int,8> GridProcRanges;

  DecomposeDomain(NumPointsPerGrid, NumWorldProcs, GridProcRanges);

  bool BackgroundIsLocal = WorldRank >= GridProcRanges[0] && WorldRank < GridProcRanges[1];
  bool Blob1IsLocal = WorldRank >= GridProcRanges[2] && WorldRank < GridProcRanges[3];
  bool Blob2IsLocal = WorldRank >= GridProcRanges[4] && WorldRank < GridProcRanges[5];
  bool Blob3IsLocal = WorldRank >= GridProcRanges[6] && WorldRank < GridProcRanges[7];

  MPI_Comm BackgroundComm, Blob1Comm, Blob2Comm, Blob3Comm;

  MPI_Comm_split(MPI_COMM_WORLD, BackgroundIsLocal ? 0 : MPI_UNDEFINED, WorldRank, &BackgroundComm);
  MPI_Comm_split(MPI_COMM_WORLD, Blob1IsLocal ? 0 : MPI_UNDEFINED, WorldRank, &Blob1Comm);
  MPI_Comm_split(MPI_COMM_WORLD, Blob2IsLocal ? 0 : MPI_UNDEFINED, WorldRank, &Blob2Comm);
  MPI_Comm_split(MPI_COMM_WORLD, Blob3IsLocal ? 0 : MPI_UNDEFINED, WorldRank, &Blob3Comm);

  grid_data BackgroundData, Blob1Data, Blob2Data, Blob3Data;

  if (BackgroundIsLocal) {
    grid_data &Data = BackgroundData;
    int NumGridProcs;
    MPI_Comm_size(BackgroundComm, &NumGridProcs);
    std::array<int,3> CartDims = CreateCartesianDecompDims(NumGridProcs, 2, {{0,0,1}});
    std::array<int,3> CartPeriods = {{0,0,0}};
    MPI_Cart_create(BackgroundComm, 2, CartDims.data(), CartPeriods.data(), 1, &Data.Comm);
    MPI_Comm_free(&BackgroundComm);
    Data.Size = BackgroundSize;
    Data.Periodic = {{false,false,false}};
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
  }

  if (Blob1IsLocal) {
    grid_data &Data = Blob1Data;
    int NumGridProcs;
    MPI_Comm_size(Blob1Comm, &NumGridProcs);
    std::array<int,3> CartDims = CreateCartesianDecompDims(NumGridProcs, 2, {{0,0,1}});
    if (CartDims[1] < CartDims[0]) std::swap(CartDims[0], CartDims[1]);
    std::array<int,3> CartPeriods = {{0,1,0}};
    MPI_Cart_create(Blob1Comm, 2, CartDims.data(), CartPeriods.data(), 1, &Data.Comm);
    MPI_Comm_free(&Blob1Comm);
    Data.Size = BlobSize;
    Data.Periodic = {{false,true,false}};
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
  }

  if (Blob2IsLocal) {
    grid_data &Data = Blob2Data;
    int NumGridProcs;
    MPI_Comm_size(Blob2Comm, &NumGridProcs);
    std::array<int,3> CartDims = CreateCartesianDecompDims(NumGridProcs, 2, {{0,0,1}});
    if (CartDims[1] < CartDims[0]) std::swap(CartDims[0], CartDims[1]);
    std::array<int,3> CartPeriods = {{0,1,0}};
    MPI_Cart_create(Blob2Comm, 2, CartDims.data(), CartPeriods.data(), 1, &Data.Comm);
    MPI_Comm_free(&Blob2Comm);
    Data.Size = BlobSize;
    Data.Periodic = {{false,true,false}};
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
  }

  if (Blob3IsLocal) {
    grid_data &Data = Blob3Data;
    int NumGridProcs;
    MPI_Comm_size(Blob3Comm, &NumGridProcs);
    std::array<int,3> CartDims = CreateCartesianDecompDims(NumGridProcs, 2, {{0,0,1}});
    if (CartDims[1] < CartDims[0]) std::swap(CartDims[0], CartDims[1]);
    std::array<int,3> CartPeriods = {{0,1,0}};
    MPI_Cart_create(Blob3Comm, 2, CartDims.data(), CartPeriods.data(), 1, &Data.Comm);
    MPI_Comm_free(&Blob3Comm);
    Data.Size = BlobSize;
    Data.Periodic = {{false,true,false}};
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
  }

  constexpr int BACKGROUND_ID = 1;
  constexpr int BLOB_1_ID = 2;
  constexpr int BLOB_2_ID = 3;
  constexpr int BLOB_3_ID = 4;

  std::array<int,4> GridIDs = {{
    BACKGROUND_ID,
    BLOB_1_ID,
    BLOB_2_ID,
    BLOB_3_ID
  }};

  std::array<ovk::optional<ovk::grid::params>,4> MaybeGridParams;

  if (BackgroundIsLocal) {
    const grid_data &Data = BackgroundData;
    MaybeGridParams[0] = ovk::grid::params()
      .SetName("Background")
      .SetDimension(2)
      .SetComm(Data.Comm)
      .SetGlobalRange({Data.Size})
      .SetLocalRange({&Data.LocalRange[0], &Data.LocalRange[3]});
  }

  if (Blob1IsLocal) {
    const grid_data &Data = Blob1Data;
    MaybeGridParams[1] = ovk::grid::params()
      .SetName("Blob1")
      .SetDimension(2)
      .SetComm(Data.Comm)
      .SetGlobalRange({Data.Size})
      .SetLocalRange({&Data.LocalRange[0], &Data.LocalRange[3]})
      .SetPeriodic(Data.Periodic)
      .SetPeriodicStorage(ovk::periodic_storage::UNIQUE);
  }

  if (Blob2IsLocal) {
    const grid_data &Data = Blob2Data;
    MaybeGridParams[2] = ovk::grid::params()
      .SetName("Blob2")
      .SetDimension(2)
      .SetComm(Data.Comm)
      .SetGlobalRange({Data.Size})
      .SetLocalRange({&Data.LocalRange[0], &Data.LocalRange[3]})
      .SetPeriodic(Data.Periodic)
      .SetPeriodicStorage(ovk::periodic_storage::UNIQUE);
  }

  if (Blob3IsLocal) {
    const grid_data &Data = Blob3Data;
    MaybeGridParams[3] = ovk::grid::params()
      .SetName("Blob3")
      .SetDimension(2)
      .SetComm(Data.Comm)
      .SetGlobalRange({Data.Size})
      .SetLocalRange({&Data.LocalRange[0], &Data.LocalRange[3]})
      .SetPeriodic(Data.Periodic)
      .SetPeriodicStorage(ovk::periodic_storage::UNIQUE);
  }

  Domain.CreateGrids(GridIDs, MaybeGridParams);

  constexpr int GEOMETRY_ID = 1;
  constexpr int STATE_ID = 2;
  constexpr int OVERLAP_ID = 3;
  constexpr int CONNECTIVITY_ID = 4;

  Domain.CreateComponent<ovk::geometry_component>(GEOMETRY_ID);
  Domain.CreateComponent<ovk::state_component>(STATE_ID);
  Domain.CreateComponent<ovk::overlap_component>(OVERLAP_ID);
  Domain.CreateComponent<ovk::connectivity_component>(CONNECTIVITY_ID);

  {

    auto GeometryComponentHandle = Domain.EditComponent<ovk::geometry_component>(GEOMETRY_ID);
    ovk::geometry_component &GeometryComponent = *GeometryComponentHandle;

    auto StateComponentHandle = Domain.EditComponent<ovk::state_component>(STATE_ID);
    ovk::state_component &StateComponent = *StateComponentHandle;

    std::array<ovk::optional<ovk::geometry::params>,4> MaybeGeometryParams;
    if (BackgroundIsLocal) {
      MaybeGeometryParams[0] = ovk::geometry::params()
        .SetType(ovk::geometry_type::UNIFORM);
    }

    GeometryComponent.CreateGeometries(GridIDs, MaybeGeometryParams);

    StateComponent.CreateStates(GridIDs);

    if (BackgroundIsLocal) {
      const grid_data &Data = BackgroundData;
      const std::array<int,6> &LocalRange = Data.LocalRange;
      auto GeometryHandle = GeometryComponent.EditGeometry(BACKGROUND_ID);
      ovk::geometry &Geometry = *GeometryHandle;
      auto CoordsHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsHandle;
      for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
        for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
          double U = double(i)/double(BackgroundSize[0]-1);
          double V = double(j)/double(BackgroundSize[1]-1);
          Coords(0)(i,j,0) = 2.*(U-0.5);
          Coords(1)(i,j,0) = 2.*(V-0.5);
        }
      }
    }

    constexpr double SeparationScale = 0.8;

    constexpr double Shift = 4.;

    if (Blob1IsLocal) {
      const grid_data &Data = Blob1Data;
      const std::array<int,6> &LocalRange = Data.LocalRange;
      auto GeometryHandle = GeometryComponent.EditGeometry(BLOB_1_ID);
      ovk::geometry &Geometry = *GeometryHandle;
      auto CoordsHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsHandle;
      for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
        for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
          double U = double(i)/double(BlobSize[0]-1);
          double V = double(j)/double(BlobSize[1]);
          double Theta = 2.*PI*V;
          double RMin = 0.1*(1.+0.2*std::sin(3.*Theta) + 0.1*std::sin(2.*Theta+PI/4.));
          double RMax = 0.5 + RMin;
          double RShift = (1. - std::pow(2., Shift))*RMin;
          double RMinRel = RMin - RShift;
          double RMaxRel = RMax - RShift;
          double Radius = RMinRel * std::pow(RMaxRel/RMinRel, U) + RShift;
          Coords(0)(i,j,0) = -0.425*SeparationScale + Radius*std::cos(Theta);
          Coords(1)(i,j,0) = -0.025*SeparationScale + Radius*std::sin(Theta);
        }
      }
      auto StateHandle = StateComponent.EditState(BLOB_1_ID);
      auto FlagsHandle = StateHandle->EditFlags();
      ovk::distributed_field<ovk::state_flags> &Flags = *FlagsHandle;
      if (LocalRange[0] == 0) {
        for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
          Flags(0,j,0) |= ovk::state_flags::DOMAIN_BOUNDARY;
        }
      }
    }

    if (Blob2IsLocal) {
      const grid_data &Data = Blob2Data;
      const std::array<int,6> &LocalRange = Data.LocalRange;
      auto GeometryHandle = GeometryComponent.EditGeometry(BLOB_2_ID);
      ovk::geometry &Geometry = *GeometryHandle;
      auto CoordsHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsHandle;
      for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
        for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
          double U = double(i)/double(BlobSize[0]-1);
          double V = double(j)/double(BlobSize[1]);
          double Theta = 2.*PI*V;
          double RMin = 0.1*(1.+0.2*std::sin(4.*Theta+PI/4.) + 0.1*std::sin(2.*Theta));
          double RMax = 0.5 + RMin;
          double RShift = (1. - std::pow(2., Shift))*RMin;
          double RMinRel = RMin - RShift;
          double RMaxRel = RMax - RShift;
          double Radius = RMinRel * std::pow(RMaxRel/RMinRel, U) + RShift;
          Coords(0)(i,j,0) = 0.075*SeparationScale + Radius*std::cos(Theta);
          Coords(1)(i,j,0) = 0.425*SeparationScale + Radius*std::sin(Theta);
        }
      }
      auto StateHandle = StateComponent.EditState(BLOB_2_ID);
      auto FlagsHandle = StateHandle->EditFlags();
      ovk::distributed_field<ovk::state_flags> &Flags = *FlagsHandle;
      if (LocalRange[0] == 0) {
        for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
          Flags(0,j,0) |= ovk::state_flags::DOMAIN_BOUNDARY;
        }
      }
    }

    if (Blob3IsLocal) {
      const grid_data &Data = Blob3Data;
      const std::array<int,6> &LocalRange = Data.LocalRange;
      auto GeometryHandle = GeometryComponent.EditGeometry(BLOB_3_ID);
      ovk::geometry &Geometry = *GeometryHandle;
      auto CoordsHandle = Geometry.EditCoords();
      ovk::array<ovk::distributed_field<double>> &Coords = *CoordsHandle;
      for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
        for (int i = LocalRange[0]; i < LocalRange[3]; ++i) {
          double U = double(i)/double(BlobSize[0]-1);
          double V = double(j)/double(BlobSize[1]);
          double Theta = 2.*PI*V;
          double RMin = 0.1*(1.+0.2*std::sin(5.*Theta+PI/4.) + 0.1*std::sin(3.*Theta));
          double RMax = 0.5 + RMin;
          double RShift = (1. - std::pow(2., Shift))*RMin;
          double RMinRel = RMin - RShift;
          double RMaxRel = RMax - RShift;
          double Radius = RMinRel * std::pow(RMaxRel/RMinRel, U) + RShift;
          Coords(0)(i,j,0) = 0.375*SeparationScale + Radius*std::cos(Theta);
          Coords(1)(i,j,0) = -0.375*SeparationScale + Radius*std::sin(Theta);
        }
      }
      auto StateHandle = StateComponent.EditState(BLOB_3_ID);
      auto FlagsHandle = StateHandle->EditFlags();
      ovk::distributed_field<ovk::state_flags> &Flags = *FlagsHandle;
      if (LocalRange[0] == 0) {
        for (int j = LocalRange[1]; j < LocalRange[4]; ++j) {
          Flags(0,j,0) |= ovk::state_flags::DOMAIN_BOUNDARY;
        }
      }
    }

  }

  ovk::assembler Assembler = ovk::CreateAssembler(Context);

  Assembler.Bind(Domain, ovk::assembler::bindings()
    .SetGeometryComponentID(GEOMETRY_ID)
    .SetStateComponentID(STATE_ID)
    .SetOverlapComponentID(OVERLAP_ID)
    .SetConnectivityComponentID(CONNECTIVITY_ID)
  );

  {
    auto OptionsHandle = Assembler.EditOptions();
    ovk::assembler::options &Options = *OptionsHandle;
    Options.SetOverlappable({ovk::ALL_GRIDS,ovk::ALL_GRIDS}, true);
    Options.SetInferBoundaries(ovk::ALL_GRIDS, true);
    Options.SetCutBoundaryHoles({ovk::ALL_GRIDS,ovk::ALL_GRIDS}, true);
    Options.SetOccludes({ovk::ALL_GRIDS,ovk::ALL_GRIDS}, ovk::occludes::COARSE);
    Options.SetOccludes({ovk::ALL_GRIDS,1}, ovk::occludes::ALL);
    Options.SetEdgePadding({ovk::ALL_GRIDS,ovk::ALL_GRIDS}, 3);
    Options.SetEdgeSmoothing(ovk::ALL_GRIDS, 2);
    Options.SetConnectionType({ovk::ALL_GRIDS,ovk::ALL_GRIDS}, ovk::connection_type::LINEAR);
    Options.SetFringeSize(ovk::ALL_GRIDS, 2);
    Options.SetMinimizeOverlap({ovk::ALL_GRIDS,ovk::ALL_GRIDS}, true);
  }

  Assembler.Assemble();

#ifdef OVK_HAVE_XDMF
  std::array<xdmf_grid_meta,4> XDMFGrids = {{
    {"Background", BackgroundSize},
    {"Blob1", BlobSize},
    {"Blob2", BlobSize},
    {"Blob3", BlobSize},
  }};

  std::array<xdmf_attribute_meta,1> XDMFAttributes = {{
    {"State", xdmf_attribute_type::INT}
  }};

  CreateXDMF("Blobs.xmf", 2, MPI_COMM_WORLD, std::move(XDMFGrids), std::move(XDMFAttributes));

  auto CreateOutputState = [&Domain, STATE_ID, CONNECTIVITY_ID](int GridID) -> ovk::field<int> {
    auto &StateComponent = Domain.Component<ovk::state_component>(STATE_ID);
    auto &ConnectivityComponent = Domain.Component<ovk::connectivity_component>(CONNECTIVITY_ID);
    const ovk::grid &Grid = Domain.Grid(GridID);
    const ovk::state &State = StateComponent.State(GridID);
    const ovk::distributed_field<ovk::state_flags> &Flags = State.Flags();
    const ovk::range &LocalRange = Grid.LocalRange();
    ovk::field<int> OutputState(LocalRange, 1);
    for (int k = LocalRange.Begin(2); k < LocalRange.End(2); ++k) {
      for (int j = LocalRange.Begin(1); j < LocalRange.End(1); ++j) {
        for (int i = LocalRange.Begin(0); i < LocalRange.End(0); ++i) {
          ovk::tuple<int> Point = {i,j,k};
          if ((Flags(Point) & ovk::state_flags::ACTIVE) == ovk::state_flags::NONE) {
            OutputState(Point) = 0;
          }
        }
      }
    }
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

  auto &GeometryComponent = Domain.Component<ovk::geometry_component>(GEOMETRY_ID);

  if (BackgroundIsLocal) {
    xdmf XDMF = OpenXDMF("Blobs.xmf", BackgroundData.Comm);
    const ovk::geometry &Geometry = GeometryComponent.Geometry(BACKGROUND_ID);
    const ovk::array<ovk::distributed_field<double>> &Coords = Geometry.Coords();
    ovk::field<int> OutputState = CreateOutputState(BACKGROUND_ID);
    for (int iDim = 0; iDim < 2; ++iDim) {
      XDMF.WriteGeometry("Background", iDim, Coords(iDim), Coords(iDim).LocalRange());
    }
    XDMF.WriteAttribute("Background", "State", OutputState);
  }

  if (Blob1IsLocal) {
    xdmf XDMF = OpenXDMF("Blobs.xmf", Blob1Data.Comm);
    const ovk::geometry &Geometry = GeometryComponent.Geometry(BLOB_1_ID);
    const ovk::array<ovk::distributed_field<double>> &Coords = Geometry.Coords();
    ovk::field<int> OutputState = CreateOutputState(BLOB_1_ID);
    for (int iDim = 0; iDim < 2; ++iDim) {
      XDMF.WriteGeometry("Blob1", iDim, Coords(iDim), Coords(iDim).LocalRange());
    }
    XDMF.WriteAttribute("Blob1", "State", OutputState);
  }

  if (Blob2IsLocal) {
    xdmf XDMF = OpenXDMF("Blobs.xmf", Blob2Data.Comm);
    const ovk::geometry &Geometry = GeometryComponent.Geometry(BLOB_2_ID);
    const ovk::array<ovk::distributed_field<double>> &Coords = Geometry.Coords();
    ovk::field<int> OutputState = CreateOutputState(BLOB_2_ID);
    for (int iDim = 0; iDim < 2; ++iDim) {
      XDMF.WriteGeometry("Blob2", iDim, Coords(iDim), Coords(iDim).LocalRange());
    }
    XDMF.WriteAttribute("Blob2", "State", OutputState);
  }

  if (Blob3IsLocal) {
    xdmf XDMF = OpenXDMF("Blobs.xmf", Blob3Data.Comm);
    const ovk::geometry &Geometry = GeometryComponent.Geometry(BLOB_3_ID);
    const ovk::array<ovk::distributed_field<double>> &Coords = Geometry.Coords();
    ovk::field<int> OutputState = CreateOutputState(BLOB_3_ID);
    for (int iDim = 0; iDim < 2; ++iDim) {
      XDMF.WriteGeometry("Blob3", iDim, Coords(iDim), Coords(iDim).LocalRange());
    }
    XDMF.WriteAttribute("Blob3", "State", OutputState);
  }
#endif

}

}
