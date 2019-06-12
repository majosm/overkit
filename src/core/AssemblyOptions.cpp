// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/AssemblyOptions.hpp"

#include "ovk/core/ArrayView.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/IDMap.hpp"
#include "ovk/core/IDSet.hpp"

namespace ovk {

namespace {

template <typename T> const T &GetOption(const id_map<1,T> &Option, int GridID);
template <typename T> const T &GetOption(const id_map<2,T> &Option, int MGridID, int NGridID);
template <typename T> void SetOption(id_map<1,T> &Option, int GridID, T Value);
template <typename T> void SetOption(id_map<2,T> &Option, int MGridID, int NGridID, T Value);

}

assembly_options::assembly_options(int NumDims, const id_set<1> &GridIDs):
  NumDims_(NumDims),
  NumGrids_(GridIDs.Count())
{

  OVK_DEBUG_ASSERT(NumDims_ == 2 || NumDims_ == 3, "Invalid dimension.");
  OVK_DEBUG_ASSERT(NumGrids_ >= 0, "Invalid grid count.");

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      elem<int,2> IDPair = {MGridID,NGridID};
      Overlappable_.Insert(IDPair, false);
      OverlapTolerance_.Insert(IDPair, 1.e-12);
    }
    OverlapAccelDepthAdjust_.Insert(MGridID, 0.);
    OverlapAccelResolutionAdjust_.Insert(MGridID, 0.);
  }

  for (int MGridID : GridIDs) {
    InferBoundaries_.Insert(MGridID, false);
    for (int NGridID : GridIDs) {
      elem<int,2> IDPair = {MGridID,NGridID};
      CutBoundaryHoles_.Insert(IDPair, false);
    }
  }

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      elem<int,2> IDPair = {MGridID,NGridID};
      Occludes_.Insert(IDPair, occludes::NONE);
      EdgePadding_.Insert(IDPair, 0);
    }
    EdgeSmoothing_.Insert(MGridID, 0);
  }

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      elem<int,2> IDPair = {MGridID,NGridID};
      ConnectionType_.Insert(IDPair, connection_type::NONE);
      MinimizeOverlap_.Insert(IDPair, false);
    }
    FringeSize_.Insert(MGridID, 0);
  }

}

bool assembly_options::Overlappable(int MGridID, int NGridID) const {

  return GetOption(Overlappable_, MGridID, NGridID);

}

assembly_options &assembly_options::SetOverlappable(int MGridID, int NGridID, bool Overlappable) {

  SetOption(Overlappable_, MGridID, NGridID, Overlappable);

  return *this;

}

double assembly_options::OverlapTolerance(int MGridID, int NGridID) const {

  return GetOption(OverlapTolerance_, MGridID, NGridID);

}

assembly_options &assembly_options::SetOverlapTolerance(int MGridID, int NGridID, double
  OverlapTolerance) {

  OVK_DEBUG_ASSERT(OverlapTolerance >= 0., "Invalid overlap tolerance value.");

  SetOption(OverlapTolerance_, MGridID, NGridID, OverlapTolerance);

  return *this;

}

double assembly_options::OverlapAccelDepthAdjust(int MGridID) const {

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");

  return GetOption(OverlapAccelDepthAdjust_, MGridID);

}

assembly_options &assembly_options::SetOverlapAccelDepthAdjust(int MGridID, double
  OverlapAccelDepthAdjust) {

  SetOption(OverlapAccelDepthAdjust_, MGridID, OverlapAccelDepthAdjust);

  return *this;

}

double assembly_options::OverlapAccelResolutionAdjust(int MGridID) const {

  return GetOption(OverlapAccelResolutionAdjust_, MGridID);

}

assembly_options &assembly_options::SetOverlapAccelResolutionAdjust(int MGridID, double
  OverlapAccelResolutionAdjust) {

  SetOption(OverlapAccelResolutionAdjust_, MGridID, OverlapAccelResolutionAdjust);

  return *this;

}

bool assembly_options::InferBoundaries(int GridID) const {

  return GetOption(InferBoundaries_, GridID);

}

assembly_options &assembly_options::SetInferBoundaries(int GridID, bool InferBoundaries) {

  SetOption(InferBoundaries_, GridID, InferBoundaries);

  return *this;

}

bool assembly_options::CutBoundaryHoles(int MGridID, int NGridID) const {

  return GetOption(CutBoundaryHoles_, MGridID, NGridID);

}

assembly_options &assembly_options::SetCutBoundaryHoles(int MGridID, int NGridID, bool
  CutBoundaryHoles) {

  SetOption(CutBoundaryHoles_, MGridID, NGridID, CutBoundaryHoles);

  return *this;

}

occludes assembly_options::Occludes(int MGridID, int NGridID) const {

  return GetOption(Occludes_, MGridID, NGridID);

}

assembly_options &assembly_options::SetOccludes(int MGridID, int NGridID, occludes Occludes) {

  OVK_DEBUG_ASSERT(ValidOccludes(Occludes), "Invalid occludes value.");

  SetOption(Occludes_, MGridID, NGridID, Occludes);

  return *this;

}

int assembly_options::EdgePadding(int MGridID, int NGridID) const {

  return GetOption(EdgePadding_, MGridID, NGridID);

}

assembly_options &assembly_options::SetEdgePadding(int MGridID, int NGridID, int EdgePadding) {

  OVK_DEBUG_ASSERT(EdgePadding >= 0, "Invalid edge padding value.");

  SetOption(EdgePadding_, MGridID, NGridID, EdgePadding);

  return *this;

}

int assembly_options::EdgeSmoothing(int NGridID) const {

  return GetOption(EdgeSmoothing_, NGridID);

}

assembly_options &assembly_options::SetEdgeSmoothing(int NGridID, int EdgeSmoothing) {

  OVK_DEBUG_ASSERT(EdgeSmoothing >= 0, "Invalid edge smoothing value.");

  SetOption(EdgeSmoothing_, NGridID, EdgeSmoothing);

  return *this;

}

connection_type assembly_options::ConnectionType(int MGridID, int NGridID) const {

  return GetOption(ConnectionType_, MGridID, NGridID);

}

assembly_options &assembly_options::SetConnectionType(int MGridID, int NGridID, connection_type
  ConnectionType) {

  OVK_DEBUG_ASSERT(ValidConnectionType(ConnectionType), "Invalid connection type.");

  SetOption(ConnectionType_, MGridID, NGridID, ConnectionType);

  return *this;

}

int assembly_options::FringeSize(int NGridID) const {

  return GetOption(FringeSize_, NGridID);

}

assembly_options &assembly_options::SetFringeSize(int NGridID, int FringeSize) {

  OVK_DEBUG_ASSERT(FringeSize >= 0, "Invalid fringe size.");

  SetOption(FringeSize_, NGridID, FringeSize);

  return *this;

}

bool assembly_options::MinimizeOverlap(int MGridID, int NGridID) const {

  return GetOption(MinimizeOverlap_, MGridID, NGridID);

}

assembly_options &assembly_options::SetMinimizeOverlap(int MGridID, int NGridID, bool
  MinimizeOverlap) {

  SetOption(MinimizeOverlap_, MGridID, NGridID, MinimizeOverlap);

  return *this;

}

void assembly_options::PrintOptions_() {

  id_set<1> GridIDs = OverlapAccelDepthAdjust_.Keys();

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      std::printf("Overlappable(%i,%i) = %c\n", MGridID, NGridID, Overlappable_(MGridID,NGridID) ?
        'T' : 'F');
    }
  }

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      std::printf("OverlapTolerance(%i,%i) = %16.8f\n", MGridID, NGridID, OverlapTolerance_(MGridID,
        NGridID));
    }
  }

  for (int MGridID : GridIDs) {
    std::printf("OverlapAccelDepthAdjust(%i) = %16.8f\n", MGridID, OverlapAccelDepthAdjust_(
      MGridID));
  }

  for (int MGridID : GridIDs) {
    std::printf("OverlapAccelResolutionAdjust(%i) = %16.8f\n", MGridID,
      OverlapAccelResolutionAdjust_(MGridID));
  }

  for (int MGridID : GridIDs) {
    std::printf("InferBoundaries(%i) = %c\n", MGridID, InferBoundaries_(MGridID) ? 'T' : 'F');
  }

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      std::printf("CutBoundaryHoles(%i,%i) = %c\n", MGridID, NGridID, CutBoundaryHoles_(MGridID,
        NGridID) ? 'T' : 'F');
    }
  }

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      std::string OccludesString;
      switch (Occludes_(MGridID,NGridID)) {
      case occludes::NONE:
        OccludesString = "NONE";
        break;
      case occludes::ALL:
        OccludesString = "ALL";
        break;
      case occludes::COARSE:
        OccludesString = "COARSE";
        break;
      }
      std::printf("Occludes(%i,%i) = %s\n", MGridID, NGridID, OccludesString.c_str());
    }
  }

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      std::printf("EdgePadding(%i,%i) = %i\n", MGridID, NGridID, EdgePadding_(MGridID,NGridID));
    }
  }

  for (int MGridID : GridIDs) {
    std::printf("EdgeSmoothing(%i) = %i\n", MGridID, EdgeSmoothing_(MGridID));
  }

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      std::string ConnectionTypeString;
      switch (ConnectionType_(MGridID,NGridID)) {
      case connection_type::NONE:
        ConnectionTypeString = "NONE";
        break;
      case connection_type::NEAREST:
        ConnectionTypeString = "NEAREST";
        break;
      case connection_type::LINEAR:
        ConnectionTypeString = "LINEAR";
        break;
      case connection_type::CUBIC:
        ConnectionTypeString = "CUBIC";
        break;
      }
      std::printf("ConnectionType(%i,%i) = %s\n", MGridID, NGridID, ConnectionTypeString.c_str());
    }
  }

  for (int MGridID : GridIDs) {
    for (int NGridID : GridIDs) {
      std::printf("MinimizeOverlap(%i,%i) = %c\n", MGridID, NGridID, MinimizeOverlap_(MGridID,
        NGridID) ? 'T' : 'F');
    }
  }

  for (int MGridID : GridIDs) {
    std::printf("FringeSize(%i) = %i\n", MGridID, FringeSize_(MGridID));
  }

}

namespace {

template <typename T> const T &GetOption(const id_map<1,T> &Option, int GridID) {

  OVK_DEBUG_ASSERT(GridID >= 0, "Invalid grid ID.");
  OVK_DEBUG_ASSERT(Option.Contains(GridID), "No option for grid %i.", GridID);

  return Option(GridID);

}

template <typename T> const T &GetOption(const id_map<2,T> &Option, int MGridID, int NGridID)
  {

  OVK_DEBUG_ASSERT(MGridID >= 0, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0, "Invalid N grid ID.");
  OVK_DEBUG_ASSERT(Option.Contains({MGridID,NGridID}), "No option for grid pair (%i,%i).", MGridID,
    NGridID);

  return Option(MGridID,NGridID);

}

template <typename T> void SetOption(id_map<1,T> &Option, int GridID, T Value) {

  OVK_DEBUG_ASSERT(GridID >= 0 || GridID == ALL_GRIDS, "Invalid grid ID.");

  if (GridID == ALL_GRIDS) {
    for (auto &Entry : Option) {
      Entry.Value() = Value;
    }
  } else {
    OVK_DEBUG_ASSERT(Option.Contains(GridID), "No option for grid %i.", GridID);
    Option(GridID) = Value;
  }

}

template <typename T> void SetOption(id_map<2,T> &Option, int MGridID, int NGridID, T Value) {

  OVK_DEBUG_ASSERT(MGridID >= 0 || MGridID == ALL_GRIDS, "Invalid M grid ID.");
  OVK_DEBUG_ASSERT(NGridID >= 0 || NGridID == ALL_GRIDS, "Invalid N grid ID.");

  if (MGridID == ALL_GRIDS && NGridID == ALL_GRIDS) {
    for (auto &Entry : Option) {
      Entry.Value() = Value;
    }
  } else if (MGridID == ALL_GRIDS) {
    if (OVK_DEBUG) {
      bool Found = false;
      for (auto &Entry : Option) {
        if (Entry.Key(1) == NGridID) Found = true;
      }
      OVK_DEBUG_ASSERT(Found, "No option for grid %i.", NGridID);
    }
    for (auto &Entry : Option) {
      if (Entry.Key(1) == NGridID) {
        Entry.Value() = Value;
      }
    }
  } else if (NGridID == ALL_GRIDS) {
    if (OVK_DEBUG) {
      bool Found = false;
      for (auto &Entry : Option) {
        if (Entry.Key(0) == MGridID) Found = true;
      }
      OVK_DEBUG_ASSERT(Found, "No option for grid %i.", MGridID);
    }
    for (auto &Entry : Option) {
      if (Entry.Key(0) == MGridID) {
        Entry.Value() = Value;
      }
    }
  } else {
    OVK_DEBUG_ASSERT(Option.Contains({MGridID,NGridID}), "No option for grid pair (%i,%i).",
      MGridID, NGridID);
    Option(MGridID,NGridID) = Value;
  }

}

} 

}
