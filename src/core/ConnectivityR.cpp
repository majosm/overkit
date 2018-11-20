// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ConnectivityR.hpp"

#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/ErrorHandler.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Logger.hpp"

#include <mpi.h>

#include <vector>

namespace ovk {

namespace {

bool EditingPoints(const connectivity_r &Receivers);
bool EditingSources(const connectivity_r &Receivers);
bool EditingSourceRanks(const connectivity_r &Receivers);

void DefaultEdits(connectivity_r::edits &Edits);

}

namespace core {

void CreateConnectivityReceiverSide(connectivity_r &Receivers, const grid &Grid, int SourceGridID,
  core::logger &Logger, core::error_handler &ErrorHandler) {

  int GridID;
  int NumDims;
  MPI_Comm Comm;
  int CommSize, CommRank;
  GetGridID(Grid, GridID);
  GetGridDimension(Grid, NumDims);
  GetGridComm(Grid, Comm);
  GetGridCommSize(Grid, CommSize);
  GetGridCommRank(Grid, CommRank);

  MPI_Barrier(Comm);

  Receivers.Logger_ = &Logger;
  Receivers.ErrorHandler_ = &ErrorHandler;

  Receivers.GridID_ = GridID;
  Receivers.SourceGridID_ = SourceGridID;
  Receivers.NumDims_ = NumDims;
  Receivers.Comm_ = Comm;
  Receivers.CommSize_ = CommSize;
  Receivers.CommRank_ = CommRank;
  Receivers.Count_ = 0;

  Receivers.Grid_ = &Grid;

  DefaultEdits(Receivers.Edits_);

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers.Points_[iDim] = nullptr;
  }
  Receivers.PointsEditRefCount_ = 0;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers.Sources_[iDim] = nullptr;
  }
  Receivers.SourcesEditRefCount_ = 0;

  Receivers.SourceRanksEditRefCount_ = 0;

  MPI_Barrier(Receivers.Comm_);

}

void DestroyConnectivityReceiverSide(connectivity_r &Receivers) {

  MPI_Barrier(Receivers.Comm_);

  Receivers.PointsFlat_.clear();
  Receivers.SourcesFlat_.clear();
  Receivers.SourceRanks_.clear();

  MPI_Barrier(Receivers.Comm_);

}

}

void GetConnectivityReceiverSideGridID(const connectivity_r &Receivers, int &GridID) {

  GridID = Receivers.GridID_;

}

void GetConnectivityReceiverSideSourceGridID(const connectivity_r &Receivers, int &SourceGridID) {

  SourceGridID = Receivers.SourceGridID_;

}

void GetConnectivityReceiverSideDimension(const connectivity_r &Receivers, int &NumDims) {

  NumDims = Receivers.NumDims_;

}

void GetConnectivityReceiverSideComm(const connectivity_r &Receivers, MPI_Comm &Comm) {

  Comm = Receivers.Comm_;

}

void GetConnectivityReceiverSideCommSize(const connectivity_r &Receivers, int &CommSize) {

  CommSize = Receivers.CommSize_;

}

void GetConnectivityReceiverSideCommRank(const connectivity_r &Receivers, int &CommRank) {

  CommRank = Receivers.CommRank_;

}

void GetConnectivityReceiverSideCount(const connectivity_r &Receivers, long long &NumReceivers) {

  NumReceivers = Receivers.Count_;

}

void GetConnectivityReceiverSideGrid(const connectivity_r &Receivers, const grid *&Grid) {

  Grid = Receivers.Grid_;

}

void ResizeReceivers(connectivity_r &Receivers, long long NumReceivers) {

  OVK_DEBUG_ASSERT(NumReceivers >= 0, "Invalid receiver count.");

  MPI_Barrier(Receivers.Comm_);

  OVK_DEBUG_ASSERT(!EditingPoints(Receivers), "Cannot resize receivers while editing points.");
  OVK_DEBUG_ASSERT(!EditingSources(Receivers), "Cannot resize receivers while editing sources.");
  OVK_DEBUG_ASSERT(!EditingSourceRanks(Receivers), "Cannot resize receivers while editing source "
    "ranks.");

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers.Points_[iDim] = nullptr;
  }
  Receivers.PointsFlat_.clear();
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Receivers.Sources_[iDim] = nullptr;
  }
  Receivers.SourcesFlat_.clear();
  Receivers.SourceRanks_.clear();

  Receivers.Count_ = NumReceivers;

  if (NumReceivers > 0) {

    Receivers.PointsFlat_.resize(MAX_DIMS*NumReceivers);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Receivers.Points_[iDim] = Receivers.PointsFlat_.data() + NumReceivers*iDim;
    }
    Receivers.SourcesFlat_.resize(MAX_DIMS*NumReceivers);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Receivers.Sources_[iDim] = Receivers.SourcesFlat_.data() + NumReceivers*iDim;
    }
    Receivers.SourceRanks_.resize(NumReceivers);

    for (long long iReceiver = 0; iReceiver < NumReceivers; ++iReceiver) {
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Receivers.Points_[iDim][iReceiver] = 0;
      }
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Receivers.Sources_[iDim][iReceiver] = 0;
      }
      Receivers.SourceRanks_[iReceiver] = -1;
    }

  }

  Receivers.Edits_.Count_ = true;
  Receivers.Edits_.Points_ = true;
  Receivers.Edits_.Sources_ = true;

  MPI_Barrier(Receivers.Comm_);

}

void GetReceiverPoints(const connectivity_r &Receivers, int Dimension, const int *&Points) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  Points = Receivers.Points_[Dimension];

}

void EditReceiverPoints(connectivity_r &Receivers, int Dimension, int *&Points) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  bool StartEdit = Receivers.PointsEditRefCount_ == 0;
  ++Receivers.PointsEditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Receivers.Comm_);
  }

  Points = Receivers.Points_[Dimension];

}

void ReleaseReceiverPoints(connectivity_r &Receivers, int Dimension, int *&Points) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Points == Receivers.Points_[Dimension], "Invalid points pointer.");
  OVK_DEBUG_ASSERT(EditingPoints(Receivers), "Unable to release points; not currently being "
    "edited.");

  --Receivers.PointsEditRefCount_;
  bool EndEdit = Receivers.PointsEditRefCount_ == 0;

  Points = nullptr;

  if (EndEdit) {
    Receivers.Edits_.Points_ = true;
    MPI_Barrier(Receivers.Comm_);
  }

}

void GetReceiverSources(const connectivity_r &Receivers, int Dimension, const int *&Sources) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  Sources = Receivers.Sources_[Dimension];

}

void EditReceiverSources(connectivity_r &Receivers, int Dimension, int *&Sources) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  bool StartEdit = Receivers.SourcesEditRefCount_ == 0;
  ++Receivers.SourcesEditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Receivers.Comm_);
  }

  Sources = Receivers.Sources_[Dimension];

}

void ReleaseReceiverSources(connectivity_r &Receivers, int Dimension, int *&Sources) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Sources == Receivers.Sources_[Dimension], "Invalid sources pointer.");
  OVK_DEBUG_ASSERT(EditingSources(Receivers), "Unable to release sources; not currently being "
    "edited.");

  --Receivers.SourcesEditRefCount_;
  bool EndEdit = Receivers.SourcesEditRefCount_ == 0;

  Sources = nullptr;

  if (EndEdit) {
    Receivers.Edits_.Sources_ = true;
    MPI_Barrier(Receivers.Comm_);
  }

}

void GetReceiverSourceRanks(const connectivity_r &Receivers, const int *&SourceRanks) {

  SourceRanks = Receivers.SourceRanks_.data();

}

void EditReceiverSourceRanks(connectivity_r &Receivers, int *&SourceRanks) {

  bool StartEdit = Receivers.SourceRanksEditRefCount_ == 0;
  ++Receivers.SourceRanksEditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Receivers.Comm_);
  }

  SourceRanks = Receivers.SourceRanks_.data();

}

void ReleaseReceiverSourceRanks(connectivity_r &Receivers, int *&SourceRanks) {

  OVK_DEBUG_ASSERT(SourceRanks == Receivers.SourceRanks_.data(), "Invalid source ranks pointer.");
  OVK_DEBUG_ASSERT(EditingSourceRanks(Receivers), "Unable to release source ranks; not currently "
    "being edited.");

  --Receivers.SourceRanksEditRefCount_;
  bool EndEdit = Receivers.SourceRanksEditRefCount_ == 0;

  SourceRanks = nullptr;

  if (EndEdit) {
    Receivers.Edits_.Sources_ = true;
    MPI_Barrier(Receivers.Comm_);
  }

}

namespace {

bool EditingPoints(const connectivity_r &Receivers) {

  return Receivers.PointsEditRefCount_ > 0;

}

bool EditingSources(const connectivity_r &Receivers) {

  return Receivers.SourcesEditRefCount_ > 0;

}

bool EditingSourceRanks(const connectivity_r &Receivers) {

  return Receivers.SourceRanksEditRefCount_ > 0;

}

}

namespace core {

void GetConnectivityReceiverSideEdits(const connectivity_r &Receivers, const connectivity_r::edits
  *&Edits) {

  Edits = &Receivers.Edits_;

}

void ResetConnectivityReceiverSideEdits(connectivity_r &Receivers) {

  MPI_Barrier(Receivers.Comm_);

  DefaultEdits(Receivers.Edits_);

  MPI_Barrier(Receivers.Comm_);

}

}

namespace {

void DefaultEdits(connectivity_r::edits &Edits) {

  Edits.Count_ = false;
  Edits.Points_ = false;
  Edits.Sources_ = false;

}

}

}
