// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONNECTIVITY_R_HPP_INCLUDED
#define OVK_CORE_CONNECTIVITY_R_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/ErrorHandler.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Logger.hpp>

#include <mpi.h>

namespace ovk {

struct connectivity_r {

  struct edits {
    bool Count_;
    bool Points_;
    bool Sources_;
  };

  mutable core::logger *Logger_;
  mutable core::error_handler *ErrorHandler_;
  int NumDims_;
  core::comm_view Comm_;
  const grid *Grid_;
  int SourceGridID_;
  long long Count_;
  edits Edits_;
  array<int,2> Points_;
  int PointsEditRefCount_;
  array<int,2> Sources_;
  int SourcesEditRefCount_;
  array<int> SourceRanks_;
  int SourceRanksEditRefCount_;

};

namespace core {
void CreateConnectivityReceiverSide(connectivity_r &Receivers, const grid &Grid, int SourceGridID,
  core::logger &Logger, core::error_handler &ErrorHandler);
void DestroyConnectivityReceiverSide(connectivity_r &Receivers);
}

void GetConnectivityReceiverSideGridID(const connectivity_r &Receivers, int &GridID);
void GetConnectivityReceiverSideSourceGridID(const connectivity_r &Receivers, int &SourceGridID);
void GetConnectivityReceiverSideDimension(const connectivity_r &Receivers, int &NumDims);
void GetConnectivityReceiverSideComm(const connectivity_r &Receivers, MPI_Comm &Comm);
void GetConnectivityReceiverSideCommSize(const connectivity_r &Receivers, int &CommSize);
void GetConnectivityReceiverSideCommRank(const connectivity_r &Receivers, int &CommRank);
void GetConnectivityReceiverSideCount(const connectivity_r &Receivers, long long &NumReceivers);
void GetConnectivityReceiverSideGrid(const connectivity_r &Receivers, const grid *&Grid);

void ResizeReceivers(connectivity_r &Receivers, long long NumReceivers);

void GetReceiverPoints(const connectivity_r &Receivers, int Dimension, const int *&Points);
void EditReceiverPoints(connectivity_r &Receivers, int Dimension, int *&Points);
void ReleaseReceiverPoints(connectivity_r &Receivers, int Dimension, int *&Points);

void GetReceiverSources(const connectivity_r &Receivers, int iDim, const int *&Sources);
void EditReceiverSources(connectivity_r &Receivers, int iDim, int *&Sources);
void ReleaseReceiverSources(connectivity_r &Receivers, int iDim, int *&Sources);

void GetReceiverSourceRanks(const connectivity_r &Receivers, const int *&SourceRanks);
void EditReceiverSourceRanks(connectivity_r &Receivers, int *&SourceRanks);
void ReleaseReceiverSourceRanks(connectivity_r &Receivers, int *&SourceRanks);

namespace core {
void GetConnectivityReceiverSideEdits(const connectivity_r &Receivers, const connectivity_r::edits
  *&Edits);
void ResetConnectivityReceiverSideEdits(connectivity_r &Receivers);
}

}

#endif
