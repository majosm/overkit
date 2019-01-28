// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/ConnectivityD.hpp"

#include "ovk/core/Comm.hpp"
#include "ovk/core/Constants.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/ErrorHandler.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Grid.hpp"
#include "ovk/core/Logger.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

#include <vector>

namespace ovk {

namespace {

bool EditingExtents(const connectivity_d &Donors);
bool EditingCoords(const connectivity_d &Donors);
bool EditingInterpCoefs(const connectivity_d &Donors);
bool EditingDestinations(const connectivity_d &Donors);
bool EditingDestinationRanks(const connectivity_d &Donors);

void DefaultEdits(connectivity_d::edits &Edits);

}

namespace core {

void CreateConnectivityDonorSide(connectivity_d &Donors, const grid &Grid, int DestinationGridID,
  core::logger &Logger, core::error_handler &ErrorHandler) {

  Donors.Comm_ = core::GetGridComm(Grid);

  MPI_Barrier(Donors.Comm_);

  Donors.Logger_ = &Logger;
  Donors.ErrorHandler_ = &ErrorHandler;

  GetGridID(Grid, Donors.GridID_);
  GetGridDimension(Grid, Donors.NumDims_);

  Donors.DestinationGridID_ = DestinationGridID;
  Donors.Count_ = 0;
  Donors.MaxSize_ = 1;

  Donors.Grid_ = &Grid;

  DefaultEdits(Donors.Edits_);

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors.Extents_[0][iDim] = nullptr;
    Donors.Extents_[1][iDim] = nullptr;
  }
  Donors.ExtentsEditRefCount_ = 0;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors.Coords_[iDim] = nullptr;
  }
  Donors.CoordsEditRefCount_ = 0;

  Donors.InterpCoefPtrsFlat_.resize(MAX_DIMS, nullptr);
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors.InterpCoefs_[iDim] = Donors.InterpCoefPtrsFlat_.data() + iDim;
  }
  Donors.InterpCoefsEditRefCount_ = 0;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors.Destinations_[iDim] = nullptr;
  }
  Donors.DestinationsEditRefCount_ = 0;

  Donors.DestinationRanksEditRefCount_ = 0;

  MPI_Barrier(Donors.Comm_);

}

void DestroyConnectivityDonorSide(connectivity_d &Donors) {

  MPI_Barrier(Donors.Comm_);

  Donors.ExtentsFlat_.clear();
  Donors.CoordsFlat_.clear();
  Donors.InterpCoefPtrsFlat_.clear();
  Donors.InterpCoefsFlat_.clear();
  Donors.DestinationsFlat_.clear();
  Donors.DestinationRanks_.clear();

  MPI_Barrier(Donors.Comm_);

  Donors.Comm_.Reset();

}

}

void GetConnectivityDonorSideGridID(const connectivity_d &Donors, int &GridID) {

  GridID = Donors.GridID_;

}

void GetConnectivityDonorSideDestinationGridID(const connectivity_d &Donors, int
  &DestinationGridID) {

  DestinationGridID = Donors.DestinationGridID_;

}

void GetConnectivityDonorSideDimension(const connectivity_d &Donors, int &NumDims) {

  NumDims = Donors.NumDims_;

}

void GetConnectivityDonorSideComm(const connectivity_d &Donors, MPI_Comm &Comm) {

  Comm = Donors.Comm_.Get();

}

void GetConnectivityDonorSideCommSize(const connectivity_d &Donors, int &CommSize) {

  CommSize = Donors.Comm_.Size();

}

void GetConnectivityDonorSideCommRank(const connectivity_d &Donors, int &CommRank) {

  CommRank = Donors.Comm_.Rank();

}

void GetConnectivityDonorSideCount(const connectivity_d &Donors, long long &NumDonors) {

  NumDonors = Donors.Count_;

}

void GetConnectivityDonorSideMaxSize(const connectivity_d &Donors, int &MaxSize) {

  MaxSize = Donors.MaxSize_;

}

void GetConnectivityDonorSideGrid(const connectivity_d &Donors, const grid *&Grid) {

  Grid = Donors.Grid_;

}

void ResizeDonors(connectivity_d &Donors, long long NumDonors, int MaxSize) {

  OVK_DEBUG_ASSERT(NumDonors >= 0, "Invalid donor count.");
  OVK_DEBUG_ASSERT(MaxSize > 0, "Invalid max size.");

  MPI_Barrier(Donors.Comm_);

  OVK_DEBUG_ASSERT(!EditingExtents(Donors), "Cannot resize donors while editing extents.");
  OVK_DEBUG_ASSERT(!EditingCoords(Donors), "Cannot resize donors while editing coords.");
  OVK_DEBUG_ASSERT(!EditingInterpCoefs(Donors), "Cannot resize donors while editing interp coefs.");
  OVK_DEBUG_ASSERT(!EditingDestinations(Donors), "Cannot resize donors while editing destinations.");
  OVK_DEBUG_ASSERT(!EditingDestinationRanks(Donors), "Cannot resize donors while editing "
    "destination ranks.");

  // Pretty sure this isn't needed anymore, since MaxSize is now required to be greater than 0, and
  // edit only blocks on the first call
  // if (OVK_DEBUG) {
  //   // Needed because editing interp coefs blocks on comm -- if max size was 0 on some ranks,
  //   // calling the edit function in a loop from 0 to max size would result in it not being called
  //   // on those ranks
  //   int GlobalMaxSize;
  //   MPI_Allreduce(&MaxSize, &GlobalMaxSize, 1, MPI_INT, MPI_MAX, Donors->comm);
  //   OVK_DEBUG_ASSERT(MaxSize == GlobalMaxSize, "Max size must be the same on all connectivity "
  //     "processes.");
  // }

  int NumDims = Donors.NumDims_;

  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors.Extents_[0][iDim] = nullptr;
    Donors.Extents_[1][iDim] = nullptr;
  }
  Donors.ExtentsFlat_.clear();
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors.Coords_[iDim] = nullptr;
  }
  Donors.CoordsFlat_.clear();
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors.InterpCoefs_[iDim] = nullptr;
  }
  Donors.InterpCoefPtrsFlat_.clear();
  Donors.InterpCoefsFlat_.clear();
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors.Destinations_[iDim] = nullptr;
  }
  Donors.DestinationsFlat_.clear();
  Donors.DestinationRanks_.clear();

  Donors.Count_ = NumDonors;
  Donors.MaxSize_ = MaxSize;

  Donors.InterpCoefPtrsFlat_.resize(MaxSize*MAX_DIMS, nullptr);
  for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
    Donors.InterpCoefs_[iDim] = Donors.InterpCoefPtrsFlat_.data() + MaxSize*iDim;
  }

  if (NumDonors > 0) {

    Donors.ExtentsFlat_.resize(2*MAX_DIMS*NumDonors);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Donors.Extents_[0][iDim] = Donors.ExtentsFlat_.data() + iDim*NumDonors;
      Donors.Extents_[1][iDim] = Donors.ExtentsFlat_.data() + (MAX_DIMS+iDim)*NumDonors;
    }
    Donors.CoordsFlat_.resize(MAX_DIMS*NumDonors);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Donors.Coords_[iDim] = Donors.CoordsFlat_.data() + iDim*NumDonors;
    }
    Donors.InterpCoefsFlat_.resize(MAX_DIMS*MaxSize*NumDonors);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      for (int iPoint = 0; iPoint < MaxSize; ++iPoint) {
        Donors.InterpCoefs_[iDim][iPoint] = Donors.InterpCoefsFlat_.data() + (iDim*MaxSize + iPoint)
          * NumDonors;
      }
    }
    Donors.DestinationsFlat_.resize(MAX_DIMS*NumDonors);
    for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
      Donors.Destinations_[iDim] = Donors.DestinationsFlat_.data() + iDim*NumDonors;
    }
    Donors.DestinationRanks_.resize(NumDonors);

    for (long long iDonor = 0; iDonor < NumDonors; ++iDonor) {
      for (int iDim = 0; iDim < NumDims; ++iDim) {
        Donors.Extents_[0][iDim][iDonor] = 0;
        Donors.Extents_[1][iDim][iDonor] = 0;
      }
      for (int iDim = NumDims; iDim < MAX_DIMS; ++iDim) {
        Donors.Extents_[0][iDim][iDonor] = 0;
        Donors.Extents_[1][iDim][iDonor] = 1;
      }
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Donors.Coords_[iDim][iDonor] = 0.;
      }
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        for (int iPoint = 0; iPoint < MaxSize; ++iPoint) {
          Donors.InterpCoefs_[iDim][iPoint][iDonor] = 0.;
        }
      }
      for (int iDim = 0; iDim < MAX_DIMS; ++iDim) {
        Donors.Destinations_[iDim][iDonor] = 0;
      }
      Donors.DestinationRanks_[iDonor] = -1;
    }

  }

  Donors.Edits_.Count_ = true;
  Donors.Edits_.Extents_ = true;
  Donors.Edits_.Coords_ = true;
  Donors.Edits_.InterpCoefs_ = true;
  Donors.Edits_.Destinations_ = true;

  MPI_Barrier(Donors.Comm_);

}

void GetDonorExtents(const connectivity_d &Donors, int Dimension, const int *&Begins, const int
  *&Ends) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  Begins = Donors.Extents_[0][Dimension];
  Ends = Donors.Extents_[1][Dimension];

}

void EditDonorExtents(connectivity_d &Donors, int Dimension, int *&Begins, int *&Ends) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  bool StartEdit = Donors.ExtentsEditRefCount_ == 0;
  ++Donors.ExtentsEditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Donors.Comm_);
  }

  Begins = Donors.Extents_[0][Dimension];
  Ends = Donors.Extents_[1][Dimension];

}

void ReleaseDonorExtents(connectivity_d &Donors, int Dimension, int *&Begins, int *&Ends) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Begins == Donors.Extents_[0][Dimension], "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(Ends == Donors.Extents_[1][Dimension], "Invalid ends pointer.");
  OVK_DEBUG_ASSERT(EditingExtents(Donors), "Unable to release extents; not currently being edited.");

  --Donors.ExtentsEditRefCount_;
  bool EndEdit = Donors.ExtentsEditRefCount_ == 0;

  Begins = nullptr;
  Ends = nullptr;

  if (EndEdit) {
    Donors.Edits_.Extents_ = true;
    MPI_Barrier(Donors.Comm_);
  }

}

void GetDonorCoords(const connectivity_d &Donors, int Dimension, const double *&Coords) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  Coords = Donors.Coords_[Dimension];

}

void EditDonorCoords(connectivity_d &Donors, int Dimension, double *&Coords) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  bool StartEdit = Donors.CoordsEditRefCount_ == 0;
  ++Donors.CoordsEditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Donors.Comm_);
  }

  Coords = Donors.Coords_[Dimension];

}

void ReleaseDonorCoords(connectivity_d &Donors, int Dimension, double *&Coords) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Coords == Donors.Coords_[Dimension], "Invalid coords pointer.");
  OVK_DEBUG_ASSERT(EditingCoords(Donors), "Unable to release coords; not currently being edited.");

  --Donors.CoordsEditRefCount_;
  bool EndEdit = Donors.CoordsEditRefCount_ == 0;

  Coords = nullptr;

  if (EndEdit) {
    Donors.Edits_.Coords_ = true;
    MPI_Barrier(Donors.Comm_);
  }

}

void GetDonorInterpCoefs(const connectivity_d &Donors, int Dimension, int Point, const double
  *&InterpCoefs) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Point >= 0 && Point < Donors.MaxSize_, "Invalid point.");

  InterpCoefs = Donors.InterpCoefs_[Dimension][Point];

}

void EditDonorInterpCoefs(connectivity_d &Donors, int Dimension, int Point, double *&InterpCoefs) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Point >= 0 && Point < Donors.MaxSize_, "Invalid point.");

  bool StartEdit = Donors.InterpCoefsEditRefCount_ == 0;
  ++Donors.InterpCoefsEditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Donors.Comm_);
  }

  InterpCoefs = Donors.InterpCoefs_[Dimension][Point];

}

void ReleaseDonorInterpCoefs(connectivity_d &Donors, int Dimension, int Point, double *&InterpCoefs)
  {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Point >= 0 && Point < Donors.MaxSize_, "Invalid point.");
  OVK_DEBUG_ASSERT(InterpCoefs == Donors.InterpCoefs_[Dimension][Point], "Invalid interp coefs "
    "pointer.");
  OVK_DEBUG_ASSERT(EditingInterpCoefs(Donors), "Unable to release interp coefs; not currently "
    "being edited.");

  --Donors.InterpCoefsEditRefCount_;
  bool EndEdit = Donors.InterpCoefsEditRefCount_ == 0;

  InterpCoefs = nullptr;

  if (EndEdit) {
    Donors.Edits_.InterpCoefs_ = true;
    MPI_Barrier(Donors.Comm_);
  }

}

void GetDonorDestinations(const connectivity_d &Donors, int Dimension, const int *&Destinations) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  Destinations = Donors.Destinations_[Dimension];

}

void EditDonorDestinations(connectivity_d &Donors, int Dimension, int *&Destinations) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");

  bool StartEdit = Donors.DestinationsEditRefCount_ == 0;
  ++Donors.DestinationsEditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Donors.Comm_);
  }

  Destinations = Donors.Destinations_[Dimension];

}

void ReleaseDonorDestinations(connectivity_d &Donors, int Dimension, int *&Destinations) {

  OVK_DEBUG_ASSERT(Dimension >= 0 && Dimension < MAX_DIMS, "Invalid dimension.");
  OVK_DEBUG_ASSERT(Destinations == Donors.Destinations_[Dimension], "Invalid destinations pointer.");
  OVK_DEBUG_ASSERT(EditingDestinations(Donors), "Unable to release destinations; not currently "
    "being edited.");

  --Donors.DestinationsEditRefCount_;
  bool EndEdit = Donors.DestinationsEditRefCount_ == 0;

  Destinations = nullptr;

  if (EndEdit) {
    Donors.Edits_.Destinations_ = true;
    MPI_Barrier(Donors.Comm_);
  }


}

void GetDonorDestinationRanks(const connectivity_d &Donors, const int *&DestinationRanks) {

  DestinationRanks = Donors.DestinationRanks_.data();

}

void EditDonorDestinationRanks(connectivity_d &Donors, int *&DestinationRanks) {

  bool StartEdit = Donors.DestinationRanksEditRefCount_ == 0;
  ++Donors.DestinationRanksEditRefCount_;

  if (StartEdit) {
    MPI_Barrier(Donors.Comm_);
  }

  DestinationRanks = Donors.DestinationRanks_.data();

}

void ReleaseDonorDestinationRanks(connectivity_d &Donors, int *&DestinationRanks) {

  OVK_DEBUG_ASSERT(DestinationRanks == Donors.DestinationRanks_.data(), "Invalid destination ranks "
    "pointer.");
  OVK_DEBUG_ASSERT(EditingDestinationRanks(Donors), "Unable to release destination ranks; not "
    "currently being edited.");

  --Donors.DestinationRanksEditRefCount_;
  bool EndEdit = Donors.DestinationRanksEditRefCount_ == 0;

  DestinationRanks = nullptr;

  if (EndEdit) {
    Donors.Edits_.Destinations_ = true;
    MPI_Barrier(Donors.Comm_);
  }

}

namespace {

bool EditingExtents(const connectivity_d &Donors) {

  return Donors.ExtentsEditRefCount_ > 0;

}

bool EditingCoords(const connectivity_d &Donors) {

  return Donors.CoordsEditRefCount_ > 0;

}

bool EditingInterpCoefs(const connectivity_d &Donors) {

  return Donors.InterpCoefsEditRefCount_ > 0;

}

bool EditingDestinations(const connectivity_d &Donors) {

  return Donors.DestinationsEditRefCount_ > 0;

}

bool EditingDestinationRanks(const connectivity_d &Donors) {

  return Donors.DestinationRanksEditRefCount_ > 0;

}

}

namespace core {

void GetConnectivityDonorSideEdits(const connectivity_d &Donors, const connectivity_d::edits
  *&Edits) {

  Edits = &Donors.Edits_;

}

void ResetConnectivityDonorSideEdits(connectivity_d &Donors) {

  MPI_Barrier(Donors.Comm_);

  DefaultEdits(Donors.Edits_);

  MPI_Barrier(Donors.Comm_);

}

}

namespace {

void DefaultEdits(connectivity_d::edits &Edits) {

  Edits.Count_ = false;
  Edits.Extents_ = false;
  Edits.Coords_ = false;
  Edits.InterpCoefs_ = false;
  Edits.Destinations_ = false;

}

}

}
