// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_CONNECTIVITY_D_HPP_INCLUDED
#define OVK_CORE_CONNECTIVITY_D_HPP_INCLUDED

#include <ovk/core/Array.hpp>
#include <ovk/core/Comm.hpp>
#include <ovk/core/Constants.hpp>
#include <ovk/core/ErrorHandler.hpp>
#include <ovk/core/Global.hpp>
#include <ovk/core/Grid.hpp>
#include <ovk/core/Logger.hpp>

#include <mpi.h>

namespace ovk {

struct connectivity_d {

  struct edits {
    bool Count_;
    bool Extents_;
    bool Coords_;
    bool InterpCoefs_;
    bool Destinations_;
  };

  mutable core::logger *Logger_;
  mutable core::error_handler *ErrorHandler_;
  int GridID_;
  int DestinationGridID_;
  int NumDims_;
  core::comm Comm_;
  long long Count_;
  int MaxSize_;
  const grid *Grid_;
  edits Edits_;
  array<int,3> Extents_;
  int ExtentsEditRefCount_;
  array<double,2> Coords_;
  int CoordsEditRefCount_;
  array<double,3> InterpCoefs_;
  int InterpCoefsEditRefCount_;
  array<int,2> Destinations_;
  int DestinationsEditRefCount_;
  array<int> DestinationRanks_;
  int DestinationRanksEditRefCount_;

};

namespace core {
void CreateConnectivityDonorSide(connectivity_d &Donors, const grid &Grid, int DestinationGridID,
  core::logger &Logger, core::error_handler &ErrorHandler);
void DestroyConnectivityDonorSide(connectivity_d &Donors);
}

void GetConnectivityDonorSideGridID(const connectivity_d &Donors, int &GridID);
void GetConnectivityDonorSideDestinationGridID(const connectivity_d &Donors, int &GridID);
void GetConnectivityDonorSideDimension(const connectivity_d &Donors, int &NumDims);
void GetConnectivityDonorSideComm(const connectivity_d &Donors, MPI_Comm &Comm);
void GetConnectivityDonorSideCommSize(const connectivity_d &Donors, int &CommSize);
void GetConnectivityDonorSideCommRank(const connectivity_d &Donors, int &CommRank);
void GetConnectivityDonorSideCount(const connectivity_d &Donors, long long &NumDonors);
void GetConnectivityDonorSideMaxSize(const connectivity_d &Donors, int &MaxSize);
void GetConnectivityDonorSideGrid(const connectivity_d &Donors, const grid *&Grid);

void ResizeDonors(connectivity_d &Donors, long long NumDonors, int MaxSize);

void GetDonorExtents(const connectivity_d &Donors, int Dimension, const int *&Begins, const int
  *&Ends);
void EditDonorExtents(connectivity_d &Donors, int Dimension, int *&Begins, int *&Ends);
void ReleaseDonorExtents(connectivity_d &Donors, int Dimension, int *&Begins, int *&Ends);

void GetDonorCoords(const connectivity_d &Donors, int Dimension, const double *&Coords);
void EditDonorCoords(connectivity_d &Donors, int Dimension, double *&Coords);
void ReleaseDonorCoords(connectivity_d &Donors, int Dimension, double *&Coords);

void GetDonorInterpCoefs(const connectivity_d &Donors, int Dimension, int Point,
  const double *&InterpCoefs);
void EditDonorInterpCoefs(connectivity_d &Donors, int Dimension, int Point, double *&InterpCoefs);
void ReleaseDonorInterpCoefs(connectivity_d &Donors, int Dimension, int Point, double *&InterpCoefs);

void GetDonorDestinations(const connectivity_d &Donors, int Dimension, const int *&Destinations);
void EditDonorDestinations(connectivity_d &Donors, int Dimension, int *&Destinations);
void ReleaseDonorDestinations(connectivity_d &Donors, int Dimension, int *&Destinations);

void GetDonorDestinationRanks(const connectivity_d &Donors, const int *&DestinationRanks);
void EditDonorDestinationRanks(connectivity_d &Donors, int *&DestinationRanks);
void ReleaseDonorDestinationRanks(connectivity_d &Donors, int *&DestinationRanks);

namespace core {
void GetConnectivityDonorSideEdits(const connectivity_d &Donors, const connectivity_d::edits
  *&Edits);
void ResetConnectivityDonorSideEdits(connectivity_d &Donors);
}

}

#endif
