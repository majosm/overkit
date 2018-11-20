// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core-c/ConnectivityD.h"

#include "ovk/core-c/Global.h"
#include "ovk/core-c/Grid.h"
#include "ovk/core/ConnectivityD.hpp"
#include "ovk/core/Debug.hpp"
#include "ovk/core/Global.hpp"

#include <mpi.h>

void ovkGetConnectivityDonorSideGridID(const ovk_connectivity_d *Donors, int *GridID) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(GridID, "Invalid grid ID pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetConnectivityDonorSideGridID(DonorsCPP, *GridID);

}

void ovkGetConnectivityDonorSideDestinationGridID(const ovk_connectivity_d *Donors, int
  *DestinationGridID) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(DestinationGridID, "Invalid destination grid ID pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetConnectivityDonorSideDestinationGridID(DonorsCPP, *DestinationGridID);

}

void ovkGetConnectivityDonorSideDimension(const ovk_connectivity_d *Donors, int *NumDims) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(NumDims, "Invalid num dims pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetConnectivityDonorSideDimension(DonorsCPP, *NumDims);

}

void ovkGetConnectivityDonorSideComm(const ovk_connectivity_d *Donors, MPI_Comm *Comm) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Comm, "Invalid comm pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetConnectivityDonorSideComm(DonorsCPP, *Comm);

}

void ovkGetConnectivityDonorSideCommSize(const ovk_connectivity_d *Donors, int *CommSize) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(CommSize, "Invalid comm size pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetConnectivityDonorSideCommSize(DonorsCPP, *CommSize);

}

void ovkGetConnectivityDonorSideCommRank(const ovk_connectivity_d *Donors, int *CommRank) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(CommRank, "Invalid comm rank pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetConnectivityDonorSideCommRank(DonorsCPP, *CommRank);

}

void ovkGetConnectivityDonorSideCount(const ovk_connectivity_d *Donors, long long *NumDonors) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(NumDonors, "Invalid num donors pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetConnectivityDonorSideCount(DonorsCPP, *NumDonors);

}

void ovkGetConnectivityDonorSideMaxSize(const ovk_connectivity_d *Donors, int *MaxSize) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(MaxSize, "Invalid max size pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetConnectivityDonorSideMaxSize(DonorsCPP, *MaxSize);

}

void ovkGetConnectivityDonorSideGrid(const ovk_connectivity_d *Donors, const ovk_grid **Grid) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Grid, "Invalid donor grid pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);

  const ovk::grid *GridCPPPtr;
  ovk::GetConnectivityDonorSideGrid(DonorsCPP, GridCPPPtr);

  *Grid = reinterpret_cast<const ovk_grid *>(GridCPPPtr);

}

void ovkResizeDonors(ovk_connectivity_d *Donors, long long NumDonors, int MaxSize) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::ResizeDonors(DonorsCPP, NumDonors, MaxSize);

}

void ovkGetDonorExtents(const ovk_connectivity_d *Donors, int Dimension, const int **Begins,
  const int **Ends) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Begins, "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(Ends, "Invalid ends pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetDonorExtents(DonorsCPP, Dimension, *Begins, *Ends);

}

void ovkEditDonorExtents(ovk_connectivity_d *Donors, int Dimension, int **Begins, int **Ends) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Begins, "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(Ends, "Invalid ends pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::EditDonorExtents(DonorsCPP, Dimension, *Begins, *Ends);

}

void ovkReleaseDonorExtents(ovk_connectivity_d *Donors, int Dimension, int **Begins, int **Ends) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Begins, "Invalid begins pointer.");
  OVK_DEBUG_ASSERT(Ends, "Invalid ends pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::ReleaseDonorExtents(DonorsCPP, Dimension, *Begins, *Ends);

}

void ovkGetDonorCoords(const ovk_connectivity_d *Donors, int Dimension, const double **Coords) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetDonorCoords(DonorsCPP, Dimension, *Coords);

}

void ovkEditDonorCoords(ovk_connectivity_d *Donors, int Dimension, double **Coords) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::EditDonorCoords(DonorsCPP, Dimension, *Coords);

}

void ovkReleaseDonorCoords(ovk_connectivity_d *Donors, int Dimension, double **Coords) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Coords, "Invalid coords pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::ReleaseDonorCoords(DonorsCPP, Dimension, *Coords);

}

void ovkGetDonorInterpCoefs(const ovk_connectivity_d *Donors, int Dimension, int Point,
  const double **InterpCoefs) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(InterpCoefs, "Invalid interp coefs pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetDonorInterpCoefs(DonorsCPP, Dimension, Point, *InterpCoefs);

}

void ovkEditDonorInterpCoefs(ovk_connectivity_d *Donors, int Dimension, int Point,
  double **InterpCoefs) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(InterpCoefs, "Invalid interp coefs pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::EditDonorInterpCoefs(DonorsCPP, Dimension, Point, *InterpCoefs);

}

void ovkReleaseDonorInterpCoefs(ovk_connectivity_d *Donors, int Dimension, int Point, double
  **InterpCoefs) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(InterpCoefs, "Invalid donor coefs pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::ReleaseDonorInterpCoefs(DonorsCPP, Dimension, Point, *InterpCoefs);

}

void ovkGetDonorDestinations(const ovk_connectivity_d *Donors, int Dimension, const int
  **Destinations) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetDonorDestinations(DonorsCPP, Dimension, *Destinations);

}

void ovkEditDonorDestinations(ovk_connectivity_d *Donors, int Dimension, int **Destinations) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::EditDonorDestinations(DonorsCPP, Dimension, *Destinations);

}

void ovkReleaseDonorDestinations(ovk_connectivity_d *Donors, int Dimension, int **Destinations) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(Destinations, "Invalid destinations pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::ReleaseDonorDestinations(DonorsCPP, Dimension, *Destinations);

}

void ovkGetDonorDestinationRanks(const ovk_connectivity_d *Donors, const int **DestinationRanks) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  auto &DonorsCPP = *reinterpret_cast<const ovk::connectivity_d *>(Donors);
  ovk::GetDonorDestinationRanks(DonorsCPP, *DestinationRanks);

}

void ovkEditDonorDestinationRanks(ovk_connectivity_d *Donors, int **DestinationRanks) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::EditDonorDestinationRanks(DonorsCPP, *DestinationRanks);

}

void ovkReleaseDonorDestinationRanks(ovk_connectivity_d *Donors, int **DestinationRanks) {

  OVK_DEBUG_ASSERT(Donors, "Invalid donors pointer.");
  OVK_DEBUG_ASSERT(DestinationRanks, "Invalid destination ranks pointer.");

  auto &DonorsCPP = *reinterpret_cast<ovk::connectivity_d *>(Donors);
  ovk::ReleaseDonorDestinationRanks(DonorsCPP, *DestinationRanks);

}
