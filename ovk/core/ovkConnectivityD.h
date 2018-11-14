// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CONNECTIVITY_D_INCLUDED
#define OVK_CORE_PUBLIC_CONNECTIVITY_D_INCLUDED

#include <ovk/core/ovkGlobal.h>
#include <ovk/core/ovkGrid.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity_d;
typedef struct ovk_connectivity_d ovk_connectivity_d;

void ovkGetConnectivityDonorSideGridID(const ovk_connectivity_d *Donors, int *GridID);
void ovkGetConnectivityDonorSideDestinationGridID(const ovk_connectivity_d *Donors, int *GridID);
void ovkGetConnectivityDonorSideDimension(const ovk_connectivity_d *Donors, int *NumDims);
void ovkGetConnectivityDonorSideComm(const ovk_connectivity_d *Donors, MPI_Comm *Comm);
void ovkGetConnectivityDonorSideCommSize(const ovk_connectivity_d *Donors, int *CommSize);
void ovkGetConnectivityDonorSideCommRank(const ovk_connectivity_d *Donors, int *CommRank);
void ovkGetConnectivityDonorSideCount(const ovk_connectivity_d *Donors, size_t *NumDonors);
void ovkGetConnectivityDonorSideMaxSize(const ovk_connectivity_d *Donors, int *MaxSize);

void ovkResizeDonors(ovk_connectivity_d *Donors, size_t NumDonors, int MaxSize);

void ovkGetDonorExtents(const ovk_connectivity_d *Donors, int Dimension, const int **Begins,
  const int **Ends);
void ovkEditDonorExtents(ovk_connectivity_d *Donors, int Dimension, int **Begins, int **Ends);
void ovkReleaseDonorExtents(ovk_connectivity_d *Donors, int Dimension, int **Begins, int **Ends);

void ovkGetDonorCoords(const ovk_connectivity_d *Donors, int Dimension, const double **Coords);
void ovkEditDonorCoords(ovk_connectivity_d *Donors, int Dimension, double **Coords);
void ovkReleaseDonorCoords(ovk_connectivity_d *Donors, int Dimension, double **Coords);

void ovkGetDonorInterpCoefs(const ovk_connectivity_d *Donors, int Dimension, int Point,
  const double **InterpCoefs);
void ovkEditDonorInterpCoefs(ovk_connectivity_d *Donors, int Dimension, int Point,
  double **InterpCoefs);
void ovkReleaseDonorInterpCoefs(ovk_connectivity_d *Donors, int Dimension, int Point,
  double **InterpCoefs);

void ovkGetDonorDestinations(const ovk_connectivity_d *Donors, int Dimension,
  const int **Destinations);
void ovkEditDonorDestinations(ovk_connectivity_d *Donors, int Dimension, int **Destinations);
void ovkReleaseDonorDestinations(ovk_connectivity_d *Donors, int Dimension, int **Destinations);

void ovkGetDonorDestinationRanks(const ovk_connectivity_d *Donors, const int **DestinationRanks);
void ovkEditDonorDestinationRanks(ovk_connectivity_d *Donors, int **DestinationRanks);
void ovkReleaseDonorDestinationRanks(ovk_connectivity_d *Donors, int **DestinationRanks);

void ovkGetConnectivityDonorSideGrid(const ovk_connectivity_d *Donors, const ovk_grid **DonorGrid);

#ifdef __cplusplus
}
#endif

#endif
