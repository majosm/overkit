// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CONNECTIVITY_D_INCLUDED
#define OVK_CORE_PUBLIC_CONNECTIVITY_D_INCLUDED

#include <ovk/core/ovkGlobal.h>
#include <ovk/core/ovkGrid.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity_d_properties;
typedef struct ovk_connectivity_d_properties ovk_connectivity_d_properties;

struct ovk_connectivity_d;
typedef struct ovk_connectivity_d ovk_connectivity_d;

void ovkGetConnectivityDonorSideProperties(const ovk_connectivity_d *Donors,
  const ovk_connectivity_d_properties **Properties);

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

void ovkGetConnectivityDonorSidePropertyGridID(const ovk_connectivity_d_properties *Properties,
  int *GridID);
void ovkGetConnectivityDonorSidePropertyDestinationGridID(const ovk_connectivity_d_properties
  *Properties, int *GridID);
void ovkGetConnectivityDonorSidePropertyDimension(const ovk_connectivity_d_properties *Properties,
  int *NumDims);
void ovkGetConnectivityDonorSidePropertyComm(const ovk_connectivity_d_properties *Properties,
  MPI_Comm *Comm);
void ovkGetConnectivityDonorSidePropertyCommSize(const ovk_connectivity_d_properties *Properties,
  int *CommSize);
void ovkGetConnectivityDonorSidePropertyCommRank(const ovk_connectivity_d_properties *Properties,
  int *CommRank);
void ovkGetConnectivityDonorSidePropertyCount(const ovk_connectivity_d_properties *Properties,
  size_t *NumDonors);
void ovkGetConnectivityDonorSidePropertyMaxSize(const ovk_connectivity_d_properties *Properties,
  int *MaxSize);

#ifdef __cplusplus
}
#endif

#endif
