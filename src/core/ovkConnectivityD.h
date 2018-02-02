// Copyright (c) 2017 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_PUBLIC_CONNECTIVITY_D_INCLUDED
#define OVK_CORE_PUBLIC_CONNECTIVITY_D_INCLUDED

#include <ovkGlobal.h>
#include <ovkGrid.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ovk_connectivity_d_properties;
typedef struct ovk_connectivity_d_properties ovk_connectivity_d_properties;

struct ovk_connectivity_d;
typedef struct ovk_connectivity_d ovk_connectivity_d;

void ovkGetConnectivityDonorSideProperties(const ovk_connectivity_d *Donors,
  const ovk_connectivity_d_properties **Properties);

void ovkResizeDonors(ovk_connectivity_d *Donors, size_t NumDonors, int MaxDonorSize);
void ovkEditDonorExtents(ovk_connectivity_d *Donors, int Dimension, int **Begins, int **Ends);
void ovkReleaseDonorExtents(ovk_connectivity_d *Donors, int Dimension, int **Begins, int **Ends);
void ovkEditDonorCoords(ovk_connectivity_d *Donors, int Dimension, double **Coords);
void ovkReleaseDonorCoords(ovk_connectivity_d *Donors, int Dimension, double **Coords);
void ovkEditDonorInterpCoefs(ovk_connectivity_d *Donors, int Dimension, int Index,
  double **InterpCoefs);
void ovkReleaseDonorInterpCoefs(ovk_connectivity_d *Donors, int Dimension, int Index,
  double **InterpCoefs);
void ovkEditDonorDestinations(ovk_connectivity_d *Donors, int Dimension, int **Destinations);
void ovkReleaseDonorDestinations(ovk_connectivity_d *Donors, int Dimension, int **Destinations);
void ovkEditDonorDestinationRanks(ovk_connectivity_d *Donors, int **DestinationRanks);
void ovkReleaseDonorDestinationRanks(ovk_connectivity_d *Donors, int **DestinationRanks);

void ovkGetConnectivityDonorSideGrid(const ovk_connectivity_d *Donors, const ovk_grid **DonorGrid);

void ovkGetConnectivityDonorSidePropertyGridID(const ovk_connectivity_d_properties *Properties,
  int *GridID);
void ovkGetConnectivityDonorSidePropertyDimension(const ovk_connectivity_d_properties *Properties,
  int *NumDims);
void ovkGetConnectivityDonorSidePropertyComm(const ovk_connectivity_d_properties *Properties,
  MPI_Comm *Comm);
void ovkGetConnectivityDonorSidePropertyCommSize(const ovk_connectivity_d_properties *Properties,
  int *CommSize);
void ovkGetConnectivityDonorSidePropertyCommRank(const ovk_connectivity_d_properties *Properties,
  int *CommRank);
void ovkGetConnectivityDonorSidePropertyDonorCount(const ovk_connectivity_d_properties *Properties,
  size_t *NumDonors);
void ovkGetConnectivityDonorSidePropertyMaxDonorSize(const ovk_connectivity_d_properties *Properties,
  int *MaxDonorSize);

#ifdef __cplusplus
}
#endif

#endif
