// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_C_XINTOUT_H_INCLUDED
#define OVK_EXTRAS_C_XINTOUT_H_INCLUDED

#include <ovk/extras-c/Constants.h>
#include <ovk/extras-c/Global.h>
#include <ovk/core-c/Constants.h>
#include <ovk/core-c/Domain.h>

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

void ovkImportXINTOUT(ovk_domain *Domain, int ConnectivityComponentID, const char *HOPath, const
  char *XPath, int ReadGranularityAdjust, MPI_Info MPIInfo, ovk_error *Error);
void ovkExportXINTOUT(const ovk_domain *Domain, int ConnectivityComponentID, const char *HOPath,
  const char *XPath, ovk_xintout_format Format, ovk_endian Endian, int WriteGranularityAdjust,
  MPI_Info MPIInfo, ovk_error *Error);

#ifdef __cplusplus
}
#endif

#endif
