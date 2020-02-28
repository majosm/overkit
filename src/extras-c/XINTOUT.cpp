// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/extras-c/XINTOUT.h"

#include "ovk/extras-c/Global.h"
#include "ovk/extras/Global.hpp"
#include "ovk/extras/XINTOUT.hpp"
#include "ovk/core-c/Global.h"
#include "ovk/core/Debug.hpp"

#include <string>

extern "C" {

void ovkImportXINTOUT(ovk_domain *Domain, int ConnectivityComponentID, const char *HOPath, const
  char *XPath, int ReadGranularityAdjust, MPI_Info MPIInfo, ovk_error *Error) {

  OVK_DEBUG_ASSERT(Domain, "Invalid exchange pointer.");
  OVK_DEBUG_ASSERT(HOPath, "Invalid HO path pointer.");
  OVK_DEBUG_ASSERT(XPath, "Invalid X path pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  ovk::captured_error ErrorCPP;
  ovk::ImportXINTOUT(DomainCPP, ConnectivityComponentID, HOPath, XPath, ReadGranularityAdjust,
    MPIInfo, ErrorCPP);

  *Error = ovk_error(ErrorCPP.Code());

}

void ovkExportXINTOUT(const ovk_domain *Domain, int ConnectivityComponentID, const char *HOPath,
  const char *XPath, ovk_xintout_format Format, ovk_endian Endian, int WriteGranularityAdjust,
  MPI_Info MPIInfo, ovk_error *Error) {

  OVK_DEBUG_ASSERT(Domain, "Invalid exchange pointer.");
  OVK_DEBUG_ASSERT(HOPath, "Invalid HO path pointer.");
  OVK_DEBUG_ASSERT(XPath, "Invalid X path pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  ovk::captured_error ErrorCPP;
  ovk::ExportXINTOUT(DomainCPP, ConnectivityComponentID, HOPath, XPath, ovk::xintout_format(Format),
    ovk::endian(Endian), WriteGranularityAdjust, MPIInfo, ErrorCPP);

  *Error = ovk_error(ErrorCPP.Code());

}

}
