// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/extras-c/XINTOUT.h"

#include "ovk/extras-c/Constants.h"
#include "ovk/extras-c/Global.h"
#include "ovk/extras/Constants.hpp"
#include "ovk/extras/Global.hpp"
#include "ovk/extras/XINTOUT.hpp"
#include "ovk/core-c/Constants.h"
#include "ovk/core-c/Global.h"
#include "ovk/core/Debug.hpp"

#include <string>

extern "C" {

ovk_error ovkImportXINTOUT(ovk_domain *Domain, const char *HOPath, const char *XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo) {

  OVK_DEBUG_ASSERT(Domain, "Invalid exchange pointer.");
  OVK_DEBUG_ASSERT(HOPath, "Invalid ho path pointer.");
  OVK_DEBUG_ASSERT(XPath, "Invalid x path pointer.");

  auto &DomainCPP = *reinterpret_cast<ovk::domain *>(Domain);
  std::string HOPathCPP = HOPath;
  std::string XPathCPP = XPath;

  ovk::error Error = ovk::ImportXINTOUT(DomainCPP, HOPathCPP, XPathCPP, ReadGranularityAdjust,
    MPIInfo);

  return ovk_error(Error);

}

ovk_error ovkExportXINTOUT(const ovk_domain *Domain, const char *HOPath, const char *XPath,
  ovk_xintout_format Format, ovk_endian Endian, int WriteGranularityAdjust, MPI_Info MPIInfo) {

  OVK_DEBUG_ASSERT(Domain, "Invalid exchange pointer.");
  OVK_DEBUG_ASSERT(HOPath, "Invalid ho path pointer.");
  OVK_DEBUG_ASSERT(XPath, "Invalid x path pointer.");

  auto &DomainCPP = *reinterpret_cast<const ovk::domain *>(Domain);
  std::string HOPathCPP = HOPath;
  std::string XPathCPP = XPath;

  ovk::error Error = ovk::ExportXINTOUT(DomainCPP, HOPathCPP, XPathCPP, ovk::xintout_format(Format),
    ovk::endian(Endian), WriteGranularityAdjust, MPIInfo);

  return ovk_error(Error);

}

}
