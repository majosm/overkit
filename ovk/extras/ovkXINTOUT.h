// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_PUBLIC_XINTOUT_INCLUDED
#define OVK_EXTRAS_PUBLIC_XINTOUT_INCLUDED

#include <ovk/extras/ovkGlobal.h>
#include <ovk/core/ovkDomain.h>

#ifdef __cplusplus
extern "C" {
#endif

ovk_error ovkEXTImportXINTOUT(ovk_domain *Domain, const char *HOPath, const char *XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo);
ovk_error ovkEXTExportXINTOUT(const ovk_domain *Domain, const char *HOPath, const char *XPath,
  ovk_ext_xintout_format Format, ovk_ext_endian Endian, int WriteGranularityAdjust,
  MPI_Info MPIInfo);

#ifdef __cplusplus
}
#endif

#endif
