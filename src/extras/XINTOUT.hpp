// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_XINTOUT_HPP_INCLUDED
#define OVK_EXTRAS_XINTOUT_HPP_INCLUDED

#include <ovk/extras/Global.hpp>
#include <ovk/extras/XINTOUTBase.h>
#include <ovk/core/Domain.hpp>
#include <ovk/core/Error.hpp>

#include <mpi.h>

#include <string>

namespace ovk {

enum class xintout_format {
  STANDARD = OVK_XINTOUT_STANDARD,
  EXTENDED = OVK_XINTOUT_EXTENDED
};

inline bool ValidXINTOUTFormat(xintout_format Format) {
  return ovkValidXINTOUTFormat(ovk_xintout_format(Format));
}

void ImportXINTOUT(domain &Domain, int ConnectivityComponentID, const std::string &HOPath, const
  std::string &XPath, int ReadGranularityAdjust, MPI_Info MPIInfo);
void ImportXINTOUT(domain &Domain, int ConnectivityComponentID, const std::string &HOPath, const
  std::string &XPath, int ReadGranularityAdjust, MPI_Info MPIInfo, error &Error);

void ExportXINTOUT(const domain &Domain, int ConnectivityComponentID, const std::string &HOPath,
  const std::string &XPath, xintout_format Format, endian Endian, int WriteGranularityAdjust,
  MPI_Info MPIInfo);
void ExportXINTOUT(const domain &Domain, int ConnectivityComponentID, const std::string &HOPath,
  const std::string &XPath, xintout_format Format, endian Endian, int WriteGranularityAdjust,
  MPI_Info MPIInfo, error &Error);

}

#endif
