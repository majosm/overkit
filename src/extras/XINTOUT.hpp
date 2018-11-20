// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_EXTRAS_XINTOUT_HPP_INCLUDED
#define OVK_EXTRAS_XINTOUT_HPP_INCLUDED

#include <ovk/extras/Constants.hpp>
#include <ovk/extras/Global.hpp>
#include <ovk/core/Domain.hpp>

#include <mpi.h>

#include <string>

namespace ovk {

error ImportXINTOUT(domain &Domain, const std::string &HOPath, const std::string &XPath,
  int ReadGranularityAdjust, MPI_Info MPIInfo);
error ExportXINTOUT(const domain &Domain, const std::string &HOPath, const std::string &XPath,
  xintout_format Format, endian Endian, int WriteGranularityAdjust, MPI_Info MPIInfo);

}

#endif
