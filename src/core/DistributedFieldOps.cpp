// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/DistributedFieldOps.hpp"

#include "ovk/core/Comm.hpp"
#include "ovk/core/DistributedField.hpp"
#include "ovk/core/Global.hpp"
#include "ovk/core/Range.hpp"

#include <mpi.h>

namespace ovk {
namespace core {

long long CountDistributedMask(const distributed_field<bool> &Mask) {

  long long Count = 0;

  for (int k = Mask.LocalRange().Begin(2); k < Mask.LocalRange().End(2); ++k) {
    for (int j = Mask.LocalRange().Begin(1); j < Mask.LocalRange().End(1); ++j) {
      for (int i = Mask.LocalRange().Begin(0); i < Mask.LocalRange().End(0); ++i) {
        Count += (long long)(Mask(i,j,k));
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &Count, 1, MPI_LONG_LONG, MPI_SUM, Mask.Comm());

  return Count;

}

}}
