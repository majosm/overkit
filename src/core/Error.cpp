// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#include "ovk/core/Error.hpp"

#include "ovk/core/Global.hpp"

#include <limits>
#include <type_traits>

namespace ovk {
namespace core {

void SyncError(error &Error, comm_view Comm) {

  error MaxError = error(std::numeric_limits<typename std::underlying_type<error>::type>::max());

  int MinErrorInt = int(Error != error::NONE ? Error : MaxError);
  MPI_Allreduce(MPI_IN_PLACE, &MinErrorInt, 1, MPI_INT, MPI_MIN, Comm);
  error MinError = error(MinErrorInt);

  if (MinError != MaxError) {
    Error = MinError;
  }

}

void CheckError(error Error) {

  if (Error != error::NONE) {
    ThrowError(Error);
  }

}

void ThrowError(error Error) {

  switch (Error) {
  case error::FILE_OPEN:
    throw file_open_error();
  case error::FILE_READ:
    throw file_read_error();
  case error::FILE_WRITE:
    throw file_write_error();
  default:
    break;
  }

}

}}
