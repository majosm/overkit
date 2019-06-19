// Copyright (c) 2019 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ERROR_HPP_INCLUDED
#define OVK_CORE_ERROR_HPP_INCLUDED

#include <ovk/core/Comm.hpp>
#include <ovk/core/Error.h>
#include <ovk/core/Global.hpp>

#include <mpi.h>

#include <stdexcept>

namespace ovk {

enum class error {
  NONE = OVK_ERROR_NONE,
  FILE_OPEN = OVK_ERROR_FILE_OPEN,
  FILE_READ = OVK_ERROR_FILE_READ,
  FILE_WRITE = OVK_ERROR_FILE_WRITE
};

namespace error_internal {

inline const char *ExceptionString(error Error) {

  switch (Error) {
  case error::FILE_OPEN:
    return "ovk::error::FILE_OPEN";
  case error::FILE_READ:
    return "ovk::error::FILE_READ";
  case error::FILE_WRITE:
    return "ovk::error::FILE_WRITE";
  default:
    return "";
  }

}

}

class exception : public std::runtime_error {

public:

  explicit exception(error Error):
    runtime_error(error_internal::ExceptionString(Error)),
    Error_(Error)
  {}

  error Error() const noexcept { return Error_; }

private:

  error Error_;

};

class io_error : public exception {
public:
  explicit io_error(error Error):
    exception(Error)
  {}
};

class file_open_error : public io_error {
public:
  file_open_error():
    io_error(error::FILE_OPEN)
  {}
};

class file_read_error : public io_error {
public:
  file_read_error():
    io_error(error::FILE_READ)
  {}
};

class file_write_error : public io_error {
public:
  file_write_error():
    io_error(error::FILE_WRITE)
  {}
};

namespace core {
void SyncError(error &Error, comm_view Comm);
void CheckError(error Error);
void ThrowError(error Error);
}

}

#endif