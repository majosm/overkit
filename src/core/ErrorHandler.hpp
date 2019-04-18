// Copyright (c) 2018 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_CORE_ERROR_HANDLER_HPP_INCLUDED
#define OVK_CORE_ERROR_HANDLER_HPP_INCLUDED

#include <ovk/core/Constants.hpp>
#include <ovk/core/Global.hpp>

#include <mpi.h>

namespace ovk {
namespace core {

class error_handler {

public:

  // Get rid of this when context no longer needs it
  error_handler();
  explicit error_handler(error_handler_type Type);

  error_handler_type Type() const { return Type_; }

  void SetType(error_handler_type Type);

private:

  error_handler_type Type_;

};

#define OVK_EH_CHECK_SKIP_TO(Handler, Error, Label) \
  do { \
    if ((Error) != error::NONE) { \
      goto Label; \
    } \
  } while (0)

#define OVK_EH_CHECK(Handler, Error) \
  do { \
    ovk::error OVK_EH_Error = (Error); \
    if (OVK_EH_Error != error::NONE) { \
      ovk::core::error_handler &OVK_EH_Handler = (Handler); \
      switch (OVK_EH_Handler.Type()) { \
      case ovk::error_handler_type::ABORT: \
        MPI_Abort(MPI_COMM_WORLD, int(OVK_EH_Error)); \
      case ovk::error_handler_type::RETURN: \
        return OVK_EH_Error; \
      } \
    } \
  } while (0)

#define OVK_EH_SYNC(Handler, Error, Comm) \
  do { \
    ovk::error &OVK_EH_Error = (Error); \
    int OVK_EH_MinErrorInt = int(OVK_EH_Error != ovk::error::NONE ? \
      OVK_EH_Error : ovk::error::MAX); \
    MPI_Allreduce(MPI_IN_PLACE, &OVK_EH_MinErrorInt, 1, MPI_INT, MPI_MIN, (Comm)); \
    ovk::error OVK_EH_MinError = ovk::error(OVK_EH_MinErrorInt); \
    if (OVK_EH_MinError != ovk::error::MAX) { \
      OVK_EH_Error = OVK_EH_MinError; \
    } \
  } while (0)

}}

#endif
