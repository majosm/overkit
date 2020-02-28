// Copyright (c) 2020 Matthew J. Smith and Overkit contributors
// License: MIT (http://opensource.org/licenses/MIT)

#ifndef OVK_SUPPORT_XDMF_H_LOADED
#define OVK_SUPPORT_XDMF_H_LOADED

#include <ovk/core/Global.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef OVK_HAVE_HDF5

#define OVK_HAVE_XDMF 1

typedef enum {
  SUPPORT_XDMF_ERROR_NONE = 0,
  SUPPORT_XDMF_ERROR_FILE_CREATE,
  SUPPORT_XDMF_ERROR_FILE_OPEN,
  SUPPORT_XDMF_ERROR_FILE_READ,
  SUPPORT_XDMF_ERROR_FILE_WRITE
} support_xdmf_error;

typedef enum {
  SUPPORT_XDMF_ATTRIBUTE_TYPE_INT = 0,
  SUPPORT_XDMF_ATTRIBUTE_TYPE_LONG_LONG,
  SUPPORT_XDMF_ATTRIBUTE_TYPE_DOUBLE
} support_xdmf_attribute_type;

#endif

#ifdef __cplusplus
}
#endif

#endif
