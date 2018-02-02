# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

if("${BINARY_DIR}" STREQUAL "")
  message(FATAL_ERROR "Incorrectly set BINARY_DIR.")
endif()

file(GLOB DISTCLEAN_FILES "${BINARY_DIR}/*")
file(REMOVE_RECURSE ${DISTCLEAN_FILES})
