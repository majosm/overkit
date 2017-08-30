# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

if("${BINARY_DIR}" STREQUAL "")
  message(FATAL_ERROR "Incorrectly set BINARY_DIR.")
endif()

# Remove coverage data files
file(GLOB_RECURSE GCDA_FILES "${BINARY_DIR}/*.gcda")
file(REMOVE_RECURSE ${GCDA_FILES})

# Remove coverage report directory
file(REMOVE_RECURSE "${BINARY_DIR}/coverage")
