# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

if(COVERAGE)

  # Only available in debug mode
  if(NOT DEBUG)
    message(FATAL_ERROR "Must build in debug mode to enable coverage analysis.")
  endif()

endif()
