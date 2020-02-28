# Copyright (c) 2020 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

separate_arguments(INCLUDE_HEADERS)

if(${SYSTEM})
  set(INCLUDES)
  foreach(HEADER ${INCLUDE_HEADERS})
    set(INCLUDES "${INCLUDES}#include <${HEADER}>\n")
  endforeach()
else()
  set(INCLUDES)
  foreach(HEADER ${INCLUDE_HEADERS})
    set(INCLUDES "${INCLUDES}#include \"${HEADER}\"\n")
  endforeach()
endif()

configure_file(${HEADER_IN} ${HEADER_OUT})
