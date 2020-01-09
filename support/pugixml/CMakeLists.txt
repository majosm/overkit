# Copyright (c) 2019 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#=========================
# pugixml for XML parsing
#=========================

#--------------
# Source files
#--------------

set(HEADERS
  pugiconfig.hpp
  pugixml.cpp
  pugixml.hpp
)

#----------------
# Header targets
#----------------

set(BUILT_HEADER_PREFIX support/pugixml)
set(BUILT_HEADERS)
foreach(HEADER_NAME ${HEADERS})
  set(BUILT_HEADER ${CMAKE_BINARY_DIR}/include/${BUILT_HEADER_PREFIX}/${HEADER_NAME})
  set(HEADER ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER_NAME})
  add_custom_command(OUTPUT ${BUILT_HEADER}
    MAIN_DEPENDENCY ${HEADER}
    COMMAND "${CMAKE_COMMAND}" -E copy "${HEADER}" "${BUILT_HEADER}"
  )
  list(APPEND BUILT_HEADERS ${BUILT_HEADER})
endforeach()
set_source_files_properties(${BUILT_HEADERS} PROPERTIES GENERATED TRUE)

add_custom_target(pugixml-headers ALL DEPENDS ${BUILT_HEADERS})
list(APPEND LOCAL_TARGETS pugixml-headers)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
