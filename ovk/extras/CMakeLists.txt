# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#=================
# Overkit library
#=================

set(OVERKIT_EXTRAS_SOURCES
  Global.c
)

add_library(overkit-extras ${OVERKIT_EXTRAS_SOURCES})

# Link to Overkit
target_link_libraries(overkit-extras overkit)

# Link C math library
target_link_libraries(overkit-extras m)

#================
# Compiler flags
#================

# Nothing here at the moment

#==============
# Installation
#==============

# Version information
set_property(TARGET overkit-extras PROPERTY VERSION ${Overkit_VERSION})
set_property(TARGET overkit-extras PROPERTY SOVERSION ${Overkit_VERSION_MAJOR})
set_property(TARGET overkit-extras PROPERTY INTERFACE_Overkit_VERSION_MAJOR ${Overkit_VERSION_MAJOR})
set_property(TARGET overkit-extras APPEND PROPERTY COMPATIBLE_INTERFACE_STRING Overkit_VERSION_MAJOR)

# Install library
install(TARGETS overkit-extras
  EXPORT OverkitTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include
)

# Install headers
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/" DESTINATION include/ovk/extras
  FILES_MATCHING PATTERN "ovk*.h" PATTERN "ovk*.inl")
