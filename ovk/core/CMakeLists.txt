# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#=================
# Overkit library
#=================

set(OVERKIT_SOURCES
  Connectivity.c
  ConnectivityD.c
  ConnectivityR.c
  Context.c
  Domain.c
  ErrorHandler.c
  Exchange.c
  Global.c
  Grid.c
  Logger.c
  MiscUtils.cpp
  MPIUtils.c
  PartitionHash.c
  ProfileUtils.c
)

add_library(overkit ${OVERKIT_SOURCES})

# Link C math library
target_link_libraries(overkit m)

#================
# Compiler flags
#================

# Define OVERKIT_DEBUG flag in debug & coverage build types
set_property(TARGET overkit APPEND PROPERTY COMPILE_DEFINITIONS
  $<$<OR:$<CONFIG:Debug>,$<CONFIG:Coverage>>:OVERKIT_DEBUG>
)

#==============
# Installation
#==============

# Version information
set_property(TARGET overkit PROPERTY VERSION ${Overkit_VERSION})
set_property(TARGET overkit PROPERTY SOVERSION ${Overkit_VERSION_MAJOR})
set_property(TARGET overkit PROPERTY INTERFACE_Overkit_VERSION_MAJOR ${Overkit_VERSION_MAJOR})
set_property(TARGET overkit APPEND PROPERTY COMPATIBLE_INTERFACE_STRING Overkit_VERSION_MAJOR)

# Install library
install(TARGETS overkit
  EXPORT OverkitTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include
)

# Install headers
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/" DESTINATION include/ovk/core
  FILES_MATCHING PATTERN "ovk*.h" PATTERN "ovk*.inl")
