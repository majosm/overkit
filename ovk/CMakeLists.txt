# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#=================
# Overkit library
#=================

#----------------
# Library target
#----------------

add_library(overkit-parallel
  $<TARGET_OBJECTS:overkit-core-c>
  $<TARGET_OBJECTS:overkit-core-cxx>
  $<TARGET_OBJECTS:overkit-extras>
)
list(APPEND LOCAL_TARGETS overkit-parallel)

add_library(OverkitParallel::overkit-parallel ALIAS overkit-parallel)

set_target_properties(overkit-parallel PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

target_include_directories(overkit-parallel PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_property(TARGET overkit-parallel PROPERTY VERSION ${OverkitParallel_VERSION})
set_property(TARGET overkit-parallel PROPERTY SOVERSION ${OverkitParallel_VERSION_MAJOR})
set_property(TARGET overkit-parallel PROPERTY INTERFACE_OverkitParallel_VERSION_MAJOR ${OverkitParallel_VERSION_MAJOR})
set_property(TARGET overkit-parallel APPEND PROPERTY COMPATIBLE_INTERFACE_STRING OverkitParallel_VERSION_MAJOR)

install(TARGETS overkit-parallel
  EXPORT OverkitParallelTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include
)

#-------------------
# Compiling/linking
#-------------------

# MPI
if(EXTERNAL_MPI)
  target_include_directories(overkit-parallel SYSTEM PUBLIC ${MPI_INCLUDE_PATH})
  target_link_libraries(overkit-parallel PUBLIC ${MPI_C_LIBRARIES})
endif()

# C math library
target_link_libraries(overkit-parallel PUBLIC ${C_MATH_LIBRARY})

# Coverage
if(COVERAGE)
  target_link_libraries(overkit-parallel PUBLIC ${COVERAGE_LINK_FLAGS})
endif()

#---------------
# Header target
#---------------

set(BUILT_OVERKIT_HEADER ${CMAKE_BINARY_DIR}/include/overkit.h)
set(OVERKIT_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/overkit.h)

add_custom_command(OUTPUT ${BUILT_OVERKIT_HEADER}
  MAIN_DEPENDENCY ${OVERKIT_HEADER}
  COMMAND "${CMAKE_COMMAND}" -E copy "${OVERKIT_HEADER}" "${BUILT_OVERKIT_HEADER}"
)

install(FILES ${BUILT_OVERKIT_HEADER} DESTINATION include)

set_source_files_properties(${BUILT_OVERKIT_HEADER} PROPERTIES GENERATED TRUE)

add_custom_target(overkit-header ALL DEPENDS ${BUILT_OVERKIT_HEADER})
list(APPEND LOCAL_TARGETS overkit-header)

#================
# Subdirectories
#================

# Core functionality
add_subdirectory(core)

# Extra stuff that may be useful
add_subdirectory(extras)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
