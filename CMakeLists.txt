# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

# Initialize warning flag for reporting warnings at end of configuration
set(GLOBAL_WARNING_FLAG FALSE CACHE INTERNAL "")

# Forbid in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed. Please run CMake from a separate directory (you may need to remove CMakeFiles and CMakeCache.txt).")
endif()

#--------------------
# CMake requirements
#--------------------

cmake_minimum_required(VERSION 2.8.12)

include(CheckIncludeFile)
include(CheckSymbolExists)
include(CMakePackageConfigHelpers)

#--------------------
# Project definition
#--------------------

project(OverkitParallel C CXX)
set(OverkitParallel_VERSION_MAJOR 0)
set(OverkitParallel_VERSION_MINOR 1)
set(OverkitParallel_VERSION_PATCH 0)
set(OverkitParallel_VERSION ${OverkitParallel_VERSION_MAJOR}.${OverkitParallel_VERSION_MINOR}.${OverkitParallel_VERSION_PATCH})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/config/modules")
set(C_ENABLED TRUE)
set(CXX_ENABLED TRUE)
set(Fortran_ENABLED FALSE)

#=========================
# Build types and options
#=========================

# Build types
include(DebugBuildType)
include(ReleaseBuildType)
if(NOT CMAKE_CONFIGURATION_TYPES)
  # Single-configuration generator
  if(CMAKE_BUILD_TYPE)
    set(BUILD_TYPE_VALUE ${CMAKE_BUILD_TYPE})
  else()
    message(STATUS "No build type specified; setting to 'Release'.")
    set(BUILD_TYPE_VALUE Release)
  endif()
  set(CMAKE_BUILD_TYPE ${BUILD_TYPE_VALUE} CACHE STRING
    "Choose the type of build, options are: None Debug Release." FORCE)
endif()

# Options
option(BUILD_SHARED_LIBS "Build libraries as shared" OFF)
option(EXAMPLES "Enable examples" OFF)
option(TESTS "Enable tests" OFF)
option(COVERAGE "Enable coverage analysis" OFF)
option(PROFILE "Enable profiling" OFF)
option(TIMERS "Enable internal timers" OFF)
option(XPACC "Enable XPACC-specific extras" OFF)

if(NOT DEFINED SUBPROJECT)
  set(SUBPROJECT FALSE CACHE INTERNAL "")
endif()

#=======================
# External dependencies
#=======================

#-----
# MPI
#-----

if(NOT COMPILER_BUILT_IN_MPI AND NOT EXTERNAL_MPI)
  set(COMPILER_BUILT_IN_MPI FALSE CACHE INTERNAL "")
  set(EXTERNAL_MPI FALSE CACHE INTERNAL "")
  # If the user specified MPI wrappers as compilers, use them
  check_include_file("mpi.h" MPI_HEADER_FOUND)
  if(MPI_HEADER_FOUND)
    check_symbol_exists(MPI_Init "mpi.h" MPI_WORKS)
    if(MPI_WORKS)
      set(COMPILER_BUILT_IN_MPI TRUE CACHE INTERNAL "")
      set(MPI_INCLUDE_PATH "" CACHE INTERNAL "")
      set(MPI_C_LIBRARIES "" CACHE INTERNAL "")
    endif()
  endif()
  # Otherwise, try to find an MPI installation
  if(NOT COMPILER_BUILT_IN_MPI)
    message(STATUS "Not using MPI compiler wrappers; searching for external MPI installation")
    find_package(MPI REQUIRED)
    set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
    check_include_file(mpi.h FOUND_MPI_HEADER)
    if(FOUND_MPI_HEADER)
      # TODO: figure out why this doesn't work on some machines (e.g. Ray)
#       set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
#       check_symbol_exists(MPI_Init "mpi.h" MPI_WORKS)
#       if(MPI_WORKS)
        set(EXTERNAL_MPI TRUE CACHE INTERNAL "")
#       endif()
    endif()
    if(NOT EXTERNAL_MPI)
      message(FATAL_ERROR "Unable to build with detected MPI installation.")
    endif()
  endif()
elseif(EXTERNAL_MPI)
  find_package(MPI REQUIRED)
endif()

# Check if MPI_Ibarrier is supported
set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
check_symbol_exists(MPI_Ibarrier "mpi.h" HAVE_MPI_IBARRIER)
if(NOT HAVE_MPI_IBARRIER)
  message(WARNING "MPI non-blocking collectives aren't supported. Performance may be degraded for certain operations.")
  set(GLOBAL_WARNING_FLAG TRUE CACHE INTERNAL "")
endif()

#----------------
# C math library
#----------------

find_library(C_MATH_LIBRARY m)
if(NOT C_MATH_LIBRARY)
  message(FATAL_ERROR "Could not find C math library.")
endif()

#=======================
# Project configuration
#=======================

#---------
# Testing
#---------

if(TESTS)

  enable_testing()

  # Verbose version of 'make test'
  add_custom_target(test-v COMMAND ${CMAKE_CTEST_COMMAND} --verbose)

endif()

#----------
# Coverage
#----------

if(COVERAGE)

  # Requires GNU
  if (NOT "${CMAKE_C_COMPILER_ID}" MATCHES "GNU")
    message(FATAL_ERROR "Building with coverage analysis requires GNU C compiler.")
  endif()
  if (NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    message(FATAL_ERROR "Building with coverage analysis requires GNU C++ compiler.")
  endif()

  # Coverage report
  add_custom_target(coverage
    COMMAND "${CMAKE_COMMAND}" -DSOURCE_DIR="${CMAKE_SOURCE_DIR}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
      -DCOVERAGE_DIRS="ovk/core;ovk/extras" -P "${CMAKE_SOURCE_DIR}/config/scripts/gen-coverage-report.cmake"
  )

  # Reset coverage data
  add_custom_target(coverage-clean
    COMMAND "${CMAKE_COMMAND}" -DBINARY_DIR="${CMAKE_BINARY_DIR}" -P
      "${CMAKE_SOURCE_DIR}/config/scripts/coverage-clean.cmake"
  )

endif()

#-------
# RPATH
#-------

if(BUILD_SHARED_LIBS)

  # Enable RPATH on Mac OS X
  if(APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)
  endif()

  # Add install lib dir to any executables' RPATH, as well as any paths to external libraries
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

endif()

#-------
# Misc.
#-------

# Distclean
add_custom_target(distclean COMMAND "${CMAKE_COMMAND}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
  -P "${CMAKE_SOURCE_DIR}/config/scripts/distclean.cmake")

# Shut ranlib up about no symbols
if(APPLE)
  set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
endif()

# API headers end up in the include directory of the build tree
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/include")
include_directories("${CMAKE_BINARY_DIR}/include")

# Internal headers end up in the include-internal directory of the build tree
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/include-internal")
include_directories("${CMAKE_BINARY_DIR}/include-internal")

# set(EXTRA_C_FLAGS "")
# set(EXTRA_CXX_FLAGS "")
# set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${EXTRA_C_FLAGS}" CACHE STRING
#   "Flags used by the C compiler during debug builds." FORCE)
# set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_C_FLAGS}" CACHE STRING
#   "Flags used by the C compiler during release builds." FORCE)
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${EXTRA_CXX_FLAGS}" CACHE STRING
#   "Flags used by the C++ compiler during debug builds." FORCE)
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_CXX_FLAGS}" CACHE STRING
#   "Flags used by the C++ compiler during release builds." FORCE)

#==================
# Pre-build checks
#==================

# Extra build configuration checks that can't be done at configure time
add_custom_target(config-check
  COMMAND "${CMAKE_COMMAND}"
    -DDEBUG=$<CONFIG:Debug>
    -DCOVERAGE=${COVERAGE}
    -P "${CMAKE_SOURCE_DIR}/config/scripts/check-config.cmake"
)

add_custom_target(pre-build ALL)
add_dependencies(pre-build config-check)

#================
# Subdirectories
#================

# Main library code
add_subdirectory(ovk)

# Supporting code for tests, examples, utils
add_subdirectory(support)

# Utility programs
add_subdirectory(utils)

# Test code
if(TESTS)
  add_subdirectory(tests)
endif()

# Example programs
if(EXAMPLES)
  add_subdirectory(examples)
endif()

#==============================
# Package configuration/export
#==============================

#-------------------
# Package variables
#-------------------

# CMake
if(SUBPROJECT)
  set(OverkitParallel_CMAKE_DIR "${CMAKE_BINARY_DIR}")
  set(OverkitParallel_CMAKE_INCLUDES "${CMAKE_BINARY_DIR}/include")
else()
  set(OverkitParallel_CMAKE_DIR "${CMAKE_INSTALL_PREFIX}")
  set(OverkitParallel_CMAKE_INCLUDES "${CMAKE_INSTALL_PREFIX}/include")
endif()
set(OverkitParallel_CMAKE_LIBRARIES overkit-parallel)
set(OverkitParallel_CMAKE_HEADERS overkit.h)

set(OverkitParallel_CMAKE_OPTION_SHARED_LIBS ${BUILD_SHARED_LIBS})
set(OverkitParallel_CMAKE_OPTION_EXAMPLES ${EXAMPLES})
set(OverkitParallel_CMAKE_OPTION_TESTS ${TESTS})
set(OverkitParallel_CMAKE_OPTION_COVERAGE ${COVERAGE})
set(OverkitParallel_CMAKE_OPTION_PROFILE ${PROFILE})
set(OverkitParallel_CMAKE_OPTION_TIMERS ${TIMERS})
set(OverkitParallel_CMAKE_OPTION_XPACC ${XPACC})

configure_file("${CMAKE_SOURCE_DIR}/config/OverkitDefs.cmake.in"
  "${CMAKE_BINARY_DIR}/OverkitParallelDefs.cmake" @ONLY)

install(FILES "${CMAKE_BINARY_DIR}/OverkitParallelDefs.cmake" DESTINATION lib/OverkitParallel/CMake)

# Make
if(SUBPROJECT)
  set(OverkitParallel_MAKE_DIR "${CMAKE_BINARY_DIR}")
  set(OverkitParallel_MAKE_CPPFLAGS "-I\"${CMAKE_BINARY_DIR}/include\"")
  set(OverkitParallel_MAKE_LDFLAGS "-L\"${CMAKE_BINARY_DIR}/lib\"")
else()
  set(OverkitParallel_MAKE_DIR "${CMAKE_INSTALL_PREFIX}")
  set(OverkitParallel_MAKE_CPPFLAGS "-I\"${CMAKE_INSTALL_PREFIX}/include\"")
  set(OverkitParallel_MAKE_LDFLAGS "-L\"${CMAKE_INSTALL_PREFIX}/lib\"")
endif()
set(OverkitParallel_MAKE_LIBS "-loverkit-parallel")
set(OverkitParallel_MAKE_HEADERS overkit.h)

macro(convert_to_01 IN OUT)
if(${${IN}})
  set(${OUT} 1)
else()
  set(${OUT} 0)
endif()
endmacro()

convert_to_01(BUILD_SHARED_LIBS OverkitParallel_MAKE_OPTION_SHARED_LIBS)
convert_to_01(EXAMPLES OverkitParallel_MAKE_OPTION_EXAMPLES)
convert_to_01(TESTS OverkitParallel_MAKE_OPTION_TESTS)
convert_to_01(COVERAGE OverkitParallel_MAKE_OPTION_COVERAGE)
convert_to_01(PROFILE OverkitParallel_MAKE_OPTION_PROFILE)
convert_to_01(TIMERS OverkitParallel_MAKE_OPTION_TIMERS)
convert_to_01(XPACC OverkitParallel_MAKE_OPTION_XPACC)

configure_file("${CMAKE_SOURCE_DIR}/config/OverkitDefs.mk.in"
  "${CMAKE_BINARY_DIR}/OverkitParallelDefs.mk" @ONLY)

install(FILES "${CMAKE_BINARY_DIR}/OverkitParallelDefs.mk" DESTINATION lib/OverkitParallel/Make)

#--------------------------
# Exports and version info
#--------------------------

configure_file("${CMAKE_SOURCE_DIR}/config/OverkitConfig.cmake"
  "${CMAKE_BINARY_DIR}/OverkitParallelConfig.cmake" COPYONLY)

write_basic_package_version_file("${CMAKE_BINARY_DIR}/OverkitParallelConfigVersion.cmake"
  VERSION ${OverkitParallel_VERSION}
  COMPATIBILITY SameMajorVersion
)

if(CMAKE_MAJOR_VERSION GREATER 3)
  export(EXPORT OverkitParallelTargets FILE "${CMAKE_BINARY_DIR}/OverkitParallelTargets.cmake")
else()
  export(TARGETS overkit-parallel FILE "${CMAKE_BINARY_DIR}/OverkitParallelTargets.cmake")
endif()

install(
  FILES
    "${CMAKE_BINARY_DIR}/OverkitParallelConfig.cmake"
    "${CMAKE_BINARY_DIR}/OverkitParallelConfigVersion.cmake"
  DESTINATION lib/OverkitParallel/CMake
)

install(EXPORT OverkitParallelTargets DESTINATION lib/OverkitParallel/CMake)

#=============================
# Write configuration summary
#=============================

message(STATUS "")
message(STATUS "=====================")
message(STATUS "Configuration Summary")
message(STATUS "=====================")
message(STATUS "")
if(CMAKE_CONFIGURATION_TYPES)
message(STATUS "Configuration types: ${CMAKE_CONFIGURATION_TYPES}")
else()
message(STATUS "Build type:          ${CMAKE_BUILD_TYPE}")
endif()
message(STATUS "C compiler:          ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler:        ${CMAKE_CXX_COMPILER}")
message(STATUS "Shared libraries:    ${BUILD_SHARED_LIBS}")
message(STATUS "Install prefix:      ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Examples:            ${EXAMPLES}")
message(STATUS "Tests:               ${TESTS}")
message(STATUS "Coverage:            ${COVERAGE}")
message(STATUS "Profiling:           ${PROFILE}")
message(STATUS "Timers:              ${TIMERS}")
if(XPACC)
message(STATUS "XPACC:               ${XPACC}")
endif()
message(STATUS "")
if(GLOBAL_WARNING_FLAG)
  message(STATUS "**** Warnings were produced. See details above. **** ")
  message(STATUS "")
endif()
