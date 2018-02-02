# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(GLOBAL_WARNING_FLAG FALSE CACHE INTERNAL "")

#====================
# CMake requirements
#====================

cmake_minimum_required(VERSION 2.8.12)

include(CheckIncludeFile)
include(CheckSymbolExists)

#===============
# Project setup
#===============

project(Overkit C CXX)
set(Overkit_VERSION_MAJOR 0)
set(Overkit_VERSION_MINOR 1)
set(Overkit_VERSION_PATCH 0)
set(Overkit_VERSION ${Overkit_VERSION_MAJOR}.${Overkit_VERSION_MINOR}.${Overkit_VERSION_PATCH})

# Module stuff
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/config/modules")
set(C_ENABLED TRUE)
set(CXX_ENABLED TRUE)
set(Fortran_ENABLED FALSE)

# OS X RPATH for dynamic libraries
set(CMAKE_MACOSX_RPATH TRUE)

# Add install lib dir to any executables' RPATH, as well as any paths to external libraries
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Forbid in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed. Please run CMake from a separate directory (you may need to remove CMakeFiles and CMakeCache.txt).")
endif()

# Adds "make distclean" command to delete contents of build directory (safer than rm -rf *)
add_custom_target(distclean COMMAND "${CMAKE_COMMAND}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
  -P "${CMAKE_SOURCE_DIR}/config/scripts/distclean.cmake")

#=========================
# Build types and options
#=========================

# Build types
include(DebugBuildType)
include(ReleaseBuildType)
include(CoverageBuildType)
if(NOT CMAKE_CONFIGURATION_TYPES)
  # Single-configuration generator
  if(CMAKE_BUILD_TYPE)
    set(BUILD_TYPE_VALUE ${CMAKE_BUILD_TYPE})
  else()
    message(STATUS "No build type specified; setting to 'Release'.")
    set(BUILD_TYPE_VALUE Release)
  endif()
  set(CMAKE_BUILD_TYPE ${BUILD_TYPE_VALUE} CACHE STRING
    "Choose the type of build, options are: None Debug Release Coverage." FORCE)
else()
  # Multi-configuration generator
  # (Not sure what to do here)
endif()

# Build options
option(TESTS "Enable tests" OFF)
option(EXAMPLES "Enable examples" OFF)
option(BUILD_SHARED_LIBS "" OFF)

if(NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS FALSE)
endif()

#=======================
# External dependencies
#=======================

if(NOT BUILT_IN_MPI AND NOT EXTERNAL_MPI)
  set(BUILT_IN_MPI FALSE CACHE INTERNAL "")
  set(EXTERNAL_MPI FALSE CACHE INTERNAL "")
  # If the user specified MPI wrappers as compilers, use them
  check_include_file("mpi.h" MPI_HEADER_FOUND)
  if(MPI_HEADER_FOUND)
    check_symbol_exists(MPI_Init "mpi.h" MPI_WORKS)
    if(MPI_WORKS)
      set(BUILT_IN_MPI TRUE CACHE INTERNAL "")
      set(MPI_INCLUDE_PATH "" CACHE INTERNAL "")
      set(MPI_C_LIBRARIES "" CACHE INTERNAL "")
    endif()
  endif()
  # Otherwise, try to find an MPI installation
  if(NOT BUILT_IN_MPI)
    message(STATUS "Not using MPI compiler wrappers; searching for external MPI installation")
    find_package(MPI REQUIRED)
    set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
    check_include_file(mpi.h FOUND_MPI_HEADER)
    if(FOUND_MPI_HEADER)
      set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
      check_symbol_exists(MPI_Init "mpi.h" MPI_WORKS)
      if(MPI_WORKS)
        set(EXTERNAL_MPI TRUE CACHE INTERNAL "")
      endif()
    endif()
    if(NOT EXTERNAL_MPI)
      message(FATAL_ERROR "Unable to build with detected MPI installation.")
    endif()
  endif()
elseif(EXTERNAL_MPI)
  find_package(MPI REQUIRED)
endif()

include_directories(SYSTEM ${MPI_INCLUDE_PATH})
link_libraries(${MPI_C_LIBRARIES})

# Check if MPI_Ibarrier is supported
set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
check_symbol_exists(MPI_Ibarrier "mpi.h" IBARRIER_SUPPORTED)
if(IBARRIER_SUPPORTED)
  add_definitions(-DHAVE_IBARRIER)
else()
  message(WARNING "MPI non-blocking collectives aren't supported. Performance may be degraded for certain operations.")
  set(GLOBAL_WARNING_FLAG TRUE CACHE INTERNAL "")
endif()

#======================
# Extra compiler flags
#======================

set(EXTRA_C_FLAGS "")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${EXTRA_C_FLAGS}" CACHE STRING
  "Flags used by the C compiler during debug builds." FORCE)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_C_FLAGS}" CACHE STRING
  "Flags used by the C compiler during release builds." FORCE)
set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_C_FLAGS_COVERAGE} ${EXTRA_C_FLAGS}" CACHE STRING
  "Flags used by the C compiler during coverage builds." FORCE)

#======================
# Testing and coverage
#======================

# Enable testing
if(TESTS)
  enable_testing()
endif()

# Verbose version of 'make test'
add_custom_target(test-v COMMAND ${CMAKE_CTEST_COMMAND} --verbose)

# Coverage report
add_custom_target(coverage
  COMMAND "${CMAKE_COMMAND}" -DSOURCE_DIR="${CMAKE_SOURCE_DIR}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -DCOVERAGE_DIR="src/core" -P "${CMAKE_SOURCE_DIR}/config/scripts/gen-coverage-data.cmake"
  COMMAND "${CMAKE_COMMAND}" -DSOURCE_DIR="${CMAKE_SOURCE_DIR}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -DCOVERAGE_DIR="src/support" -P "${CMAKE_SOURCE_DIR}/config/scripts/gen-coverage-data.cmake"
  COMMAND "${CMAKE_COMMAND}" -DSOURCE_DIR="${CMAKE_SOURCE_DIR}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -DCOVERAGE_DIR="src/utils" -P "${CMAKE_SOURCE_DIR}/config/scripts/gen-coverage-data.cmake"
  COMMAND "${CMAKE_COMMAND}" -DSOURCE_DIR="${CMAKE_SOURCE_DIR}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -P "${CMAKE_SOURCE_DIR}/config/scripts/gen-coverage-report.cmake"
)

# Reset coverage data
add_custom_target(coverageclean
  COMMAND "${CMAKE_COMMAND}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -P "${CMAKE_SOURCE_DIR}/config/scripts/coverageclean.cmake"
)

#=====================
# Other configuration
#=====================

# Shut ranlib up about no symbols
if(APPLE)
  set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
endif()

#================
# Subdirectories
#================

add_subdirectory(examples)
add_subdirectory(src)
add_subdirectory(tests)

#==============
# Installation
#==============

include(CMakePackageConfigHelpers)

write_basic_package_version_file("${CMAKE_BINARY_DIR}/OverkitConfigVersion.cmake"
  VERSION ${Overkit_VERSION}
  COMPATIBILITY SameMajorVersion
)

export(EXPORT OverkitTargets
  FILE "${CMAKE_BINARY_DIR}/OverkitTargets.cmake"
)

configure_file(config/OverkitConfig.cmake
  "${CMAKE_BINARY_DIR}/OverkitConfig.cmake"
  COPYONLY
)

install(
  FILES
    "${CMAKE_BINARY_DIR}/OverkitConfig.cmake"
    "${CMAKE_BINARY_DIR}/OverkitConfigVersion.cmake"
  DESTINATION lib/Overkit/CMake
)

install(EXPORT OverkitTargets
  DESTINATION lib/Overkit/CMake
)

#=============================
# Write configuration summary
#=============================

message(STATUS "")
message(STATUS "=====================")
message(STATUS "Configuration Summary")
message(STATUS "=====================")
message(STATUS "")
message(STATUS "Build type:         ${CMAKE_BUILD_TYPE}")
message(STATUS "C compiler:         ${CMAKE_C_COMPILER}")
message(STATUS "C++ compiler:       ${CMAKE_CXX_COMPILER}")
message(STATUS "Shared libraries:   ${BUILD_SHARED_LIBS}")
message(STATUS "Install prefix:     ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Examples:           ${EXAMPLES}")
message(STATUS "Tests:              ${TESTS}")
message(STATUS "")
if(GLOBAL_WARNING_FLAG)
  message(STATUS "**** Warnings were produced. See details above. **** ")
  message(STATUS "")
endif()
