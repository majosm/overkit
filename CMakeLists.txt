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

#===============
# Project setup
#===============

project(Overkit C)
set(Overkit_VERSION_MAJOR 0)
set(Overkit_VERSION_MINOR 1)
set(Overkit_VERSION_PATCH 0)
set(Overkit_VERSION ${Overkit_VERSION_MAJOR}.${Overkit_VERSION_MINOR}.${Overkit_VERSION_PATCH})

# Module stuff
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/config/modules")
set(Fortran_ENABLED FALSE)
set(C_ENABLED TRUE)
set(CXX_ENABLED FALSE)

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
option(EXAMPLES "Build examples" OFF)
option (BUILD_SHARED_LIBS "" OFF)

if(NOT DEFINED BUILD_SHARED_LIBS)
  set(BUILD_SHARED_LIBS FALSE)
endif()

#=======================
# External dependencies
#=======================

# MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
link_libraries(${MPI_C_LIBRARIES})

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
message(STATUS "Shared libraries:   ${BUILD_SHARED_LIBS}")
message(STATUS "Install prefix:     ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Examples:           ${EXAMPLES}")
message(STATUS "Tests:              ${TESTS}")
message(STATUS "")
if(GLOBAL_WARNING_FLAG)
  message(STATUS "**** Warnings were produced. See details above. **** ")
  message(STATUS "")
endif()
