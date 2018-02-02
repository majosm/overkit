# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

if(NOT DEFINED C_ENABLED)
  message(FATAL_ERROR "Must set C_ENABLED for ReleaseBuildType")
endif()
if(NOT DEFINED CXX_ENABLED)
  message(FATAL_ERROR "Must set CXX_ENABLED for ReleaseBuildType")
endif()
if(NOT DEFINED Fortran_ENABLED)
  message(FATAL_ERROR "Must set Fortran_ENABLED for ReleaseBuildType")
endif()

if(NOT DEFINED COVERAGEBUILDTYPE)
  if(C_ENABLED)
    set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_C_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage"
      CACHE STRING "Flags used by the C compiler during coverage builds." FORCE)
    mark_as_advanced(CMAKE_C_FLAGS_COVERAGE)
  endif()
  if(CXX_ENABLED)
    set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage"
      CACHE STRING "Flags used by the C++ compiler during coverage builds." FORCE)
    mark_as_advanced(CMAKE_CXX_FLAGS_COVERAGE)
  endif()
  if(Fortran_ENABLED)
    set(CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage"
      CACHE STRING "Flags used by the Fortran compiler during coverage builds." FORCE)
    mark_as_advanced(CMAKE_Fortran_FLAGS_COVERAGE)
  endif()
  set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS_DEBUG}" CACHE STRING
    "Flags used for linking binaries during coverage builds." FORCE)
  mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_COVERAGE)
  set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}" CACHE STRING
    "Flags used by the shared libraries linker during coverage builds." FORCE)
  mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_COVERAGE)
  set(COVERAGEBUILDTYPE TRUE CACHE INTERNAL "")
endif()
