# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#=========================
# Directory-wide settings
#=========================

#-------------------
# Compiling/linking
#-------------------

# MPI
if(EXTERNAL_MPI)
  include_directories(SYSTEM "${MPI_INCLUDE_PATH}")
  link_libraries(${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif()

# C math library
link_libraries("${C_MATH_LIBRARY}")

# Overkit library
link_libraries(overkit)

#==========
# Examples
#==========

#------
# Test
#------

add_executable(Test Test.c)
list(APPEND LOCAL_TARGETS Test)

#-------
# Debug
#-------

# add_executable(Debug Debug.c)
# list(APPEND LOCAL_TARGETS Debug)

#------------
# ReadInterp
#------------

if(XPACC)
  add_executable(ReadInterp ReadInterp.c)
  list(APPEND LOCAL_TARGETS ReadInterp)
endif()

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
