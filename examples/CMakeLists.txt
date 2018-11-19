# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#==========
# Examples
#==========

#-----------------
# Example targets
#-----------------

add_executable(Test Test.c)
list(APPEND LOCAL_TARGETS Test)

# add_executable(Debug Debug.c)
# list(APPEND LOCAL_TARGETS Debug)

if(XPACC)
  add_executable(ReadInterp ReadInterp.c)
  list(APPEND LOCAL_TARGETS ReadInterp)
endif()

#-------------------
# Compiling/linking
#-------------------

foreach(EXAMPLE ${LOCAL_TARGETS})

  # Baseline compiler flags
  set(BASE_C_FLAGS_DEBUG -O0 ${DIALECT_C_FLAGS} ${WARNING_C_FLAGS})
  set(BASE_C_FLAGS_RELEASE ${DIALECT_C_FLAGS} ${OPT_C_FLAGS})
  target_compile_options(${EXAMPLE} PRIVATE
    $<$<CONFIG:Debug>:${BASE_C_FLAGS_DEBUG}>
    $<$<CONFIG:Release>:${BASE_C_FLAGS_RELEASE}>
    $<$<CONFIG:RelWithDebInfo>:${BASE_C_FLAGS_RELEASE}>
    $<$<CONFIG:MinSizeRel>:${BASE_C_FLAGS_RELEASE}>
  )

  # Profiling
  if(PROFILE)
    target_compile_options(${EXAMPLE} PRIVATE ${PROFILE_COMPILE_FLAGS})
  endif()

  # MPI
  if(EXTERNAL_MPI)
    target_include_directories(${EXAMPLE} SYSTEM PRIVATE ${MPI_INCLUDE_PATH})
    target_link_libraries(${EXAMPLE} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
  endif()

  # C math library
  target_link_libraries(${EXAMPLE} ${C_MATH_LIBRARY})

  # Overkit library
  target_link_libraries(${EXAMPLE} overkit-parallel)

endforeach()

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
