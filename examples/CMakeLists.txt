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
list(APPEND C_TARGETS Test)

# add_executable(Debug Debug.c)
# list(APPEND LOCAL_TARGETS Debug)
# list(APPEND C_TARGETS Debug)

if(XPACC)
  add_executable(ReadInterp ReadInterp.c)
  list(APPEND LOCAL_TARGETS ReadInterp)
  list(APPEND C_TARGETS ReadInterp)
endif()

#-------------------
# Compiling/linking
#-------------------

foreach(EXAMPLE ${C_TARGETS})

  set(BASE_C_FLAGS_DEBUG ${WARNING_C_FLAGS})
  set(BASE_C_FLAGS_RELEASE ${OPT_C_FLAGS})
  target_compile_options(${EXAMPLE} PRIVATE
    $<$<CONFIG:SlowDebug>:${BASE_C_FLAGS_DEBUG}>
    $<$<CONFIG:FastDebug>:${BASE_C_FLAGS_DEBUG}>
    $<$<CONFIG:Release>:${BASE_C_FLAGS_RELEASE}>
    $<$<CONFIG:RelWithDebInfo>:${BASE_C_FLAGS_RELEASE}>
    $<$<CONFIG:MinSizeRel>:${BASE_C_FLAGS_RELEASE}>
  )

  # Language feature requirements
  if(BUILT_IN_DIALECT_SUPPORT)
    if(DIALECT_COMPILE_FEATURE_SUPPORT)
      target_compile_features(${EXAMPLE} PRIVATE c_std_99)
    else()
      set_property(TARGET ${EXAMPLE} PROPERTY C_STANDARD 99)
    endif()
  else()
    target_compile_options(${EXAMPLE} PRIVATE ${DIALECT_C_FLAGS})
  endif()

endforeach()

foreach(EXAMPLE ${CXX_TARGETS})

  set(BASE_CXX_FLAGS_DEBUG ${WARNING_CXX_FLAGS})
  set(BASE_CXX_FLAGS_RELEASE ${OPT_CXX_FLAGS})
  target_compile_options(${EXAMPLE} PRIVATE
    $<$<CONFIG:SlowDebug>:${BASE_CXX_FLAGS_DEBUG}>
    $<$<CONFIG:FastDebug>:${BASE_CXX_FLAGS_DEBUG}>
    $<$<CONFIG:Release>:${BASE_CXX_FLAGS_RELEASE}>
    $<$<CONFIG:RelWithDebInfo>:${BASE_CXX_FLAGS_RELEASE}>
    $<$<CONFIG:MinSizeRel>:${BASE_CXX_FLAGS_RELEASE}>
  )

  # Language feature requirements
  if(BUILT_IN_DIALECT_SUPPORT)
    if(DIALECT_COMPILE_FEATURE_SUPPORT)
      target_compile_features(${EXAMPLE} PRIVATE cxx_std_11)
    else()
      set_property(TARGET ${EXAMPLE} PROPERTY CXX_STANDARD 11)
    endif()
  else()
    target_compile_options(${EXAMPLE} PRIVATE ${DIALECT_CXX_FLAGS})
  endif()

endforeach()

foreach(EXAMPLE ${C_TARGETS} ${CXX_TARGETS})

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
