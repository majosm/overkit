# Copyright (c) 2020 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

# Warn if no HDF5
if(NOT HAVE_HDF5)
  message(WARNING "Not building with HDF5; will not be able to write XDMF output files in examples.")
  set(GLOBAL_WARNING_FLAG TRUE CACHE INTERNAL "")
endif()

#==========
# Examples
#==========

#----------------
# Header targets
#----------------

# C++ header
set(BUILT_COMMON_CXX_HEADER ${CMAKE_BINARY_DIR}/include/examples/Common.hpp)
set(COMMON_CXX_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/Common.hpp)
add_custom_command(OUTPUT ${BUILT_COMMON_CXX_HEADER}
  MAIN_DEPENDENCY ${COMMON_CXX_HEADER}
  COMMAND "${CMAKE_COMMAND}" -E copy "${COMMON_CXX_HEADER}" "${BUILT_COMMON_CXX_HEADER}"
)
set_source_files_properties(${BUILT_COMMON_CXX_HEADER} PROPERTIES GENERATED TRUE)
add_custom_target(examples-common-cxx-header ALL DEPENDS ${BUILT_COMMON_CXX_HEADER})
list(APPEND LOCAL_TARGETS examples-common-cxx-header)

# C header
set(BUILT_COMMON_C_HEADER ${CMAKE_BINARY_DIR}/include/examples/Common.h)
set(COMMON_C_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/Common.h)
add_custom_command(OUTPUT ${BUILT_COMMON_C_HEADER}
  MAIN_DEPENDENCY ${COMMON_C_HEADER}
  COMMAND "${CMAKE_COMMAND}" -E copy "${COMMON_C_HEADER}" "${BUILT_COMMON_C_HEADER}"
)
set_source_files_properties(${BUILT_COMMON_C_HEADER} PROPERTIES GENERATED TRUE)
add_custom_target(examples-common-c-header ALL DEPENDS ${BUILT_COMMON_C_HEADER})
list(APPEND LOCAL_TARGETS examples-common-c-header)

#-----------------
# Example targets
#-----------------

add_executable(Interface Interface.cpp)
list(APPEND LOCAL_TARGETS Interface)
list(APPEND CXX_TARGETS Interface)

add_executable(InterfaceC InterfaceC.c)
list(APPEND LOCAL_TARGETS InterfaceC)
list(APPEND C_TARGETS InterfaceC)

add_executable(Blobs Blobs.cpp)
list(APPEND LOCAL_TARGETS Blobs)
list(APPEND CXX_TARGETS Blobs)

add_executable(BoxInBox BoxInBox.cpp)
list(APPEND LOCAL_TARGETS BoxInBox)
list(APPEND CXX_TARGETS BoxInBox)

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
    target_link_libraries(${EXAMPLE} PRIVATE ${PROFILE_EXE_LINKER_FLAGS})
  endif()

  # MPI
  if(EXTERNAL_MPI)
    target_include_directories(${EXAMPLE} SYSTEM PRIVATE ${MPI_INCLUDES})
    target_link_libraries(${EXAMPLE} PRIVATE ${MPI_LIBS})
  endif()

  # C math library
  target_link_libraries(${EXAMPLE} PRIVATE ${C_MATH_LIBRARY})

  # Overkit
  target_link_libraries(${EXAMPLE} PRIVATE overkit)

  # Support library
  target_link_libraries(${EXAMPLE} PRIVATE support)

endforeach()

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
