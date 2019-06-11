# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#========================
# Core C wrapper library
#========================

#--------------
# Source files
#--------------

set(SOURCES
  AssemblyOptions.cpp
  ConnectivityComponent.cpp
  ConnectivityM.cpp
  ConnectivityN.cpp
  Context.cpp
  Domain.cpp
  Exchanger.cpp
  Grid.cpp
  Request.cpp
)

set(PUBLIC_HEADERS
  AssemblyOptions.h
  ConnectivityComponent.h
  ConnectivityM.h
  ConnectivityN.h
  Constants.h
  Context.h
  DataType.h
  Domain.h
  Error.h
  Exchanger.h
  Global.h
  Grid.h
  Request.h
)

set(INTERNAL_HEADERS
  # Nothing here yet
)

#-----------------
# Library targets
#-----------------

add_library(overkit-core-c OBJECT ${SOURCES})
list(APPEND LOCAL_TARGETS overkit-core-c)

# Build core first
add_dependencies(overkit-core-c overkit-core)

# Generate/copy headers to build tree before compiling
add_dependencies(overkit-core-c core-c-headers)

#-------------------
# Compiling/linking
#-------------------

# Baseline compiler flags
set(BASE_CXX_FLAGS_DEBUG ${WARNING_CXX_FLAGS})
set(BASE_CXX_FLAGS_RELEASE ${OPT_CXX_FLAGS})
target_compile_options(overkit-core-c PRIVATE
  $<$<CONFIG:SlowDebug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:FastDebug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:Release>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:RelWithDebInfo>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:MinSizeRel>:${BASE_CXX_FLAGS_RELEASE}>
)

# Coverage
if(COVERAGE)
  target_compile_options(overkit-core-c PRIVATE ${COVERAGE_COMPILE_FLAGS})
endif()

# Profiling
if(PROFILE)
  target_compile_options(overkit-core-c PRIVATE ${PROFILE_COMPILE_FLAGS})
endif()

# Language feature requirements
if(BUILT_IN_DIALECT_SUPPORT)
  if(DIALECT_COMPILE_FEATURE_SUPPORT)
    target_compile_features(overkit-core-c PRIVATE cxx_std_11)
  else()
    set_property(TARGET overkit-core-c PROPERTY CXX_STANDARD 11)
  endif()
else()
  target_compile_options(overkit-core-c PRIVATE ${DIALECT_CXX_FLAGS})
endif()

# Needed for linking into shared lib
if(BUILD_SHARED_LIBS)
  set_property(TARGET overkit-core-c PROPERTY POSITION_INDEPENDENT_CODE TRUE)
endif()

# MPI
if(EXTERNAL_MPI)
  target_include_directories(overkit-core-c SYSTEM PUBLIC ${MPI_INCLUDE_PATH})
endif()

#--------------------
# API header targets
#--------------------

set(BUILT_HEADER_PREFIX ovk/core-c)

set(BUILT_HEADERS)
foreach(HEADER_NAME ${PUBLIC_HEADERS} ${INTERNAL_HEADERS})
  set(BUILT_HEADER ${CMAKE_BINARY_DIR}/include/${BUILT_HEADER_PREFIX}/${HEADER_NAME})
  set(HEADER ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER_NAME})
  add_custom_command(OUTPUT ${BUILT_HEADER}
    MAIN_DEPENDENCY ${HEADER}
    COMMAND "${CMAKE_COMMAND}" -E copy "${HEADER}" "${BUILT_HEADER}"
  )
  install(FILES ${BUILT_HEADER} DESTINATION include/${BUILT_HEADER_PREFIX})
  list(APPEND BUILT_HEADERS ${BUILT_HEADER})
endforeach()
set_source_files_properties(${BUILT_HEADERS} PROPERTIES GENERATED TRUE)

# Generated header containing includes for all public headers
set(PREFIXED_HEADERS)
foreach(HEADER_NAME ${PUBLIC_HEADERS})
  list(APPEND PREFIXED_HEADERS ${BUILT_HEADER_PREFIX}/${HEADER_NAME})
endforeach()
set(BUILT_ALL_HEADER ${CMAKE_BINARY_DIR}/include/${BUILT_HEADER_PREFIX}/all.h)
set(ALL_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/all.h.in)
add_custom_command(OUTPUT ${BUILT_ALL_HEADER}
  MAIN_DEPENDENCY ${ALL_HEADER}
  COMMAND "${CMAKE_COMMAND}" -DHEADER_IN="${ALL_HEADER}" -DHEADER_OUT="${BUILT_ALL_HEADER}"
    -DINCLUDE_HEADERS="${PREFIXED_HEADERS}" -DSYSTEM=1 -P
    "${CMAKE_SOURCE_DIR}/config/scripts/gen-aggregate-header.cmake"
)
install(FILES ${BUILT_ALL_HEADER} DESTINATION include/${BUILT_HEADER_PREFIX})
set_source_files_properties(${BUILT_ALL_HEADER} PROPERTIES GENERATED TRUE)

add_custom_target(core-c-headers ALL DEPENDS ${BUILT_HEADERS} ${BUILT_ALL_HEADER})
list(APPEND LOCAL_TARGETS core-c-headers)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
