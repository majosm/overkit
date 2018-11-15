# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#================
# Extras library
#================

#--------------
# Source files
#--------------

set(SOURCES
  Global.c
)
if(XPACC)
  list(APPEND SOURCES
    XINTOUT.c
  )
endif()

set(INTERNAL_HEADERS
  Global.h
)
if(XPACC)
  list(APPEND INTERNAL_HEADERS
    XINTOUT.h
  )
endif()

set(API_HEADERS
  ovkGlobal.h
)
if(XPACC)
  list(APPEND API_HEADERS
    ovkXINTOUT.h
  )
endif()

#----------------
# Library target
#----------------

add_library(overkit-extras OBJECT ${SOURCES})
list(APPEND LOCAL_TARGETS overkit-extras)

# Build core stuff first
add_dependencies(overkit-extras overkit-core-c overkit-core-cxx)

# Generate/copy headers to build tree before compiling
add_dependencies(overkit-extras extras-api-headers)
add_dependencies(overkit-extras extras-internal-headers)

#-------------------
# Compiling/linking
#-------------------

# Baseline compiler flags
set(BASE_C_FLAGS_DEBUG -O0 ${DIALECT_C_FLAGS} ${WARNING_C_FLAGS})
set(BASE_C_FLAGS_RELEASE ${DIALECT_C_FLAGS} ${OPT_C_FLAGS})
target_compile_options(overkit-extras PRIVATE
  $<$<CONFIG:Debug>:${BASE_C_FLAGS_DEBUG}>
  $<$<CONFIG:Release>:${BASE_C_FLAGS_RELEASE}>
  $<$<CONFIG:RelWithDebInfo>:${BASE_C_FLAGS_RELEASE}>
  $<$<CONFIG:MinSizeRel>:${BASE_C_FLAGS_RELEASE}>
)

# Coverage
if(COVERAGE)
  target_compile_options(overkit-extras PRIVATE ${COVERAGE_COMPILE_FLAGS})
endif()

# Profiling
if(PROFILE)
  target_compile_options(overkit-extras PRIVATE ${PROFILE_COMPILE_FLAGS})
endif()

# Language feature requirements
if(NOT CMAKE_VERSION VERSION_LESS 3.6)
  target_compile_features(overkit-extras PRIVATE
    cxx_auto_type
    cxx_range_for
  )
endif()

# Needed for linking into shared lib
if(BUILD_SHARED_LIBS)
  set_property(TARGET overkit-extras PROPERTY POSITION_INDEPENDENT_CODE TRUE)
endif()

# MPI
if(EXTERNAL_MPI)
  target_include_directories(overkit-extras SYSTEM PUBLIC "${MPI_INCLUDE_PATH}")
endif()

#--------------------
# API header targets
#--------------------

set(API_HEADER_PREFIX ovk/extras)

set(BUILT_API_HEADERS)
foreach(HEADER_NAME ${API_HEADERS})
  set(BUILT_HEADER "${CMAKE_BINARY_DIR}/include/${API_HEADER_PREFIX}/${HEADER_NAME}")
  set(HEADER "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER_NAME}")
  add_custom_command(OUTPUT "${BUILT_HEADER}"
    MAIN_DEPENDENCY "${HEADER}"
    COMMAND "${CMAKE_COMMAND}" -E copy "${HEADER}" "${BUILT_HEADER}"
  )
  install(FILES "${BUILT_HEADER}" DESTINATION include/${API_HEADER_PREFIX})
  list(APPEND BUILT_API_HEADERS "${BUILT_HEADER}")
endforeach()
set_source_files_properties(${BUILT_API_HEADERS} PROPERTIES GENERATED TRUE)

# Generated header containing includes for all API headers
set(PREFIXED_API_HEADERS)
foreach(HEADER_NAME ${API_HEADERS})
  list(APPEND PREFIXED_API_HEADERS ${API_HEADER_PREFIX}/${HEADER_NAME})
endforeach()
set(BUILT_API_ALL_HEADER "${CMAKE_BINARY_DIR}/include/${API_HEADER_PREFIX}/all.h")
set(API_ALL_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/all.h.in")
add_custom_command(OUTPUT "${BUILT_API_ALL_HEADER}"
  MAIN_DEPENDENCY "${API_ALL_HEADER}"
  COMMAND "${CMAKE_COMMAND}" -DHEADER_IN="${API_ALL_HEADER}" -DHEADER_OUT="${BUILT_API_ALL_HEADER}"
    -DINCLUDE_HEADERS="${PREFIXED_API_HEADERS}" -DSYSTEM=1 -P
    "${CMAKE_SOURCE_DIR}/config/scripts/gen-aggregate-header.cmake"
)
install(FILES "${BUILT_API_ALL_HEADER}" DESTINATION include/${API_HEADER_PREFIX})
set_source_files_properties("${BUILT_API_ALL_HEADER}" PROPERTIES GENERATED TRUE)

add_custom_target(extras-api-headers ALL DEPENDS ${BUILT_API_HEADERS} "${BUILT_API_ALL_HEADER}")
list(APPEND LOCAL_TARGETS extras-api-headers)

#-------------------------
# Internal header targets
#-------------------------

set(INTERNAL_HEADER_PREFIX ovk/extras)

set(BUILT_INTERNAL_HEADERS)
foreach(HEADER_NAME ${INTERNAL_HEADERS})
  set(BUILT_HEADER "${CMAKE_BINARY_DIR}/include-internal/${INTERNAL_HEADER_PREFIX}/${HEADER_NAME}")
  set(HEADER "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER_NAME}")
  add_custom_command(OUTPUT "${BUILT_HEADER}"
    MAIN_DEPENDENCY "${HEADER}"
    COMMAND "${CMAKE_COMMAND}" -E copy "${HEADER}" "${BUILT_HEADER}"
  )
  list(APPEND BUILT_INTERNAL_HEADERS "${BUILT_HEADER}")
endforeach()
set_source_files_properties(${BUILT_INTERNAL_HEADERS} PROPERTIES GENERATED TRUE)
add_custom_target(extras-internal-headers ALL DEPENDS ${BUILT_INTERNAL_HEADERS})
list(APPEND LOCAL_TARGETS extras-internal-headers)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
