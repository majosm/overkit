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
  Global.cpp
)
if(XPACC)
list(APPEND SOURCES
  XINTOUT.cpp
)
endif()

set(PUBLIC_HEADERS
  Constants.hpp
  Global.hpp
)
if(XPACC)
list(APPEND PUBLIC_HEADERS
  XINTOUT.hpp
)
endif()

set(INTERNAL_HEADERS
  ConstantsBase.h
  ConstantsBase.inl
  GlobalBase.h
)

#-----------------
# Library targets
#-----------------

add_library(overkit-extras OBJECT ${SOURCES})
list(APPEND LOCAL_TARGETS overkit-extras)

# Build core stuff first
add_dependencies(overkit-extras overkit-core)

# Generate/copy headers to build tree before compiling
add_dependencies(overkit-extras extras-headers)

#-------------------
# Compiling/linking
#-------------------

set(BASE_CXX_FLAGS_DEBUG -O0 ${WARNING_CXX_FLAGS})
set(BASE_CXX_FLAGS_RELEASE ${OPT_CXX_FLAGS})
target_compile_options(overkit-extras PRIVATE
  $<$<CONFIG:Debug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:Release>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:RelWithDebInfo>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:MinSizeRel>:${BASE_CXX_FLAGS_RELEASE}>
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
if(NOT CMAKE_VERSION VERSION_LESS 3.8)
  target_compile_features(overkit-extras PRIVATE cxx_std_11)
elseif(NOT CMAKE_VERSION VERSION_LESS 3.1)
  set_property(TARGET overkit-extras PROPERTY CXX_STANDARD 11)
else()
  target_compile_options(overkit-extras PRIVATE ${DIALECT_CXX_FLAGS})
endif()

# Needed for linking into shared lib
if(BUILD_SHARED_LIBS)
  set_property(TARGET overkit-extras PROPERTY POSITION_INDEPENDENT_CODE TRUE)
endif()

# MPI
if(EXTERNAL_MPI)
  target_include_directories(overkit-extras SYSTEM PUBLIC ${MPI_INCLUDE_PATH})
endif()

#----------------
# Header targets
#----------------

set(BUILT_HEADER_PREFIX ovk/extras)

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
set(PREFIXED_PUBLIC_HEADERS)
foreach(HEADER_NAME ${PUBLIC_HEADERS})
  list(APPEND PREFIXED_PUBLIC_HEADERS ${BUILT_HEADER_PREFIX}/${HEADER_NAME})
endforeach()
set(BUILT_ALL_HEADER ${CMAKE_BINARY_DIR}/include/${BUILT_HEADER_PREFIX}/all.hpp)
set(ALL_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/all.hpp.in)
add_custom_command(OUTPUT ${BUILT_ALL_HEADER}
  MAIN_DEPENDENCY ${ALL_HEADER}
  COMMAND "${CMAKE_COMMAND}" -DHEADER_IN="${ALL_HEADER}" -DHEADER_OUT="${BUILT_ALL_HEADER}"
    -DINCLUDE_HEADERS="${PREFIXED_PUBLIC_HEADERS}" -DSYSTEM=1 -P
    "${CMAKE_SOURCE_DIR}/config/scripts/gen-aggregate-header.cmake"
)
install(FILES ${BUILT_ALL_HEADER} DESTINATION include/${BUILT_HEADER_PREFIX})
set_source_files_properties(${BUILT_ALL_HEADER} PROPERTIES GENERATED TRUE)

add_custom_target(extras-headers ALL DEPENDS ${BUILT_HEADERS} ${BUILT_ALL_HEADER})
list(APPEND LOCAL_TARGETS extras-headers)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
