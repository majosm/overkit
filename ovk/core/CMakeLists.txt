# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#==============
# Core library
#==============

#--------------
# Source files
#--------------

set(SOURCES
  AssemblyOptions.cpp
  Connectivity.c
  ConnectivityD.c
  ConnectivityR.c
  Context.c
  Domain.c
  ErrorHandler.c
  Exchange.c
  Global.c
  Grid.c
  Logger.c
  MiscUtils.cpp
  MPIUtils.c
  PartitionHash.c
  ProfileUtils.c
)

set(INTERNAL_HEADERS
  AssemblyOptions.h
  Box.h
  Cart.h
  Connectivity.h
  ConnectivityD.h
  ConnectivityR.h
  Context.h
  Domain.h
  ErrorHandler.h
  Exchange.h
  Global.h
  Grid.h
  List.h
  Logger.h
  MPIUtils.h
  MiscUtils.h
  OrderedMap.h
  PartitionHash.h
  ProfileUtils.h
  Range.h
  TextUtils.h
)

set(API_HEADERS
  ovkAssemblyOptions.h
  ovkBox.h
  ovkBox.inl
  ovkCart.h
  ovkCart.inl
  ovkConnectivity.h
  ovkConnectivityD.h
  ovkConnectivityR.h
  ovkContext.h
  ovkDomain.h
  ovkExchange.h
  ovkGlobal.h
  ovkGrid.h
  ovkRange.h
  ovkRange.inl
)

#----------------
# Library target
#----------------

add_library(overkit-core OBJECT ${SOURCES})
list(APPEND LOCAL_TARGETS overkit-core)

# Generate/copy headers to build tree before compiling
add_dependencies(overkit-core core-api-headers)
add_dependencies(overkit-core core-internal-headers)

#-------------------
# Compiling/linking
#-------------------

# C++11 standard
if(NOT CMAKE_VERSION VERSION_LESS "3.8")
  target_compile_features(overkit-core PRIVATE cxx_std_11)
else()
  target_compile_options(overkit-core PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-std=c++11>)
endif()

# Needed for linking into shared lib
if(BUILD_SHARED_LIBS)
  set_property(TARGET overkit-core PROPERTY POSITION_INDEPENDENT_CODE TRUE)
endif()

# MPI
if(EXTERNAL_MPI)
  target_include_directories(overkit-core SYSTEM PUBLIC "${MPI_INCLUDE_PATH}")
endif()

# Coverage
if(COVERAGE)
  target_compile_options(overkit-core PRIVATE --coverage)
endif()

#--------------------
# API header targets
#--------------------

set(API_HEADER_PREFIX ovk/core)

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

# Configuration header
set(BUILT_CONFIG_HEADER "${CMAKE_BINARY_DIR}/include/${API_HEADER_PREFIX}/ovkConfig.h")
set(CONFIG_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/ovkConfig.h.in")
add_custom_command(OUTPUT "${BUILT_CONFIG_HEADER}"
  COMMAND "${CMAKE_COMMAND}"
    -DCONFIG_IN="${CONFIG_HEADER}"
    -DCONFIG_OUT="${BUILT_CONFIG_HEADER}"
    -DCONFIG_OPTIONS="@ONLY"
    -DOVK_RELEASE=$<CONFIG:Release>
    -DOVK_DEBUG=$<CONFIG:Debug>
    -DOVK_EXAMPLES=${EXAMPLES}
    -DOVK_TESTS=${TESTS}
    -DOVK_COVERAGE=${COVERAGE}
    -DOVK_TIMERS=${TIMERS}
    -DOVK_XPACC=${XPACC}
    -DOVK_HAVE_MPI_IBARRIER=${HAVE_MPI_IBARRIER}
    -P "${CMAKE_SOURCE_DIR}/config/scripts/configure-file.cmake"
)
install(FILES "${BUILT_CONFIG_HEADER}" DESTINATION include/${API_HEADER_PREFIX})
list(APPEND BUILT_API_HEADERS "${BUILT_CONFIG_HEADER}")

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

add_custom_target(core-api-headers ALL DEPENDS ${BUILT_API_HEADERS} "${BUILT_API_ALL_HEADER}")
list(APPEND LOCAL_TARGETS core-api-headers)

#-------------------------
# Internal header targets
#-------------------------

set(INTERNAL_HEADER_PREFIX ovk/core)

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
set_source_files_properties("${BUILT_INTERNAL_HEADERS}" PROPERTIES GENERATED TRUE)
add_custom_target(core-internal-headers ALL DEPENDS ${BUILT_INTERNAL_HEADERS})
list(APPEND LOCAL_TARGETS core-internal-headers)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
