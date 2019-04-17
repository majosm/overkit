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
  CollectBase.cpp
  CollectCol.cpp
  CollectRow.cpp
  Connectivity.cpp
  ConnectivityD.cpp
  ConnectivityR.cpp
  Context.cpp
  Debug.cpp
  Disperse.cpp
  DisperseBase.cpp
  Domain.cpp
  ErrorHandler.cpp
  Exchange.cpp
  Grid.cpp
  Halo.cpp
  Logger.cpp
  Misc.cpp
  Partition.cpp
  PartitionHash.cpp
  Profiler.cpp
  Recv.cpp
  Request.cpp
  Send.cpp
)

set(PUBLIC_HEADERS
  Array.hpp
  ArrayTraits.hpp
  ArrayView.hpp
  AssemblyOptions.hpp
  Box.hpp
  Box.inl
  Cart.hpp
  Cart.inl
  Connectivity.hpp
  ConnectivityD.hpp
  ConnectivityR.hpp
  Constants.hpp
  Context.hpp
  DataType.hpp
  Domain.hpp
  Elem.hpp
  Elem.inl
  Exchange.hpp
  Global.hpp
  Grid.hpp
  Indexer.hpp
  Interval.hpp
  Misc.hpp
  Misc.inl
  Range.hpp
  Range.inl
  Request.hpp
  ScalarTraits.hpp
  Tuple.hpp
)

set(INTERNAL_HEADERS
  ArrayOps.hpp
  ArrayTraitsBase.hpp
  BoxBase.h
  BoxBase.inl
  Collect.hpp
  CollectAll.hpp
  CollectAny.hpp
  CollectBase.hpp
  CollectInterp.hpp
  CollectNone.hpp
  CollectNotAll.hpp
  Comm.hpp
  Comm.inl
  ConstantsBase.h
  ConstantsBase.inl
  DataTypeBase.h
  DataTypeBase.inl
  Debug.hpp
  DebugBase.h
  Disperse.hpp
  DisperseAppend.hpp
  DisperseBase.hpp
  DisperseOverwrite.hpp
  ErrorHandler.hpp
  ForEach.hpp
  GlobalBase.h
  Halo.hpp
  Halo.inl
  IntegerSequence.hpp
  IteratorTraits.hpp
  Logger.hpp
  Logger.inl
  Partition.hpp
  PartitionBase.hpp
  PartitionHash.hpp
  Profiler.hpp
  Profiler.inl
  Recv.hpp
  Requires.hpp
  ScalarOps.hpp
  ScopeGuard.hpp
  ScopeGuard.inl
  Send.hpp
  TextProcessing.hpp
  TextProcessing.inl
  TypeSequence.hpp
  TypeTraits.hpp
  Vector.hpp
)

#-----------------
# Library targets
#-----------------

add_library(overkit-core OBJECT ${SOURCES})
list(APPEND LOCAL_TARGETS overkit-core)

# Generate/copy headers to build tree before compiling
add_dependencies(overkit-core core-headers)

#-------------------
# Compiling/linking
#-------------------

set(BASE_CXX_FLAGS_DEBUG ${WARNING_CXX_FLAGS})
set(BASE_CXX_FLAGS_RELEASE ${OPT_CXX_FLAGS})
target_compile_options(overkit-core PRIVATE
  $<$<CONFIG:Debug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:SlowDebug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:Release>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:RelWithDebInfo>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:MinSizeRel>:${BASE_CXX_FLAGS_RELEASE}>
)

# Coverage
if(COVERAGE)
  target_compile_options(overkit-core PRIVATE ${COVERAGE_COMPILE_FLAGS})
endif()

# Profiling
if(PROFILE)
  target_compile_options(overkit-core PRIVATE ${PROFILE_COMPILE_FLAGS})
endif()

# Language feature requirements
if(BUILT_IN_DIALECT_SUPPORT)
  if(DIALECT_COMPILE_FEATURE_SUPPORT)
    target_compile_features(overkit-core PRIVATE cxx_std_11)
  else()
    set_property(TARGET overkit-core PROPERTY CXX_STANDARD 11)
  endif()
else()
  target_compile_options(overkit-core PRIVATE ${DIALECT_CXX_FLAGS})
endif()

# Needed for linking into shared lib
if(BUILD_SHARED_LIBS)
  set_property(TARGET overkit-core PROPERTY POSITION_INDEPENDENT_CODE TRUE)
endif()

# MPI
if(EXTERNAL_MPI)
  target_include_directories(overkit-core SYSTEM PUBLIC ${MPI_INCLUDE_PATH})
endif()

#----------------
# Header targets
#----------------

set(BUILT_HEADER_PREFIX ovk/core)

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

# Configuration header
set(BUILT_CONFIG_HEADER ${CMAKE_BINARY_DIR}/include/${BUILT_HEADER_PREFIX}/Config.h)
set(CONFIG_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/Config.h.in)
add_custom_command(OUTPUT ${BUILT_CONFIG_HEADER}
  COMMAND "${CMAKE_COMMAND}"
    -DCONFIG_IN="${CONFIG_HEADER}"
    -DCONFIG_OUT="${BUILT_CONFIG_HEADER}"
    -DCONFIG_OPTIONS=@ONLY
    -DOVK_RELEASE=$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>,$<CONFIG:MinSizeRel>>
    -DOVK_DEBUG=$<OR:$<CONFIG:Debug>,$<CONFIG:SlowDebug>>
    -DOVK_EXAMPLES=${EXAMPLES}
    -DOVK_TESTS=${TESTS}
    -DOVK_COVERAGE=${COVERAGE}
    -DOVK_TIMERS=${TIMERS}
    -DOVK_XPACC=${XPACC}
    -DOVK_HAVE_MPI_IBARRIER=${HAVE_MPI_IBARRIER}
    -P "${CMAKE_SOURCE_DIR}/config/scripts/configure-file.cmake"
)
install(FILES ${BUILT_CONFIG_HEADER} DESTINATION include/${BUILT_HEADER_PREFIX})
set_source_files_properties(${BUILT_CONFIG_HEADER} PROPERTIES GENERATED TRUE)

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

add_custom_target(core-headers ALL DEPENDS ${BUILT_HEADERS} ${BUILT_CONFIG_HEADER}
  ${BUILT_ALL_HEADER})
list(APPEND LOCAL_TARGETS core-headers)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
