# Copyright (c) 2019 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#=============
# Google Test
#=============

set(GTEST_TAR ${CMAKE_CURRENT_SOURCE_DIR}/googletest.tar.gz)
set(GTEST_SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}/gtest_src)
set(GTEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/gtest_build)

if(NOT EXISTS ${GTEST_SOURCE_DIR})
  message(STATUS "Configuring Google Test")
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E tar xzf "${GTEST_TAR}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

set(BUILD_GMOCK ON CACHE BOOL "")
set(INSTALL_GTEST OFF CACHE BOOL "")
set(gtest_force_shared_crt ON CACHE BOOL "")
set(BUILD_SHARED_LIBS_SAVED ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(${GTEST_SOURCE_DIR} ${GTEST_BINARY_DIR})
set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_SAVED})

#==============================
# Overkit unit test executable
#==============================

#--------------
# Source files
#--------------

set(SOURCES
  ArrayOpsTests.cpp
  ArrayTests.cpp
  ArrayTraitsTests.cpp
  ArrayViewTests.cpp
  BoxTests.cpp
  CartTests.cpp
  CommTests.cpp
  ContextTests.cpp
  DecompTests.cpp
  DistributedFieldTests.cpp
  ElemTests.cpp
  ExchangerTests.cpp
  ForEachTests.cpp
  GeometricPrimitiveOpsTests.cpp
  HaloTests.cpp
  IDTests.cpp
  IndexerTests.cpp
  IntervalTests.cpp
  MapTests.cpp
  MathTests.cpp
  OptionalTests.cpp
  PartitionTests.cpp
  RangeTests.cpp
  RequestTests.cpp
  SetTests.cpp
  TupleTests.cpp
  UnitTestMain.cpp
  VectorTests.cpp
  fixtures/Interface.cpp
)
if(XPACC)
list(APPEND SOURCES
  XINTOUTTests.cpp
)
endif()

set(HEADERS
  MPIPrinter.hpp
  MPITest.hpp
  fixtures/Interface.hpp
  mocks/MultidimArray.hpp
  mocks/Noncopyable.hpp
  mocks/Nondefaultconstructible.hpp
)

#-------------------
# Executable target
#-------------------

add_executable(unit-tests ${SOURCES})
list(APPEND LOCAL_TARGETS unit-tests)

add_test(NAME UnitTests COMMAND "${CMAKE_CURRENT_BINARY_DIR}/unit-tests")

# Generate/copy headers to build tree before compiling
add_dependencies(unit-tests tests-headers)

#-------------------
# Compiling/linking
#-------------------

set(BASE_CXX_FLAGS_DEBUG ${WARNING_CXX_FLAGS})
set(BASE_CXX_FLAGS_RELEASE ${OPT_CXX_FLAGS})
target_compile_options(unit-tests PRIVATE
  $<$<CONFIG:SlowDebug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:FastDebug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:Release>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:RelWithDebInfo>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:MinSizeRel>:${BASE_CXX_FLAGS_RELEASE}>
)

# Profiling
if(PROFILE)
  target_compile_options(unit-tests PRIVATE ${PROFILE_COMPILE_FLAGS})
endif()

# Language feature requirements
if(BUILT_IN_DIALECT_SUPPORT)
  if(DIALECT_COMPILE_FEATURE_SUPPORT)
    target_compile_features(unit-tests PRIVATE cxx_std_11)
  else()
    set_property(TARGET unit-tests PROPERTY CXX_STANDARD 11)
  endif()
else()
  target_compile_options(unit-tests PRIVATE ${DIALECT_CXX_FLAGS})
endif()

# Overkit
target_link_libraries(unit-tests PRIVATE overkit-parallel)

# Support library
target_link_libraries(unit-tests PRIVATE support)

# MPI
if(EXTERNAL_MPI)
  target_include_directories(unit-tests SYSTEM PRIVATE ${MPI_INCLUDES})
  target_link_libraries(unit-tests PRIVATE ${MPI_LIBS})
endif()

# C math library
target_link_libraries(unit-tests PRIVATE ${C_MATH_LIBRARY})

# Google Test
target_link_libraries(unit-tests PRIVATE gtest_main gmock_main)

#----------------
# Header targets
#----------------

set(BUILT_HEADER_PREFIX tests)
set(BUILT_HEADERS)
foreach(HEADER_NAME ${HEADERS})
  set(BUILT_HEADER ${CMAKE_BINARY_DIR}/include/${BUILT_HEADER_PREFIX}/${HEADER_NAME})
  set(HEADER ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER_NAME})
  add_custom_command(OUTPUT ${BUILT_HEADER}
    MAIN_DEPENDENCY ${HEADER}
    COMMAND "${CMAKE_COMMAND}" -E copy "${HEADER}" "${BUILT_HEADER}"
  )
  list(APPEND BUILT_HEADERS ${BUILT_HEADER})
endforeach()
set_source_files_properties(${BUILT_HEADERS} PROPERTIES GENERATED TRUE)

add_custom_target(tests-headers ALL DEPENDS ${BUILT_HEADERS})
list(APPEND LOCAL_TARGETS tests-headers)

#-------------
# Data target
#-------------

set(BUILT_DATA ${CMAKE_BINARY_DIR}/tests/data)
set(DATA ${CMAKE_CURRENT_SOURCE_DIR}/data)
# Seems to run every time if I don't include the touch command; not sure why
add_custom_command(OUTPUT ${BUILT_DATA}
  MAIN_DEPENDENCY ${DATA}
  COMMAND "${CMAKE_COMMAND}" -E copy_directory "${DATA}" "${BUILT_DATA}"
  COMMAND "${CMAKE_COMMAND}" -E touch "${BUILT_DATA}"
)

add_custom_target(tests-data ALL DEPENDS ${BUILT_DATA})
list(APPEND LOCAL_TARGETS tests-data)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
