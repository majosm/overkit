# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
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

set(TEST_SOURCES
  ArrayTraitsTests.cpp
  ElemTests.cpp
  ExampleTests.cpp
  IndexerTests.cpp
  IntervalTests.cpp
  RangeTests.cpp
)

#-------------------
# Executable target
#-------------------

add_executable(overkit-unit-tests UnitTestMain.cpp ${TEST_SOURCES})
list(APPEND LOCAL_TARGETS overkit-unit-tests)

add_dependencies(overkit-unit-tests overkit-parallel)

target_include_directories(overkit-unit-tests PRIVATE ${CMAKE_BINARY_DIR}/include)

add_test(NAME UnitTests COMMAND "${CMAKE_CURRENT_BINARY_DIR}/overkit-unit-tests")

#-------------------
# Compiling/linking
#-------------------

set(BASE_CXX_FLAGS_DEBUG -O0 ${WARNING_CXX_FLAGS})
set(BASE_CXX_FLAGS_RELEASE ${OPT_CXX_FLAGS})
target_compile_options(overkit-unit-tests PRIVATE
  $<$<CONFIG:Debug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:Release>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:RelWithDebInfo>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:MinSizeRel>:${BASE_CXX_FLAGS_RELEASE}>
)

# Profiling
if(PROFILE)
  target_compile_options(overkit-unit-tests PRIVATE ${PROFILE_COMPILE_FLAGS})
endif()

# Language feature requirements
if(BUILT_IN_DIALECT_SUPPORT)
  if(DIALECT_COMPILE_FEATURE_SUPPORT)
    target_compile_features(overkit-unit-tests PRIVATE cxx_std_11)
  else()
    set_property(TARGET overkit-unit-tests PROPERTY CXX_STANDARD 11)
  endif()
else()
  target_compile_options(overkit-unit-tests PRIVATE ${DIALECT_CXX_FLAGS})
endif()

# Overkit
target_link_libraries(overkit-unit-tests PRIVATE overkit-parallel)

# MPI
if(EXTERNAL_MPI)
  target_include_directories(overkit-unit-tests SYSTEM PRIVATE ${MPI_INCLUDE_PATH})
  target_link_libraries(overkit-unit-tests PRIVATE ${MPI_C_LIBRARIES})
endif()

# C math library
target_link_libraries(overkit-unit-tests PRIVATE ${C_MATH_LIBRARY})

# Google Test
target_link_libraries(overkit-unit-tests PRIVATE gtest_main gmock_main)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
