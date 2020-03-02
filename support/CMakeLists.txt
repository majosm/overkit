# Copyright (c) 2020 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#===========================================
# Support library for tests, examples, etc.
#===========================================

#--------------
# Source files
#--------------

set(SOURCES
  CommandArgs.cpp
  CommandArgsC.cpp
  Debug.cpp
  Decomp.cpp
)

set(HEADERS
  CommandArgs.h
  CommandArgs.hpp
  CommandArgsC.h
  Constants.h
  Constants.hpp
  Debug.h
  Debug.hpp
  Decomp.h
  Decomp.hpp
  XDMF.h
  XDMF.hpp
  XDMFC.h
)

# XDMF I/O requires HDF5
if(HAVE_HDF5)
  list(APPEND SOURCES XDMF.cpp XDMFC.cpp)
endif()

#----------------
# Library target
#----------------

add_library(support ${SOURCES})
list(APPEND LOCAL_TARGETS support)

# Generate/copy headers to build tree before compiling
add_dependencies(support pugixml-headers)
add_dependencies(support support-headers)

#-------------------
# Compiling/linking
#-------------------

set(BASE_CXX_FLAGS_DEBUG ${WARNING_CXX_FLAGS})
set(BASE_CXX_FLAGS_RELEASE ${OPT_CXX_FLAGS})
target_compile_options(support PRIVATE
  $<$<CONFIG:SlowDebug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:FastDebug>:${BASE_CXX_FLAGS_DEBUG}>
  $<$<CONFIG:Release>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:RelWithDebInfo>:${BASE_CXX_FLAGS_RELEASE}>
  $<$<CONFIG:MinSizeRel>:${BASE_CXX_FLAGS_RELEASE}>
)

# Profiling
if(PROFILE)
  target_compile_options(support PRIVATE ${PROFILE_COMPILE_FLAGS})
endif()

# Language feature requirements
if(BUILT_IN_DIALECT_SUPPORT)
  if(DIALECT_COMPILE_FEATURE_SUPPORT)
    target_compile_features(support PRIVATE cxx_std_11)
  else()
    set_property(TARGET support PROPERTY CXX_STANDARD 11)
  endif()
else()
  target_compile_options(support PRIVATE ${DIALECT_CXX_FLAGS})
endif()

# Overkit
target_link_libraries(support PRIVATE overkit)

# MPI
if(EXTERNAL_MPI)
  target_include_directories(support SYSTEM PRIVATE ${MPI_INCLUDES})
  target_link_libraries(support PRIVATE ${MPI_LIBS})
endif()

# C math library
target_link_libraries(support PRIVATE ${C_MATH_LIBRARY})

# HDF5
if(HAVE_HDF5)
  target_include_directories(support SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(support PUBLIC ${HDF5_LIBRARIES})
endif()

#----------------
# Header targets
#----------------

set(BUILT_HEADER_PREFIX support)

set(PUGIXML_TAR ${CMAKE_CURRENT_SOURCE_DIR}/pugixml.tar.gz)
set(PUGIXML_DIR ${CMAKE_BINARY_DIR}/include/${BUILT_HEADER_PREFIX}/pugixml)
add_custom_command(OUTPUT ${PUGIXML_DIR}
  MAIN_DEPENDENCY ${PUGIXML_TAR}
  COMMAND "${CMAKE_COMMAND}" -E tar xzf "${PUGIXML_TAR}"
  # Can't do the above command directly into the target directory because CMake doesn't create the
  # directory structure without the copy
  COMMAND "${CMAKE_COMMAND}" -E copy_directory pugixml "${PUGIXML_DIR}"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
set_source_files_properties(
  ${PUGIXML_DIR}/pugiconfig.hpp
  ${PUGIXML_DIR}/pugixml.cpp
  ${PUGIXML_DIR}/pugixml.hpp
  PROPERTIES GENERATED TRUE
)

add_custom_target(pugixml-headers ALL DEPENDS ${PUGIXML_DIR})
list(APPEND LOCAL_TARGETS pugixml-headers)

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

add_custom_target(support-headers ALL DEPENDS ${BUILT_HEADERS})
list(APPEND LOCAL_TARGETS support-headers)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
