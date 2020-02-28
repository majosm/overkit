# Copyright (c) 2020 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#================
# Initialization
#================

set(LOCAL_TARGETS)

#=================
# Overkit library
#=================

#----------------
# Library target
#----------------

add_library(overkit
  $<TARGET_OBJECTS:overkit-core>
  $<TARGET_OBJECTS:overkit-core-c>
  $<TARGET_OBJECTS:overkit-extras>
  $<TARGET_OBJECTS:overkit-extras-c>
)
list(APPEND LOCAL_TARGETS overkit)

add_library(Overkit::overkit ALIAS overkit)

set_property(TARGET overkit PROPERTY ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set_property(TARGET overkit PROPERTY LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

target_include_directories(overkit PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_property(TARGET overkit PROPERTY VERSION ${Overkit_VERSION})
set_property(TARGET overkit PROPERTY SOVERSION ${Overkit_VERSION_MAJOR})
set_property(TARGET overkit PROPERTY INTERFACE_Overkit_VERSION_MAJOR ${Overkit_VERSION_MAJOR})
set_property(TARGET overkit APPEND PROPERTY COMPATIBLE_INTERFACE_STRING Overkit_VERSION_MAJOR)

install(TARGETS overkit
  EXPORT OverkitTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include
)

#-------------------
# Compiling/linking
#-------------------

# MPI
if(EXTERNAL_MPI)
  target_include_directories(overkit SYSTEM PUBLIC ${MPI_INCLUDES})
  target_link_libraries(overkit PUBLIC ${MPI_LIBS})
endif()

# OpenMP
if(OPENMP)
  target_link_libraries(overkit PUBLIC ${OpenMP_EXE_LINKER_FLAGS}
    ${OpenMP_CXX_LIB_NAMES} ${OpenMP_C_LIB_NAMES})
endif()

# C math library
target_link_libraries(overkit PUBLIC ${C_MATH_LIBRARY})

# Coverage
if(COVERAGE)
  target_link_libraries(overkit PUBLIC ${COVERAGE_LINK_FLAGS})
endif()

#----------------
# Header targets
#----------------

# C++ header
set(BUILT_OVERKIT_CXX_HEADER ${CMAKE_BINARY_DIR}/include/overkit.hpp)
set(OVERKIT_CXX_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/overkit.hpp)
add_custom_command(OUTPUT ${BUILT_OVERKIT_CXX_HEADER}
  MAIN_DEPENDENCY ${OVERKIT_CXX_HEADER}
  COMMAND "${CMAKE_COMMAND}" -E copy "${OVERKIT_CXX_HEADER}" "${BUILT_OVERKIT_CXX_HEADER}"
)
install(FILES ${BUILT_OVERKIT_CXX_HEADER} DESTINATION include)
set_source_files_properties(${BUILT_OVERKIT_CXX_HEADER} PROPERTIES GENERATED TRUE)
add_custom_target(overkit-cxx-header ALL DEPENDS ${BUILT_OVERKIT_CXX_HEADER})
list(APPEND LOCAL_TARGETS overkit-cxx-header)

# C header
set(BUILT_OVERKIT_C_HEADER ${CMAKE_BINARY_DIR}/include/overkit.h)
set(OVERKIT_C_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/overkit.h)
add_custom_command(OUTPUT ${BUILT_OVERKIT_C_HEADER}
  MAIN_DEPENDENCY ${OVERKIT_C_HEADER}
  COMMAND "${CMAKE_COMMAND}" -E copy "${OVERKIT_C_HEADER}" "${BUILT_OVERKIT_C_HEADER}"
)
install(FILES ${BUILT_OVERKIT_C_HEADER} DESTINATION include)
set_source_files_properties(${BUILT_OVERKIT_C_HEADER} PROPERTIES GENERATED TRUE)
add_custom_target(overkit-c-header ALL DEPENDS ${BUILT_OVERKIT_C_HEADER})
list(APPEND LOCAL_TARGETS overkit-c-header)

#================
# Subdirectories
#================

# Core code
add_subdirectory(core)

# C API for core code
add_subdirectory(core-c)

# Extra stuff
add_subdirectory(extras)

# C API for extras
add_subdirectory(extras-c)

#==============
# Finalization
#==============

# Run pre-build stuff first
foreach(TARGET ${LOCAL_TARGETS})
  add_dependencies(${TARGET} pre-build)
endforeach()
