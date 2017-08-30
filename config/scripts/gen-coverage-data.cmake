# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

set(DEST_DIR "${BINARY_DIR}/coverage/${COVERAGE_DIR}")
file(MAKE_DIRECTORY "${DEST_DIR}")
foreach(SOURCE_EXT c cpp F90)
  file(GLOB_RECURSE GCNO_FILES "${BINARY_DIR}/${COVERAGE_DIR}/CMakeFiles/*.${SOURCE_EXT}.gcno")
  unset(SOURCE_FILES)
  foreach(GCNO_FILE ${GCNO_FILES})
    # Copy source/.gcno/.gcda files to corresponding directory inside coverage dir, removing source
    # extension for .gcno and .gcda files
    get_filename_component(BASE_DIR "${GCNO_FILE}" DIRECTORY)
    get_filename_component(BASE_NAME "${GCNO_FILE}" NAME_WE)
    file(COPY "${GCNO_FILE}" DESTINATION "${DEST_DIR}")
    file(RENAME "${DEST_DIR}/${BASE_NAME}.${SOURCE_EXT}.gcno" "${DEST_DIR}/${BASE_NAME}.gcno")
    set(SOURCE_FILE "${SOURCE_DIR}/${COVERAGE_DIR}/${BASE_NAME}.${SOURCE_EXT}")
    if(EXISTS "${SOURCE_FILE}")
      file(COPY "${SOURCE_FILE}" DESTINATION "${DEST_DIR}")
    endif()
    set(GCDA_FILE "${BASE_DIR}/${BASE_NAME}.${SOURCE_EXT}.gcda")
    if(EXISTS "${GCDA_FILE}")
      file(COPY "${GCDA_FILE}" DESTINATION "${DEST_DIR}")
      file(RENAME "${DEST_DIR}/${BASE_NAME}.${SOURCE_EXT}.gcda" "${DEST_DIR}/${BASE_NAME}.gcda")
    endif()
    # Run gcov on the source files that have coverage data
    if(EXISTS "${SOURCE_FILE}" AND EXISTS "${GCDA_FILE}")
      execute_process(
        COMMAND gcov -b -f "${BASE_NAME}.${SOURCE_EXT}"
        WORKING_DIRECTORY "${DEST_DIR}"
      )
      endif()
  endforeach()
endforeach()
