include_directories(${PROJECT_SOURCE_DIR}/../ ${PROJECT_SOURCE_DIR}/ ${PROJECT_SOURCE_DIR}/../../ ${PROJECT_SOURCE_DIR}/src/util/)

file(GLOB jrmath_src
    "${CMAKE_CURRENT_LIST_DIR}/*.h"
    "${CMAKE_CURRENT_LIST_DIR}/*.cc")

add_definitions(-DMATHLIB_STANDALONE)
add_definitions(-DIEEE_754)
  
add_library(jrmath
  ${jrmath_src}
  )
