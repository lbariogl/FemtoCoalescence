#cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(DreamModel)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
  set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG)
include(${ROOT_USE_FILE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${binFolder}/DreamModel)
include_directories(${CMAKE_SOURCE_DIR} ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})
SET(ROOFOLD "$ENV{ROOUNFOLD_ROOT}")

include_directories(${CATS_INCLUDE})
include_directories("${ROOFOLD}")
include_directories(${PROJECT_SOURCE_DIR})

set(SRCS
  ${PROJECT_SOURCE_DIR}/Model.cxx
  ${PROJECT_SOURCE_DIR}/WaveFunction.cxx
  ${PROJECT_SOURCE_DIR}/DLM_CppTools.cpp
  ${PROJECT_SOURCE_DIR}/DLM_Integration.cpp
  )
set(HEADERS
  ${PROJECT_SOURCE_DIR}/Model.h
  ${PROJECT_SOURCE_DIR}/WaveFunction.h
  ${PROJECT_SOURCE_DIR}/DLM_CppTools.h
  ${PROJECT_SOURCE_DIR}/DLM_Integration.h
  ${PROJECT_SOURCE_DIR}/DLM_Histo.h
  ${PROJECT_SOURCE_DIR}/ReadWaveFunction.h
  )
add_library(DreamModel STATIC ${SRCS})
set_target_properties(DreamModel PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(DreamModel ${ROOT_LIBRARIES})

add_executable(excuteTest ${PROJECT_SOURCE_DIR}/Scripts/ExecuteTest.C)
target_link_libraries(excuteTest DreamModel)

