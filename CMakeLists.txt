cmake_minimum_required(VERSION 2.6)

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

project(SubdivisionQuadEdgeMeshFilter)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

include(${PROJECT_SOURCE_DIR}/Utilities/KWStyle/KWStyle.cmake)
include_directories(BEFORE ${ITK_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/Source)

set(${PROJECT_BINARY_DIR}/bin CACHE INTERNAL "Where to put the executables") 
set(${PROJECT_BINARY_DIR}/bin CACHE INTERNAL "Where to put the libraries") 

include(CTest)


if(WIN32)
  option(DISABLE_MSVC8_DEPRECATED_WARNINGS
         "Disable Visual Studio 8 deprecated warnings" ON)
  mark_as_advanced(FORCE DISABLE_MSVC8_DEPRECATED_WARNINGS)  
  if(DISABLE_MSVC8_DEPRECATED_WARNINGS)
    add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
  endif(DISABLE_MSVC8_DEPRECATED_WARNINGS)
endif(WIN32)

add_subdirectory( Testing )

