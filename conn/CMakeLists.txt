# This is the root ITK CMakeList file: 
cmake_minimum_required(VERSION 2.6)

PROJECT(conn)

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
IF (ITK_USE_FILE)
  INCLUDE(${ITK_USE_FILE})
ENDIF (ITK_USE_FILE)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
)


ADD_EXECUTABLE(conn test.cxx)

#TARGET_LINK_LIBRARIES(conn
# ITKAlgorithms
# ITKBasicFilters
# ITKCommon
# ITKIO
# ITKMetaIO
# ITKNumerics
# itkvnl
# itkvnl_algo
# itknetlib
#)

TARGET_LINK_LIBRARIES(conn ${ITK_LIBRARIES})
