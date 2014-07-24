
PROJECT(basicimg)

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
IF (USE_ITK_FILE)
  INCLUDE(${USE_ITK_FILE})
ENDIF (USE_ITK_FILE)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
)

SET(ITK_LIBRARIES
  ITKAlgorithms
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKMetaIO
  ITKNumerics
  itkvnl
  itkvnl_algo
  itknetlib
)

ADD_EXECUTABLE(splitaxis splitaxis.cxx)
ADD_EXECUTABLE(testcrop testcrop.cxx)
ADD_EXECUTABLE(testit testit.cxx)

TARGET_LINK_LIBRARIES(splitaxis ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(testcrop ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(testit ${ITK_LIBRARIES})
