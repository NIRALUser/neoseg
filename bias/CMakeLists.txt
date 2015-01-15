
PROJECT(bias)

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
IF (ITK_USE_FILE)
  INCLUDE(${ITK_USE_FILE})
ENDIF (ITK_USE_FILE)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
)

ADD_EXECUTABLE(bias test.cxx)

TARGET_LINK_LIBRARIES(bias
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
