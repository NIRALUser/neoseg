
PROJECT(xmlio)

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
IF (USE_ITK_FILE)
  INCLUDE(${USE_ITK_FILE})
ENDIF (USE_ITK_FILE)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../brainseg
  ${CMAKE_CURRENT_SOURCE_DIR}/../neoseg
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
)

ADD_EXECUTABLE(xmlio
  ../brainseg/EMSParameters.cxx
  ../neoseg/NeoSegParameters.cxx
  EMSParametersXMLFile.cxx
  NeoSegParametersXMLFile.cxx
  test.cxx
)

TARGET_LINK_LIBRARIES(xmlio
  ITKAlgorithms
  ITKBasicFilters
  ITKCommon
  ITKEXPAT
  ITKIO
  ITKMetaIO
  ITKNumerics
  itkvnl
  itkvnl_algo
  itknetlib
)
