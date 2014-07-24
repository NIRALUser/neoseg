
PROJECT(spr)

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
IF (USE_ITK_FILE)
  INCLUDE(${USE_ITK_FILE})
ENDIF (USE_ITK_FILE)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
) 

ADD_EXECUTABLE(spr
  ../common/Log.cxx
  ../common/MersenneTwisterRNG.cxx
  KernelWidthEstimator.cxx
  KMeansEstimator.cxx
  KNNClassifier.cxx
  ReducedSetDensityEstimator.cxx
  test.cxx
)

TARGET_LINK_LIBRARIES(spr
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
