# This is the root ITK CMakeList file: 
cmake_minimum_required(VERSION 2.6)

PROJECT(NEOSEG)

FIND_PACKAGE(ITK)
IF (ITK_FOUND)
  INCLUDE(${USE_ITK_FILE})
ELSE (ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../aniso
  ${CMAKE_CURRENT_SOURCE_DIR}/../basicimg
  ${CMAKE_CURRENT_SOURCE_DIR}/../bias
  ${CMAKE_CURRENT_SOURCE_DIR}/../brainseg
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
  ${CMAKE_CURRENT_SOURCE_DIR}/../conn
  ${CMAKE_CURRENT_SOURCE_DIR}/../register
  ${CMAKE_CURRENT_SOURCE_DIR}/../robust
  ${CMAKE_CURRENT_SOURCE_DIR}/../spr
  ${CMAKE_CURRENT_SOURCE_DIR}/../xmlio
)

SET(NEOSEG_SRCS
  ../brainseg/EMSParameters.cxx
  ../brainseg/filterFloatImages.cxx
  ../common/Log.cxx
  ../common/MersenneTwisterRNG.cxx
  ../common/Timer.cxx
  ../common/muFile.cxx
  ../register/AtlasRegistrationMethod_float+float.cxx
  ../register/ChainedAffineTransform3D.cxx
  ../register/PairRegistrationMethod_float.cxx
  ../register/PowellOptimizer.cxx
  ../robust/FastMCDSampleFilter.cxx
  ../robust/MSTClusteringProcess.cxx
  ../spr/KernelWidthEstimator.cxx
  ../spr/KMeansEstimator.cxx
  ../spr/KNNClassifier.cxx
  ../spr/ReducedSetDensityEstimator.cxx
  ../xmlio/EMSParametersXMLFile.cxx
  ../xmlio/NeoSegParametersXMLFile.cxx
  NeonateSegmentationFilter.cxx
  NeoSegParameters.cxx
  main.cxx
  runNeonateSeg.cxx
)

ADD_EXECUTABLE(neoseg ${NEOSEG_SRCS})

#TARGET_LINK_LIBRARIES(neoseg
# ITKAlgorithms
# ITKBasicFilters
# ITKCommon
# ITKIO
# ITKMetaIO
# ITKNumerics
# ITKStatistics
# itkvnl
# itkvnl_algo
# itknetlib
#)

TARGET_LINK_LIBRARIES(neoseg ${ITK_LIBRARIES})
