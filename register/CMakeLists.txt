
PROJECT(register)

INCLUDE (${CMAKE_ROOT}/Modules/FindITK.cmake)
IF (ITK_USE_FILE)
  INCLUDE(${ITK_USE_FILE})
ENDIF (ITK_USE_FILE)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../common
  ${CMAKE_CURRENT_SOURCE_DIR}/../spr
)

ADD_EXECUTABLE(mireg_affine
  ../common/muFile.cxx
  ../common/Log.cxx
  ../common/Timer.cxx
  ../common/MersenneTwisterRNG.cxx
  ../spr/KMeansEstimator.cxx
  ChainedAffineTransform3D.cxx
  GradientDescentOptimizer.cxx
  PowellOptimizer.cxx
  mireg_affine.cxx
)

ADD_EXECUTABLE(resample_label
  ChainedAffineTransform3D.cxx
  resample_label.cxx
)

ADD_EXECUTABLE(test_chaff
  ChainedAffineTransform3D.cxx
  test_chaff.cxx
)

SET(ITK_LIBRARIES
  ITKAlgorithms
  ITKBasicFilters
  ITKCommon
  ITKIO
  ITKMetaIO
  ITKNumerics
  itkvnl_algo
  itkvnl
  itknetlib
)

TARGET_LINK_LIBRARIES(mireg_affine ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(resample_label ${ITK_LIBRARIES})
TARGET_LINK_LIBRARIES(test_chaff ${ITK_LIBRARIES})
