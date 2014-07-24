
////////////////////////////////////////////////////////////////////////////////
//
// Robust estimation using the fast MCD algorithm
//
// Rousseeuw, P.J., Van Driessen, K.: A fast algorithm for the minimum
// covariance determinant estimator. Technometrics 41 (1999) 212-223
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 3/2004

#ifndef _FastMCDSampleFilter_h
#define _FastMCDSampleFilter_h

#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include "DynArray.h"

class FastMCDSampleFilter
{

public:

  typedef vnl_matrix<float> MatrixType;
  typedef vnl_matrix_inverse<float> MatrixInverseType;

  typedef DynArray<unsigned int> IndexList;

  FastMCDSampleFilter();
  ~FastMCDSampleFilter();

  MatrixType GetInliers(const MatrixType& samples);

  void GetRobustEstimate(MatrixType& mean, MatrixType& covariance,
    const MatrixType& samples);

  void SetChangeTolerance(float f) { m_ChangeTolerance = f; }

  void SetMaxCStepIterations(unsigned int n);

  void SetMahalanobisThreshold(float t);

private:

  float _CSteps(const MatrixType& samples, const IndexList& indices,
    unsigned int maxIters);

  void _ComputeEllipsoid(MatrixType& mean, MatrixType& covariance,
    const MatrixType& samples, const IndexList& indices);

  float m_ChangeTolerance;

  float m_MahalanobisThreshold;

  unsigned int m_MaxCStepIterations;

  unsigned int m_NumberOfStarts;

};

#endif
