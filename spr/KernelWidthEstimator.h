
////////////////////////////////////////////////////////////////////////////////
//
// Estimate the width for kernel density estimation
//
// Heuristic:
// Find width that minimizes distance to Gaussian dist and entropy
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 11/2005

// NOTE: minimize both entropy and match to a Gaussian dist
// Match cannot be KL since KL is only a proper metric if integrated over
// all state space

#ifndef _KernelWidthEstimator_h
#define _KernelWidthEstimator_h

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include "DynArray.h"

class KernelWidthEstimator
{

public:

  typedef vnl_matrix<float> MatrixType;
  typedef vnl_vector<float> VectorType;

  KernelWidthEstimator();
  ~KernelWidthEstimator();

  void SetStepLength(double d)
  { m_StepLength = d; }

  void SetLearningRate(double r)
  { m_LearningRate = r; }

  void SetMaximumIterations(unsigned int n)
  { m_MaximumIterations = n; }

  void SetNumberOfTestPoints(unsigned int n)
  { m_NumberOfTestPoints = n; }

  void SetGaussianParameters(const VectorType& mu, const MatrixType& cov);

  double GetKernelWidth(const MatrixType& samples, double h0);

private:

  double ComputeDensity(double h, const VectorType& x);
  double ComputeDensityDerivative(double h, const VectorType& x);

  double ComputeGaussianMatch(double h);
  double ComputeGaussianMatchDerivative(double h);

  unsigned int m_MaximumIterations;

  unsigned int m_NumberOfTestPoints;

  MatrixType m_Samples;
  MatrixType m_TestSamples;

  double m_LearningRate;

  double m_StepLength;

  VectorType m_GaussianMean;

  MatrixType m_GaussianCovariance;
  MatrixType m_InverseGaussianCovariance;
  double m_DeterminantGaussianCovariance;

};

#endif
