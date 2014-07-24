
////////////////////////////////////////////////////////////////////////////////
//
// Reduced Set Density Estimation
// Fast Parzen density by prototype selection
//
// Girolami, He "Probability Density Estimation from Optimally Condensed Data
// Samples" IEEE PAMI Oct 2003 Vol 5 no 10
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 8/2004

#ifndef _ReducedSetDensityEstimator_h
#define _ReducedSetDensityEstimator_h

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

class ReducedSetDensityEstimator
{

public:

  typedef vnl_matrix<float> MatrixType;
  typedef vnl_vector<float> VectorType;

  ReducedSetDensityEstimator();
  ~ReducedSetDensityEstimator();

  void SetDimension(unsigned int d);
  inline unsigned int GetDimension() const { return m_Dimension; }

  void SetInputSet(const MatrixType& m);

  inline unsigned int GetNumberOfTrainingPoints() const
  { return m_NumberOfTrainingPoints; }

  // Get / set Gaussian kernel width (std dev)
  void SetKernelWidth(double h);
  inline double GetKernelWidth() const { return m_KernelWidth; }

  double Evaluate(const VectorType& x);
  double EvaluateWithoutReduce(const VectorType& x);

  void Update();

protected:

  VectorType SequentialMinimalOpt(MatrixType& Q, MatrixType& D);

  double _DensityFunction(const MatrixType& set, const VectorType& x);

private:

  double m_KernelWidth;

  unsigned int m_Dimension;

  unsigned int m_NumberOfTrainingPoints;

  MatrixType m_InputSet;
  MatrixType m_ReducedSet;

  bool m_Modified;

};

#endif
