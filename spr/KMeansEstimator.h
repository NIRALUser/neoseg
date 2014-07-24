
////////////////////////////////////////////////////////////////////////////////
//
// K-means estimator
//
// Uses multiple random subset selection to initialize the process
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 3/2005

#ifndef _KMeansEstimator_h
#define _KMeansEstimator_h

#include "DynArray.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

class KMeansEstimator
{

public:

  typedef vnl_matrix<float> MatrixType;
  typedef vnl_vector<float> VectorType;

  KMeansEstimator();
  ~KMeansEstimator();

  void SetInput(MatrixType& samples);

  void SetNumberOfClusters(unsigned int c);
  inline unsigned int GetNumberOfClusters() const { return m_NumberOfClusters; }

  void SetNumberOfStarts(unsigned int s);
  inline unsigned int GetNumberOfStarts() const { return m_NumberOfStarts; }

  void SetMaximumIterations(unsigned int n);
  inline unsigned int GetMaximumIterations() const
  { return m_MaximumIterations; }

  MatrixType GetMeans();
  DynArray<unsigned int> GetLabels();

  void Update();

protected:

  double _Guess(DynArray<unsigned int>& labels, MatrixType& means);

private:

  MatrixType m_Samples;

  MatrixType m_Means;
  DynArray<unsigned int> m_Labels;

  unsigned int m_MaximumIterations;

  unsigned int m_NumberOfClusters;
  unsigned int m_NumberOfStarts;

  bool m_Modified;

};

#endif
