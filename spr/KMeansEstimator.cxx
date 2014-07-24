
#include "KMeansEstimator.h"

#include "DynArray.h"
#include "Heap.h"
#include "Log.h"
#include "MersenneTwisterRNG.h"

#include "muException.h"

#include "vnl/vnl_math.h"

KMeansEstimator
::KMeansEstimator()
{

  m_Samples = MatrixType(0, 1);

  m_MaximumIterations = 10;

  m_NumberOfClusters = 2;

  m_NumberOfStarts = 10;

  m_Modified = false;


}

KMeansEstimator
::~KMeansEstimator()
{

}

void
KMeansEstimator
::SetInput(MatrixType& samples)
{

  m_Samples = samples;

  m_Modified = true;
}

void
KMeansEstimator
::SetNumberOfClusters(unsigned int k)
{
  if (k == 0)
    muExceptionMacro(<< "Cannot find 0 clusters");

  m_NumberOfClusters = k;

  m_Modified = true;
}

void
KMeansEstimator
::SetNumberOfStarts(unsigned int s)
{
  if (s == 0)
    muExceptionMacro(<< "Cannot have zero starts");

  m_NumberOfStarts = s;

  m_Modified = true;
}

void
KMeansEstimator
::SetMaximumIterations(unsigned int n)
{
  if (n == 0)
    muExceptionMacro(<< "Cannot have zero max iterations");

  m_MaximumIterations = n;

  m_Modified = true;
}

void
KMeansEstimator
::Update()
{

  unsigned int k = m_NumberOfClusters;
  unsigned int n = m_Samples.rows();

  unsigned int dim = m_Samples.cols();

  if (k > n)
    muExceptionMacro(<<" Cannot find more clusters than data points");

//std::cout << "Kmeans: " << n << " samples, k = " << k << std::endl;

  DynArray<unsigned int> tempLabels;
  tempLabels.Initialize(n, 0);

  MatrixType tempMeans(k, dim, 0);

  double minsum = vnl_huge_val(1.0);

  for (unsigned int istart = 0; istart < m_NumberOfStarts; istart++)
  {
    double sumdist = this->_Guess(tempLabels, tempMeans);

//std::cout << "Start " << istart << ", sum distance = " << sumdist << std::endl;

    if (sumdist < minsum) 
    {
      minsum = sumdist;

      m_Labels = tempLabels;
      m_Means = tempMeans;
    }

  }

//std::cout << "Final within-cluster sum distance = " << minsum << std::endl;

  m_Modified = false;

}

DynArray<unsigned int>
KMeansEstimator
::GetLabels()
{
  if (m_Modified)
    this->Update();

  return m_Labels;
}

KMeansEstimator::MatrixType
KMeansEstimator
::GetMeans()
{
  if (m_Modified)
    this->Update();

  return m_Means;
}

double
KMeansEstimator
::_Guess(DynArray<unsigned int>& labels, MatrixType& means)
{
  unsigned int k = m_NumberOfClusters;
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.cols();

  if (labels.GetSize() != n)
    labels.Initialize(n, 0);

  if (means.rows() != k || means.cols() != dim)
    means.set_size(k, dim);

  // Select sample points at random as initial guess for the means
  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();
  unsigned int* guessIndices = rng->GenerateIntegerSequence(k, n-1);

  means.set_size(k, dim);
  for (unsigned int i = 0; i < k; i++)
    means.set_row(i, m_Samples.get_row(guessIndices[i]));

  delete [] guessIndices;

  //double prevSum = 1.0;

  unsigned int iter = 0;

  // Loop for max likelihood fit
  while (true)
  {
    iter++;

//std::cout << "[KMeansEstimator::_Guess] Iteration " << iter << std::endl;

    bool noLabelChange = true;

    // Classify sample points based on the means
    for (unsigned int isample = 0; isample < n; isample++)
    {
      VectorType x = m_Samples.get_row(isample);

      unsigned int which = 0;
      double mindist = vnl_huge_val(1.0);

      for (unsigned int iclass = 0; iclass < k; iclass++)
      {
        VectorType d = x - means.get_row(iclass);
        double dist = d.squared_magnitude();
        if (dist < mindist)
        {
          which = iclass;
          mindist = dist;
        }
      }

      if (labels[isample] != which)
        noLabelChange = false;

      labels[isample] = which;
    }

    // Check convergence
    if (iter >= m_MaximumIterations)
      break;
    if (noLabelChange)
      break;

/*
    double changeSum = fabs(prevSum - sumdist) / prevSum;
    if (changeSum < 1e-5)
      break;

    prevSum = sumdist;
*/

    // Recompute means
    for (unsigned int iclass = 0; iclass < k; iclass++)
    {
      VectorType mu(dim);
      mu.fill(0);

      unsigned int nclass = 0;

      for (unsigned int isample = 0; isample < n; isample++)
      {
        if (labels[isample] == iclass)
        {
          mu += m_Samples.get_row(isample);
          nclass++;
        }
      }

      if (nclass > 0)
        mu /= nclass;

      means.set_row(iclass, mu);
    }

  }

  // Compute sum of within-cluster distances
  double sumdist = 0;
  for (unsigned int iclass = 0; iclass < k; iclass++)
  {
    VectorType y = means.get_row(iclass);
    for (unsigned int isample = 0; isample < n; isample++)
    {
      if (labels[isample] == iclass)
      {
        VectorType d = y - m_Samples.get_row(isample);
        sumdist += d.squared_magnitude();
      }
    }
  }

  return sumdist;

}
