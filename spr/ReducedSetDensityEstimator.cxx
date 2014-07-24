
#include "Log.h"
#include "ReducedSetDensityEstimator.h"

#include "muException.h"

#include "vnl/vnl_math.h"

#include <float.h>
#include <math.h>


ReducedSetDensityEstimator
::ReducedSetDensityEstimator()
{

  m_Dimension = 1;

  m_InputSet = MatrixType(0, m_Dimension);
  m_ReducedSet = MatrixType(0, m_Dimension);

  m_NumberOfTrainingPoints = 0;

  m_KernelWidth = 1.0;

  m_Modified = false;

}

ReducedSetDensityEstimator
::~ReducedSetDensityEstimator()
{

}

void
ReducedSetDensityEstimator
::SetDimension(unsigned int d)
{
  if (d == m_Dimension)
    return;

  m_Dimension = d;
  m_InputSet = MatrixType(0, d);
  m_ReducedSet = MatrixType(0, d);
  m_Modified = true;
}

void
ReducedSetDensityEstimator
::SetInputSet(const MatrixType& m)
{
  if (m.columns() != m_Dimension)
    muExceptionMacro(
      << "[RSDE::SetInputSet] Dimension of input set does not match "
      << "specified dimension");

  m_InputSet = m;
  m_NumberOfTrainingPoints = m.rows();
  m_Modified = true;
}

void
ReducedSetDensityEstimator
::SetKernelWidth(double h)
{
  if (h == m_KernelWidth)
    return;

  m_KernelWidth = h;
  m_Modified = true;
}

ReducedSetDensityEstimator::VectorType
ReducedSetDensityEstimator
::SequentialMinimalOpt(MatrixType& Q, MatrixType& D)
{

  unsigned int n = Q.rows();

  // Initialize alpha values
  VectorType alpha(n);
  double sumD = FLT_EPSILON;
  for (unsigned int i = 0; i < n; i++)
    sumD += D(i, 0);
  for (unsigned int i = 0; i < n; i++)
    alpha[i] = D(i,0) / sumD;

  VectorType alphaBackup = alpha;

  double alphaTol = 1e-8;
  double errorTol = 1e-6;

  unsigned int iter = 0;

  unsigned int i1 = 0;
  unsigned int i2 = 0;

  double a1 = 0;
  double a2 = 0;

  double E = 0;
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < n; j++)
    {
      E += alpha[i]*alpha[j]*Q(i,j);
    }
  E *= 0.5;
  for (unsigned int i = 0; i < n; i++)
    E -= alpha[i]*D(i, 0);

  muLogMacro(<< "RSDE::SMO initial E = " << E << "\n");

  double old_E = E;

  bool update = false;

  unsigned char* marks = new unsigned char[n];
  for (unsigned int i = 0; i < n; i++)
    if (alpha[i] >= alphaTol)
      marks[i] = 1;
    else
      marks[i] = 0;

  while (true)
  {

    // Find i2, max of current alphas
    unsigned int imax = n;
    double amax = -vnl_huge_val(1.0);
    for (unsigned int i = 0; i < n; i++)
    {
      if ((marks[i] != 0) && (alpha[i] > amax))
      {
        imax = i;
        amax = alpha[i];
      }
    }
    i2 = imax;
    a2 = amax;

    if (a2 >= alphaTol)
    {

      //marks[i2] = 0;
      for (unsigned int j = 0; j < n; j++)
        if (alpha[j] == a2)
          marks[j] = 0;

      // Compute W (I in paper)
      double w2 = -1.0 * D(i2, 0);
      for (unsigned int j = 0; j < n; j++)
      {
        w2 += alpha[j]*Q(i2, j);
      }
      VectorType W(n);
      for (unsigned int i = 0; i < n; i++)
      {
        double s_i = 0;
        for (unsigned int j = 0; j < n; j++)
        {
          s_i += alpha[j]*Q(j, i);
        }
        W[i] = s_i - D(i, 0);
      }
/*
      VectorType W(n);
      for (unsigned int i = 0; i < n; i++)
      {
        double s_i = 0;
        for (unsigned int j = 0; j < n; j++)
        {
          if (j != i1 && j != i2)
            s_i += alpha[j] * Q(i, j);
        }
        W[i] = Q(i1, i)*a1 + Q(i2, i)*a2 + s_i - D(i, 0);
      }
*/

      // Find i1
      unsigned int jmax = n;
      double wmax = -vnl_huge_val(1.0);
      for (unsigned int j = 0; j < n; j++)
      {
        double v = fabs(W[j] - w2);
        //double v = W[j];
        if ((alpha[j] > alphaTol) && (v > wmax))
        {
          jmax = j;
          wmax = v;
        }
      }

      i1 = jmax;
      if (i1 == n)
        a1 = 0;
      else
        a1 = alpha[i1];

      double w1 = W[i1];
      //double w2 = W[i2];

      if (i1 != i2 && a1 >= alphaTol)
      {
        //marks[i1] = 0;
        for (unsigned int j = 0; j < n; j++)
          if (alpha[j] == a1)
            marks[j] = 0;

        // Update a2
        a2 = a2 + (w1 - w2) / (Q(i1, i1) - 2*Q(i1, i2) + Q(i2, i2));

        if (a2 < 0)
          a2 = 0;

        // Update a1
        double d = alpha[i1] + alpha[i2];
        a1 = d - a2;
        if (a1 < 0)
        {
          a1 = 0;
          a2 = d;
        }

        alpha[i1] = a1;
        alpha[i2] = a2;

        update = true;

      } // if a1

    } // if a2

    // Update search marks
    bool allDone = true;
    for (unsigned int i = 0; i < n; i++)
      if (marks[i] != 0)
        allDone = false;

    // Check terminating criterion
    if (allDone)
    {
      // Terminate? (case 2)
      if (!update)
        break;

      // Clear update flag for next time
      update = false;

      // Objective function
      E = 0.0;
      for (unsigned int i = 0; i < n; i++)
        for (unsigned int j = 0; j < n; j++)
        {
          E += alpha[i]*alpha[j]*Q(i,j);
        }
      E *= 0.5;
      for (unsigned int i = 0; i < n; i++)
        E -= alpha[i]*D(i, 0);

      muLogMacro(<< "RSDE::SMO iteration " << iter << " E = " << E << "\n");

      // Terminate? (case 1)
      if (E > old_E)
      {
        // Restore
        alpha = alphaBackup;
        break;
      }
      if (fabs(old_E-E) < errorTol)
        break;

      old_E = E;

      alphaBackup = alpha;

      // Reset marks
      for (unsigned int i = 0; i < n; i++)
        if (alpha[i] >= alphaTol)
          marks[i] = 1;
        else
          marks[i] = 0;

    } // if allDone

    iter++;

  } // while

  delete [] marks;

  return alpha;

}

void
ReducedSetDensityEstimator
::Update()
{

  unsigned int n = m_InputSet.rows();
  unsigned int d = m_InputSet.columns();

  double v = m_KernelWidth*m_KernelWidth;
  double v2 = 2.0*v;

  double scale = 1.0 / pow(2*vnl_math::pi*v, d/2.0);
  double scale2 = 1.0 / pow(2*vnl_math::pi*v2, d/2.0);

  MatrixType Q(n, n);
  MatrixType D(n, 1);
  for (unsigned int i = 0; i < n; i++)
  {
    VectorType x_i = m_InputSet.get_row(i);

    double k_i = 0.0;
    for (unsigned int j = 0; j < n; j++)
    {
      VectorType x_j = m_InputSet.get_row(j);

      VectorType d = x_i - x_j;

      double t = -0.5 * d.squared_magnitude();

      k_i += exp(t / v);
      Q(i, j) = scale2 * exp(t / v2);
    }
    k_i /= n;
    k_i *= scale;

    D(i, 0) = k_i;
  }

  // Compute feature weights
  VectorType alpha_list = this->SequentialMinimalOpt(Q, D);

  // Count number of non-zero weights
  unsigned int numNonZero = 0;
  for (unsigned int k = 0; k < alpha_list.size(); k++)
    if (alpha_list[k] > 0)
      numNonZero++;

  // Create the reduced set
  m_ReducedSet = MatrixType(numNonZero, d);
  unsigned int recent_i = 0;
  for (unsigned int i = 0; i < numNonZero; i++)
  {
    unsigned int other_i = 0;
    unsigned int k;
    for (k = recent_i; k < alpha_list.size(); k++)
    {
      if (alpha_list[k] > 0)
      {
        other_i = k;
        break;
      }
    }
    m_ReducedSet.set_row(i, m_InputSet.get_row(other_i));
    recent_i = k;
  }

  muLogMacro(
    << "RSDE selected " << numNonZero << " / " << n << " training features\n");

  m_Modified = false;

}

double
ReducedSetDensityEstimator
::Evaluate(const VectorType& x)
{
  if (m_Modified)
    this->Update();

  return this->_DensityFunction(m_ReducedSet, x);
}

double
ReducedSetDensityEstimator
::EvaluateWithoutReduce(const VectorType& x)
{
  // No need to update, only use input training set
  //if (m_Modified)
  //  this->Update();

  return this->_DensityFunction(m_InputSet, x);
}

double
ReducedSetDensityEstimator
::_DensityFunction(const MatrixType& set, const VectorType& x)
{

  if (x.size() != m_Dimension)
    muExceptionMacro(<< "[RSDE::_DensityFunction] Invalid input dimension");

  double fx = 0.0;

  unsigned int n = set.rows();

  double v = m_KernelWidth*m_KernelWidth;

  // Don't bother doing exp for samples > 4 std dev away
  double distThres = 4.0*4.0 * v;

  for (unsigned int i = 0; i < n; i++)
  {
    VectorType x_i = set.get_row(i);

    VectorType d = x - x_i;

    double dist_sq = d.squared_magnitude();

    if (dist_sq > distThres)
      continue;

    fx += exp(-0.5 * dist_sq / v);
  }

  fx /= n * pow(2*v*vnl_math::pi, m_Dimension/2.0);

  return fx;

}
