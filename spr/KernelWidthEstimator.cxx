
#include "KernelWidthEstimator.h"

#include "Log.h"
#include "MersenneTwisterRNG.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_qr.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

#define MIN_PROB 1e-10
#define MIN_WIDTH 1e-10

KernelWidthEstimator
::KernelWidthEstimator()
{
  m_LearningRate = 0.1;

  m_MaximumIterations = 1000;

  m_StepLength = 1e-2;

  m_DeterminantGaussianCovariance = 0.0;

  m_NumberOfTestPoints = 100;
}

KernelWidthEstimator
::~KernelWidthEstimator()
{

}

void
KernelWidthEstimator
::SetGaussianParameters(const VectorType& mu, const MatrixType& cov)
{

  // Update Gaussian parameters
  m_GaussianMean = mu;

  m_GaussianCovariance = cov;

  m_InverseGaussianCovariance = vnl_matrix_inverse<float>(cov);

  vnl_qr<float> qr(cov);
  m_DeterminantGaussianCovariance = qr.determinant();

  unsigned int dim = mu.size();

  // Compute square root of cov
  vnl_symmetric_eigensystem<float> eig(m_GaussianCovariance);
  for (unsigned int i = 0; i < dim; i++)
    if (eig.D(i, i) < 1e-10)
      eig.D(i, i) = 1e-10;
  MatrixType sqrtCov = eig.square_root();

  // Generate test points from Gaussian dist
  unsigned int n_test = m_NumberOfTestPoints;

  m_TestSamples = MatrixType(n_test, dim);

  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();
  for (unsigned int k = 0; k < n_test; k++)
  {
    VectorType x(dim);
    for (unsigned int j = 0; j < dim; j++)
      x[j] = rng->GenerateNormal(0.0, 1.0);

    VectorType y = (sqrtCov*x) + m_GaussianMean;

#if DEBUG
    std::cout << "Test sample " << k << " is " << y << std::endl;
#endif

    m_TestSamples.set_row(k, y);
  }
}

double
KernelWidthEstimator
::GetKernelWidth(const MatrixType& samples, double h0)
{

  unsigned int n = samples.rows();

  if (n == 0)
    return 0;

  unsigned int dim = samples.cols();
  
  if (dim == 0)
    return 0;

  if (dim != m_GaussianCovariance.rows())
    return 0;

  m_Samples = samples;

  // Verify sample mean
  VectorType sampleMean(dim, 0.0);
  for (unsigned int i = 0; i < n; i++)
  {
    sampleMean += m_Samples.get_row(i);
  }
  sampleMean /= n;

#if DEBUG
  std::cout << "Sample mean = " << sampleMean << std::endl;
  std::cout << "Gaussian mean = " << m_GaussianMean << std::endl;
#endif

  double h = h0;

  double dh = m_StepLength;
  double dh2 = 2*dh;
  double dh_sq = dh*dh;

  double rate = m_LearningRate;

  double lastMatch = vnl_huge_val(1.0);

  unsigned int iter = 0;

  muLogMacro(<< "Kernel bandwidth search...\n");

  while (true)
  {
    iter++;
    if (iter > m_MaximumIterations)
      break;

    double match = this->ComputeGaussianMatch(h);

    if (iter == 1)
    {
      muLogMacro(<< "  Initial match = " << match
        << ", width = " << h0 << "\n");
    }

#if DEBUG
  if ((iter % 10) == 1)
  {
    std::cout << "iter " << iter << " match = " << match
      << " h = " << h << std::endl;
  }
#endif

    double d_m = this->ComputeGaussianMatchDerivative(h);

    // Gradient descent move
    double mov = d_m * rate;

    // Do simplex search if move is smaller than step size
    if (mov < dh)
    {
      lastMatch = match;

      double forward_m = this->ComputeGaussianMatch(h+dh);
      double backward_m = this->ComputeGaussianMatch(h-dh);

      // Do a simple move if a minimum is found within the simplex
      if (forward_m < match && forward_m < backward_m)
      {
        h = h + dh;
        continue;
      }
      if (backward_m < match && backward_m < forward_m)
      {
        h = h - dh;
        continue;
      }

    }

    // Stop at local minima
    if (fabs(d_m) < 1e-10)
    {
      lastMatch = match;
      break;
    }

    // Update using gradient descent
    h -= mov;

    if (h < MIN_WIDTH)
    {
      h = MIN_WIDTH;
      break;
    }

    lastMatch = match;
  }

  muLogMacro(<< "  After " << iter << " iterations: match = "
    << lastMatch << ", width = " << h << "\n");

  return h;
}

double
KernelWidthEstimator
::ComputeGaussianMatch(double h)
{
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.columns();

  unsigned int n_test = m_TestSamples.rows();

  if (h < MIN_WIDTH)
  {
    h = MIN_WIDTH;
  }

  double denom =
    pow(2*vnl_math::pi, dim/2.0) * sqrt(m_DeterminantGaussianCovariance)
    + vnl_math::eps;

  double match = 0.0;
  for (unsigned int k = 0; k < n_test; k++)
  {
    VectorType x = m_TestSamples.get_row(k);

    VectorType y = m_InverseGaussianCovariance * x;
    double mahalo = 0.0;
    for (unsigned int j = 0; j < dim; j++)
      mahalo += x[j] * y[j];

    double p = exp(-0.5 * mahalo) / denom;

    double q = this->ComputeDensity(h, x) + MIN_PROB;

    double diff = p - q;
    diff = diff * diff;
    match += 1000.0 * diff;

    // Plus entropy of kernel estimate = -(qlogq / n_test)
    //match -= 0.01 * q*log(q) / n_test;

  }

  return match;
}

double
KernelWidthEstimator
::ComputeGaussianMatchDerivative(double h)
{
  unsigned int n = m_Samples.rows();
  unsigned int dim = m_Samples.columns();

  unsigned int n_test = m_TestSamples.rows();

  if (h < MIN_WIDTH)
  {
    h = MIN_WIDTH;
  }

  double denom =
    pow(2*vnl_math::pi, dim/2.0) * sqrt(m_DeterminantGaussianCovariance)
    + vnl_math::eps;

  double h_sq = h*h;
  double h_cube = h_sq*h;

  // Match = \sum{p log(p/q)}
  // where p is the "true" Gaussian dist and q is the estimate
  // deriv = d \sum{p log(p) - p log(q)}
  //   = d \sum{-p log(q)}

  double d_match = 0.0;
  for (int k = 0; k < n_test; k++)
  {
    VectorType x = m_TestSamples.get_row(k);

    VectorType y = m_InverseGaussianCovariance * x;
    double mahalo = 0.0;
    for (unsigned int j = 0; j < dim; j++)
      mahalo += x[j] * y[j];

    double p = exp(-0.5 * mahalo) / denom;

    double q = this->ComputeDensity(h, x) + MIN_PROB;

    double dq = this->ComputeDensityDerivative(h, x);

    d_match -= 2.0 * 1000.0 * (p-q) * dq;

    // Deriv of plogp - plogq - 1/na qlogq is:
    // -p/q dq - 1/na q/q dq - 1/na logq dq
    // -p/q dq - 1/na dq - 1/na logq dq
    // -p/q dq - 1/na dq (1 + logq)

    // Entropy derivative
    //d_match -= 0.01 * (1.0 + log(q)) * dq / n_test;

  }

  return d_match;
}

double
KernelWidthEstimator
::ComputeDensity(double h, const VectorType& x)
{
  double h_sq = h*h;

  unsigned int n = m_Samples.rows();

  double p = 0;
  for (int k = 0; k < n; k++)
  {
    VectorType dvec = m_Samples.get_row(k) - x;

    p += exp(-0.5 * dvec.squared_magnitude() / h_sq);
  }

  unsigned int d = m_Samples.columns();

  p /= n * pow(2*h_sq*vnl_math::pi, d/2.0);

  return p;
}

// Derivative of density w.r.t h at location x
double
KernelWidthEstimator
::ComputeDensityDerivative(double h, const VectorType& x)
{
  double h_sq = h*h;
  double h_cube = h_sq*h;

  unsigned int n = m_Samples.rows();
  unsigned int d = m_Samples.columns();

  double dp = 0;
  for (int k = 0; k < n; k++)
  {
    VectorType dvec = m_Samples.get_row(k) - x;

    double dist_sq = dvec.squared_magnitude();

    double e = exp(-0.5 * dist_sq / h_sq);

    dp += e * (dist_sq / h_cube - d/h);
  }

  dp /= n * pow(2*h_sq*vnl_math::pi, d/2.0);

  return dp;
}
