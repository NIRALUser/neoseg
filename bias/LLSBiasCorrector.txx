
#ifndef _LLSBiasCorrector_txx
#define _LLSBiasCorrector_txx

#include "LLSBiasCorrector.h"

#include "vnl/vnl_math.h"

#include <float.h>
#include <math.h>

#include <iostream>

//#define EXPP(x) (exp((x)/100) - 1)
//#define LOGP(x) (100 * log((x)+1))
#define EXPP(x) (exp(x) - 1)
#define LOGP(x) (log((x)+1))

static inline
double
mypow(double x, unsigned int n)
{
  double p = 1.0;
  for (unsigned int i = 0; i < n; i++)
    p *= x;
  return p;
}

template <class TInputImage, class TProbabilityImage>
LLSBiasCorrector <TInputImage, TProbabilityImage>
::LLSBiasCorrector()
{

  m_InputData = 0;

  m_MaxDegree = 4;

  m_SampleSpacing = 4.0;
  m_WorkingSpacing = 1.0;

  m_ClampBias = false;
  m_MaximumBiasMagnitude = 5.0;

  m_XMu[0] = 0.0;
  m_XMu[1] = 0.0;
  m_XMu[2] = 0.0;

  m_XStd[0] = 1.0;
  m_XStd[1] = 1.0;
  m_XStd[2] = 1.0;

}

template <class TInputImage, class TProbabilityImage>
LLSBiasCorrector <TInputImage, TProbabilityImage>
::~LLSBiasCorrector()
{
  m_Mask = 0;
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector <TInputImage, TProbabilityImage>
::CheckInput()
{

  if (m_InputData.IsNull())
  {
    itkExceptionMacro(<< "Input image not initialized");
  }

  if (m_InputData->GetImageDimension() != 3)
  {
    itkExceptionMacro(<< "Input dimension invalid: only supports 3D images");
  }

  if (m_Probabilities.GetSize() < 1)
  {
    itkExceptionMacro(<< "Must have one or more class probabilities");
  }

  InputImageSizeType size =
    m_InputData->GetLargestPossibleRegion().GetSize();

  for (int i = 0; i < m_Probabilities.GetSize(); i++)
  {
    if (m_Probabilities[i]->GetImageDimension() != 3)
      itkExceptionMacro(<< "Probability [" << i << "] has invalid dimension: only supports 3D images");
    ProbabilityImageSizeType psize =
      m_Probabilities[i]->GetLargestPossibleRegion().GetSize();
    if (size[0] != psize[0] || size[1] != psize[1] || size[2] != psize[2])
      itkExceptionMacro(<< "Image data and probabilities 3D size mismatch");
  }

}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector <TInputImage, TProbabilityImage>
::ComputeDistributions()
{

  itkDebugMacro(<< "LLSBiasCorrector: Computing means and variances...");

  unsigned int numClasses = m_Probabilities.GetSize();

  // Allocate
  m_Means = VectorType(numClasses);
  m_Variances = VectorType(numClasses);

  InputImageSizeType size =
    m_Probabilities[0]->GetLargestPossibleRegion().GetSize();

  InputImageIndexType ind;

  // Compute skips along each dimension
  InputImageSpacingType spacing = m_InputData->GetSpacing();
  unsigned skips[3];
  skips[0] = (unsigned)floor(m_SampleSpacing / spacing[0]);
  skips[1] = (unsigned)floor(m_SampleSpacing / spacing[1]);
  skips[2] = (unsigned)floor(m_SampleSpacing / spacing[2]);

  if (skips[0] == 0)
    skips[0] = 1;
  if (skips[1] == 0)
    skips[1] = 1;
  if (skips[2] == 0)
    skips[2] = 1;

  // Compute the means
  for (unsigned int i = 0; i < numClasses; i++)
  {

    double mu = 0;
    double sumClassProb = DBL_EPSILON;

    for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
      for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
        for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
        {
          mu += 
            m_Probabilities[i]->GetPixel(ind)
            *
            LOGP(m_InputData->GetPixel(ind));
          sumClassProb += m_Probabilities[i]->GetPixel(ind);
        }

    mu /= sumClassProb;

    m_Means[i] = mu;
  }

  // Compute the variances
  for (unsigned int i = 0; i < numClasses; i++)
  {

    double var = 0;
    double diff = 0;
    double sumClassProb = DBL_EPSILON;

    for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
      for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
        for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
        {
          diff = LOGP(m_InputData->GetPixel(ind)) - m_Means[i];
          var += m_Probabilities[i]->GetPixel(ind) * (diff*diff);
          sumClassProb += m_Probabilities[i]->GetPixel(ind);
        }

    var /= sumClassProb;

    m_Variances[i] = var + 1e-20;
  }

#if 1
  // Use relative variance measures (rescale sum to one)
  double sumv = DBL_EPSILON;
  for (int i = 0; i < m_Probabilities.GetSize(); i++)
    sumv += m_Variances[i];
  for (int i = 0; i < m_Probabilities.GetSize(); i++)
    m_Variances[i] = m_Variances[i] / sumv;
#endif

  itkDebugMacro(<< "Means:" << std::endl << m_Means);
  itkDebugMacro(<< "Variances:" << std::endl << m_Variances);

}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector <TInputImage, TProbabilityImage>
::SetMeans(const VectorType& v)
{
  m_Means = VectorType(m_Probabilities.GetSize());
  for (unsigned int k = 0; k < m_Probabilities.GetSize(); k++)
    m_Means[k] = v[k];
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector <TInputImage, TProbabilityImage>
::SetVariances(const VectorType& v)
{
  m_Variances = VectorType(m_Probabilities.GetSize());
  for (unsigned int k = 0; k < m_Probabilities.GetSize(); k++)
    m_Variances[k] = v[k];

#if 1
  // Use relative variance measures (rescale sum to one)
  double sumv = DBL_EPSILON;
  for (int i = 0; i < m_Probabilities.GetSize(); i++)
    sumv += m_Variances[i];
  for (int i = 0; i < m_Probabilities.GetSize(); i++)
    m_Variances[i] = m_Variances[i] / sumv;
#endif

}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector <TInputImage, TProbabilityImage>
::SetMaxDegree(unsigned int n)
{
  itkDebugMacro(<< "SetMaxDegree");

  m_MaxDegree = n;

  // Hack: update LHS
  if (m_Probabilities.GetSize() > 0)
    this->SetProbabilities(m_Probabilities);
}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector <TInputImage, TProbabilityImage>
::SetMask(MaskImageType* mask)
{

  m_Mask = mask;

  InputImageSizeType size =
    m_Mask->GetLargestPossibleRegion().GetSize();

  // Image index for iterations
  InputImageIndexType ind;

  // Compute skips along each dimension
  InputImageSpacingType spacing = m_Mask->GetSpacing();
  unsigned int skips[3];
  skips[0] = (unsigned int)(m_SampleSpacing / spacing[0]);
  skips[1] = (unsigned int)(m_SampleSpacing / spacing[1]);
  skips[2] = (unsigned int)(m_SampleSpacing / spacing[2]);

  if (skips[0] == 0)
    skips[0] = 1;
  if (skips[1] == 0)
    skips[1] = 1;
  if (skips[2] == 0)
    skips[2] = 1;

  itkDebugMacro(<< "Sample skips: " << skips[0] << " x " << skips[1] << " x " << skips[2]);

  unsigned int numCoefficients =
    (m_MaxDegree+1) * (m_MaxDegree+2)/2 * (m_MaxDegree+3)/3;

  // Number of pixels with non-zero weights, downsampled
  unsigned numEquations = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        if (m_Mask->GetPixel(ind) != 0)
          numEquations++;
      }

  itkDebugMacro(<< "Linear system size = " << numEquations << " x " << numCoefficients);

  // Make sure that number of equations >= number of unknowns
  if (numEquations < numCoefficients)
    itkExceptionMacro(<< "Number of unknowns exceed number of equations");

  // Create LHS matrix

  itkDebugMacro(<< "Creating LHS...");

  m_LHS.set_size(numEquations, numCoefficients);

  // Coordinate scaling and offset parameters
  m_XMu[0] = 0.0;
  m_XMu[1] = 0.0;
  m_XMu[2] = 0.0;

  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;

        m_XMu[0] += ind[0];
        m_XMu[1] += ind[1];
        m_XMu[2] += ind[2];
      }
  m_XMu[0] /= numEquations;
  m_XMu[1] /= numEquations;
  m_XMu[2] /= numEquations;

  m_XStd[0] = 0.0;
  m_XStd[1] = 0.0;
  m_XStd[2] = 0.0;

  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;

        double diff;

        diff = ind[0] - m_XMu[0];
        m_XStd[0] += diff*diff;
        diff = ind[1] - m_XMu[1];
        m_XStd[1] += diff*diff;
        diff = ind[2] - m_XMu[2];
        m_XStd[2] += diff*diff;
      }
  m_XStd[0] /= numEquations;
  m_XStd[1] /= numEquations;
  m_XStd[2] /= numEquations;

  m_XStd[0] = sqrt(m_XStd[0]);
  m_XStd[1] = sqrt(m_XStd[1]);
  m_XStd[2] = sqrt(m_XStd[2]);

  // Image coordinate values
  double xc, yc, zc;

  // Row and column indices
  unsigned int r;
  unsigned int c;

  // Fill in LHS matrix with the polynomial basis values
  r = 0;

  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        if (r >= numEquations)
          break;

        if (m_Mask->GetPixel(ind) == 0)
          continue;

        c = 0;
        for (unsigned int order = 0; order <= m_MaxDegree; order++)
          for (unsigned int xorder = 0; xorder <= order; xorder++)
            for (unsigned int yorder = 0; yorder <= (order-xorder); yorder++)
            {
              int zorder = order - xorder - yorder;

              xc = (ind[0] - m_XMu[0]) / m_XStd[0];
              yc = (ind[1] - m_XMu[1]) / m_XStd[1];
              zc = (ind[2] - m_XMu[2]) / m_XStd[2];

              m_LHS(r, c) =
                mypow(xc,xorder) * mypow(yc,yorder) * mypow(zc,zorder);
              c++;
            }

        r++;


      } // for 0

#if 0
  // Weighted LSQ using orthogonal transpose component of LHS
  itkDebugMacro(<< "Computing ortho part of LHS");

  // Note: vnl_qr gives Q mxm and R mxn for A mxn
  MatrixQRType qr(m_LHS);

  // Get economy size R (square)
  MatrixType Rfull = qr.R();
  MatrixType R(numCoefficients, numCoefficients, 0);
  for (unsigned int r = 0; r < numCoefficients; r++)
    for (unsigned int c = r; c < numCoefficients; c++)
      R(r, c) = Rfull(r, c);
  Rfull.set_size(1, 1);

  // Hack to get mxn Q from vnl_qr, Q'*Q = id nxn
  m_LHS_Qt = m_LHS * MatrixInverseType(R);
  m_LHS_Qt.inplace_transpose();
#else
  // Do this instead for ordinary weighted LSQ
  m_LHS_Qt = m_LHS.transpose();
#endif

}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector <TInputImage, TProbabilityImage>
::SetProbabilities(DynArray<ProbabilityImagePointer> probs)
{

  itkDebugMacro(<< "SetProbabilities");

  if (probs.GetSize() < 1)
    itkExceptionMacro(<<"Need one or more probabilities");

  for (int i = 0; i < probs.GetSize(); i++)
  {
    if (probs[i].IsNull())
      itkExceptionMacro(<<"One of input probabilities not initialized");
  }

  m_Probabilities = probs;

}

template <class TInputImage, class TProbabilityImage>
void
LLSBiasCorrector <TInputImage, TProbabilityImage>
::Correct(InputImagePointer input, InputImagePointer output, bool fullRes)
{

  // Obtain a reference count for the bias correction
  m_InputData = input;

  // Verify input
  this->CheckInput();

  InputImageSizeType size =
    m_InputData->GetLargestPossibleRegion().GetSize();

  // Image index for iterations
  InputImageIndexType ind;

  // Verify that the input and output have the same size
  if (size != output->GetLargestPossibleRegion().GetSize())
    itkExceptionMacro(<< "Input and output size mismatch");

  // Compute means and variances
  this->ComputeDistributions();

  // Compute skips along each dimension
  InputImageSpacingType spacing = m_InputData->GetSpacing();

  unsigned int sampleofft[3];
  sampleofft[0] = (unsigned)floor(m_SampleSpacing / spacing[0]);
  sampleofft[1] = (unsigned)floor(m_SampleSpacing / spacing[1]);
  sampleofft[2] = (unsigned)floor(m_SampleSpacing / spacing[2]);

  if (sampleofft[0] < 1)
    sampleofft[0] = 1;
  if (sampleofft[1] < 1)
    sampleofft[1] = 1;
  if (sampleofft[2] < 1)
    sampleofft[2] = 1;

  itkDebugMacro(
     << "Sample offsets: " << sampleofft[0] << " x " << sampleofft[1] << " x " << sampleofft[2]);

  unsigned int workingofft[3];
  workingofft[0] = (unsigned)floor(m_WorkingSpacing / spacing[0]);
  workingofft[1] = (unsigned)floor(m_WorkingSpacing / spacing[1]);
  workingofft[2] = (unsigned)floor(m_WorkingSpacing / spacing[2]);

  if (workingofft[0] < 1)
    workingofft[0] = 1;
  if (workingofft[1] < 1)
    workingofft[1] = 1;
  if (workingofft[2] < 1)
    workingofft[2] = 1;

  itkDebugMacro(<< "Working offsets: " << workingofft[0] << " x " << workingofft[1] << " x " << workingofft[2]);

  unsigned int numClasses = m_Probabilities.GetSize();

  unsigned int numCoefficients =
    (m_MaxDegree+1) * (m_MaxDegree+2)/2 * (m_MaxDegree+3)/3;

  // Create matrices and vectors
  // lhs = replicated basis polynomials for each channel
  // rhs = difference image between original and reconstructed mean image

  itkDebugMacro(<< "Creating matrices for LLS...");

  unsigned int numEquations = m_LHS.rows();

  itkDebugMacro(
    << numEquations << " equations, " << numCoefficients << " coefficients");

  MatrixType lhs(numEquations, numCoefficients);
  //VectorType rhs(numEquations);
  MatrixType rhs(numEquations, 1);

  // Image coordinate values
  double xc, yc, zc;

  // Row and column indices
  unsigned int r;
  unsigned int c;

  itkDebugMacro(<< "Fill rhs");

  // Compute ratio between original and flat image
  // Also compute weights for least squares fitting
  r = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += sampleofft[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += sampleofft[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += sampleofft[0])
      {

        if (r >= numEquations)
          break;

        if (m_Mask->GetPixel(ind) == 0)
          continue;

        // Compute reconstructed intensity and weight (sum of probs / variance)
        double sumProb = DBL_EPSILON;
        double recon = 0;
        double w = 0;
        for (int i = 0; i < numClasses; i++)
        {
          sumProb += m_Probabilities[i]->GetPixel(ind);
          recon += m_Probabilities[i]->GetPixel(ind) * m_Means[i];
          w += m_Probabilities[i]->GetPixel(ind) / m_Variances[i];
          //w += m_Probabilities[i]->GetPixel(ind);
        }
        recon /= sumProb;
        w /= sumProb;

        //double bias = m_InputData->GetPixel(ind) / (recon + 1e-20);
        double bias = LOGP(m_InputData->GetPixel(ind)) - recon;
        //double bias = LOGP(m_InputData->GetPixel(ind)) - LOGP(recon);

        //rhs[r] = w * bias;
        rhs(r, 0) = w * bias;

        for (c = 0; c < numCoefficients; c++)
          lhs(r, c) = w * m_LHS(r, c);

        r++;

      } // for 0

  itkDebugMacro(<< "LSQ solve");
#if 0
  // Linear least squares fit to get polynomial coefficients
  // c = inv(Q'*W*A) * (Q'*W*r)
  // where A is matrix of basis funcs 
  // r is difference
  //MatrixType coeffs = 
  //  MatrixInverseType(m_LHS_Qt * lhs) * (m_LHS_Qt * rhs);
  //MatrixType olhs_inv = MatrixInverseType(m_LHS_Qt * lhs);
  MatrixType olhs = m_LHS_Qt * lhs;
  MatrixType orhs = m_LHS_Qt * rhs;
  //MatrixType coeffs = olhs_inv * orhs;
  MatrixType coeffs = MatrixQRType(olhs).solve(orhs);
#else
  // Use VNL for LSQ
  MatrixType coeffs;
  {
    MatrixQRType qr(lhs);
    coeffs = qr.solve(rhs);
  }
#endif

  itkDebugMacro(<< "Bias field coefficients:" << std::endl  << coeffs);

  // Remove bias
  itkDebugMacro(<< "Correcting input image...");

  if (fullRes)
  {
    workingofft[0] = 1;
    workingofft[1] = 1;
    workingofft[2] = 1;
  }

  double logMax = LOGP(m_MaximumBiasMagnitude);
  double logMin = -1.0 * logMax;
  //double logMax = m_MaximumBiasMagnitude;
  //double logMin = 1.0 / logMax;

  // Compute the log transformed bias field
  InternalImagePointer biasField = InternalImageType::New();
  biasField->SetRegions(input->GetLargestPossibleRegion());
  biasField->Allocate();
  biasField->SetSpacing(input->GetSpacing());

  double maxBias = 0.0;
  double minBias = 0.0;

  for (ind[2] = 0; ind[2] < size[2]; ind[2] += workingofft[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += workingofft[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += workingofft[0])
      {

        double fit = 0.0;

        c = 0;
        for (unsigned int order = 0; order <= m_MaxDegree; order++)
          for (unsigned int xorder = 0; xorder <= order; xorder++)
            for (unsigned int yorder = 0; yorder <= (order-xorder); yorder++)
            {
              int zorder = order - xorder - yorder;

              xc = (ind[0] - m_XMu[0]) / m_XStd[0];
              yc = (ind[1] - m_XMu[1]) / m_XStd[1];
              zc = (ind[2] - m_XMu[2]) / m_XStd[2];

              double poly =
                mypow(xc,xorder) * mypow(yc,yorder) * mypow(zc,zorder);

              //fit += coeffs[c] * poly;
              fit += coeffs(c, 0) * poly;

              c++;
            }

        if (m_ClampBias)
        {
          if (fit < logMin)
            fit = logMin;
          if (fit > logMax)
            fit = logMax;
        }

        if (vnl_math_isnan(fit))
          fit = 0.0;
        if (vnl_math_isinf(fit))
          fit = 0.0;

        if (m_Mask->GetPixel(ind) != 0)
        {
          if (fit > maxBias)
            maxBias = fit;
          if (fit < minBias)
            minBias = fit;
        }

        biasField->SetPixel(ind, (InternalImagePixelType)fit);

      }

  // Correct image using (clamped) bias field)
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += workingofft[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += workingofft[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += workingofft[0])
      {
        double logb = biasField->GetPixel(ind);

        if (logb > maxBias)
          logb = maxBias;
        if (logb < minBias)
          logb = minBias;

        double logd = LOGP(m_InputData->GetPixel(ind)) - logb;
        double d = EXPP(logd);
        //double d = m_InputData->GetPixel(ind) / (logb + 1e-20);

        if (vnl_math_isnan(d))
          d = 0.0;
        if (vnl_math_isinf(d))
          d = 0.0;

        output->SetPixel(ind, (InputImagePixelType)d);
      }

  // Remove reference to input data when done
  m_InputData = 0;

}

#endif
