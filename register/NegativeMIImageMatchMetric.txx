
#ifndef _NegativeMIImageMatchMetric_txx
#define _NegativeMIImageMatchMetric_txx

#include "itkContinuousIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"

#include "vnl/vnl_math.h"

#include "NegativeMIImageMatchMetric.h"

#include "KMeansEstimator.h"
#include "MersenneTwisterRNG.h"

#include <float.h>
#include <math.h>

// TODO: ??? use MST instead of KMeans
// stop when reach target # of bins or when exceed max iters

// Image to histogram index mapping using K-means clustering
template <class TImage, class TIndexImage>
typename itk::SmartPointer<TIndexImage>
_kMeansMapIntensityToHistogramIndex(
  const TImage* img, unsigned int numBins, double sampleSpacing)
{
  typename TImage::SizeType size = img->GetLargestPossibleRegion().GetSize();
  typename TImage::SpacingType spacing = img->GetSpacing();

  typename TImage::OffsetType skips; 
  skips[0] = (unsigned)floor(sampleSpacing / spacing[0]);
  skips[1] = (unsigned)floor(sampleSpacing / spacing[1]);
  skips[2] = (unsigned)floor(sampleSpacing / spacing[2]);
  
  if (skips[0] == 0)
    skips[0] = 1;
  if (skips[1] == 0) 
    skips[1] = 1; 
  if (skips[2] == 0)
    skips[2] = 1;

  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  typename TImage::IndexType ind;

  double minv = vnl_huge_val(1.0);
  double maxv = -vnl_huge_val(1.0);

  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        double v = img->GetPixel(ind);
        if (v < minv)
          minv = v;
        if (v > maxv)
          maxv = v;
      }

  double rangev = maxv - minv;

  double t = 0.01 * rangev;
  minv += t;
  maxv -= t;
  rangev -= 2*t;

  // Find intensity clusters
  unsigned long numPossibleSamples = 0;

  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        double v = img->GetPixel(ind);
        if (v < minv || v > maxv)
          continue;
        numPossibleSamples++;
      }

  // Only use 20 % of the samples to find the unique intensities
  unsigned long numSamples = numPossibleSamples / 5;
  //if (numSamples > 50000)
  //  numSamples = 50000;

  unsigned char* selectMask = new unsigned char[numPossibleSamples];
  for (unsigned int i = 0; i < numPossibleSamples; i++)
    selectMask[i] = 0;

  unsigned int* seq =
    rng->GenerateIntegerSequence(numSamples, numPossibleSamples-1);

  for (unsigned int i = 0; i < numSamples; i++)
  {
    unsigned int which = seq[i];
    selectMask[which] = 1;
  }

  delete [] seq;

  DynArray<float> sampleList;

  unsigned int pix = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        double v = img->GetPixel(ind);

        if (v < minv || v > maxv)
          continue;

        if (selectMask[pix] != 0)
        {
          // Intensity samples have to be unique
          bool uniq = true;
          for (unsigned int k = 0; k < sampleList.GetSize(); k++)
            if (sampleList[k] == v)
            {
              uniq = false;
              break;
            }
 
          if (uniq)
            sampleList.Append(v);
        }
        ++pix;
      }

  delete [] selectMask;

  unsigned int numUniqueSamples = sampleList.GetSize();

  KMeansEstimator::MatrixType samples(numUniqueSamples, 1);
  for (unsigned int row = 0; row < numUniqueSamples; row++)
  {
    samples(row, 0) = sampleList[row];
  }

std::cout << "Num unique samples = " << numUniqueSamples << std::endl;

  KMeansEstimator::MatrixType means;

  if (numUniqueSamples <= numBins)
  {
//TODO
//    return _linearMapIntensityToHistogramIndex<TImage, TIndexImage>(
//        img, numBins, sampleSpacing);
    means = samples;
  }
  else
  {
    KMeansEstimator kmeans;
    kmeans.SetMaximumIterations(20);
    kmeans.SetNumberOfClusters(numBins);
    kmeans.SetNumberOfStarts(10);
    kmeans.SetInput(samples);

    means = kmeans.GetMeans();
  }

//std::cout << "K-Means output = \n" << means << std::endl;

  // Allocate index image
  typename itk::SmartPointer<TIndexImage> mapImg = TIndexImage::New();
  mapImg->SetRegions(img->GetLargestPossibleRegion());
  mapImg->Allocate();
  mapImg->SetOrigin(img->GetOrigin());
  mapImg->SetSpacing(img->GetSpacing());

  // Map whole image to histogram index
  typedef itk::ImageRegionConstIterator<TImage> ImageIteratorType;
  ImageIteratorType it(img, img->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<TIndexImage> IndexIteratorType;
  IndexIteratorType mapIt(mapImg, img->GetLargestPossibleRegion());

  it.GoToBegin();
  mapIt.GoToBegin();

  for (; !it.IsAtEnd(); ++it, ++mapIt)
  {
    double v = it.Get();

//TODO: ???
#if 1
    // Ignore extreme values
    if (v < minv || v > maxv)
    {
      mapIt.Set(numBins);
      continue;
    }
#endif

    // Find nearest cluster mean
    double mindist = vnl_huge_val(1.0);
    unsigned int map = numBins;
    for (unsigned int i = 0; i < numBins; i++)
    {
       double d = v - means(i, 0);
       double dist = d*d;
       if (dist < mindist)
       {
         mindist = dist;
         map = i;
       }
    }

    mapIt.Set(map);
  }

  return mapImg;
}

// Image to histogram index mapping using linear mapping
template <class TImage, class TIndexImage>
typename itk::SmartPointer<TIndexImage>
_linearMapIntensityToHistogramIndex(
  const TImage* img, unsigned int numBins, double sampleSpacing)
{
  typename TImage::SizeType size = img->GetLargestPossibleRegion().GetSize();
  typename TImage::SpacingType spacing = img->GetSpacing();

  typename TImage::OffsetType skips; 
  skips[0] = (unsigned)floor(sampleSpacing / spacing[0]);
  skips[1] = (unsigned)floor(sampleSpacing / spacing[1]);
  skips[2] = (unsigned)floor(sampleSpacing / spacing[2]);
  
  if (skips[0] == 0)
    skips[0] = 1;
  if (skips[1] == 0) 
    skips[1] = 1; 
  if (skips[2] == 0)
    skips[2] = 1;

  double minv = vnl_huge_val(1.0);
  double maxv = -vnl_huge_val(1.0);

  typename TImage::IndexType ind;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        double v = img->GetPixel(ind);
        if (v < minv)
          minv = v;
        if (v > maxv)
          maxv = v;
      }

  double rangev = maxv - minv;

// TODO
  double t = 0.01 * rangev;
#if 0
  // Assume MR image with range [0, Inf), clamp extremely bright pixels
  maxv -= t;
  rangev -= t;
#else
  minv += t;
  maxv -= t;
  rangev -= 2*t;
#endif

  // Allocate index image
  typename itk::SmartPointer<TIndexImage> mapImg = TIndexImage::New();
  mapImg->SetRegions(img->GetLargestPossibleRegion());
  mapImg->Allocate();
  mapImg->SetOrigin(img->GetOrigin());
  mapImg->SetSpacing(img->GetSpacing());

  // Map whole image to histogram index
  typedef itk::ImageRegionConstIterator<TImage> ImageIteratorType;
  ImageIteratorType it(img, img->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<TIndexImage> IndexIteratorType;
  IndexIteratorType mapIt(mapImg, img->GetLargestPossibleRegion());

  it.GoToBegin();
  mapIt.GoToBegin();

  for (; !it.IsAtEnd(); ++it, ++mapIt)
  {
    double v = it.Get();

    double u = (v - minv) / rangev;

    unsigned int map = 0;

//TODO: Clamp samples or drop them from histogram computations?
#if 0
    if (u < 0)
      u = 0;
    if (u > 1)
      u = 1;
    map = static_cast<unsigned int>(u * (numBins-1));
#else
    if (u < 0 || u > 1)
      map = numBins;
    else
      map = static_cast<unsigned int>(u * (numBins-1));
#endif

    mapIt.Set(map);
  }

  return mapImg;
}


template <class TFixedImage, class TMovingImage>
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::NegativeMIImageMatchMetric()
{
  m_NumberOfBins = 255;

  m_HistogramPointer = new HistogramType(m_NumberOfBins, m_NumberOfBins);
  m_HistogramPointer->fill(0);

  this->m_FixedImage = 0;
  this->m_MovingImage = 0;

  m_FixedIndexImage = 0;
  m_MovingIndexImage = 0;

  m_KMeansSampleSpacing = 4.0;
  m_SampleSpacing = 4.0;

  m_UseKMeans = false;

  m_Skips[0] = 1;
  m_Skips[1] = 1;
  m_Skips[2] = 1;

  m_Normalized = false;

  m_DerivativeStepLengths = ParametersType(1);
  m_DerivativeStepLengths.Fill(1e-2);
}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "NumberOfBins: ";
  os << m_NumberOfBins << std::endl;
}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::SetFixedImage(
  const typename NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
    ::FixedImageType* img)
{

  itkDebugMacro(<< "SetFixedImage");

  if (img->GetImageDimension() != 3)
    itkExceptionMacro(<< "Fixed image dimension invalid: only supports 3D");

  if (this->m_FixedImage != img)
  {
    this->m_FixedImage = img;
    this->Modified();
  }

  FixedImageSizeType size =
    this->m_FixedImage->GetLargestPossibleRegion().GetSize();
  FixedImageSpacingType spacing = this->m_FixedImage->GetSpacing();

  // Compute skips for downsampling
  m_Skips[0] = (unsigned int)(m_SampleSpacing / spacing[0]);
  m_Skips[1] = (unsigned int)(m_SampleSpacing / spacing[1]);
  m_Skips[2] = (unsigned int)(m_SampleSpacing / spacing[2]);

  if (m_Skips[0] == 0)
    m_Skips[0] = 1;
  if (m_Skips[1] == 0)
    m_Skips[1] = 1;
  if (m_Skips[2] == 0)
    m_Skips[2] = 1;

  this->MapFixedImage();

}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::SetSampleSpacing(double s)
{
  m_SampleSpacing = s;

  if (!this->m_FixedImage.IsNull())
  {
    FixedImageSpacingType spacing = this->m_FixedImage->GetSpacing();

    // Compute skips for downsampling
    m_Skips[0] = (unsigned int)(m_SampleSpacing / spacing[0]);
    m_Skips[1] = (unsigned int)(m_SampleSpacing / spacing[1]);
    m_Skips[2] = (unsigned int)(m_SampleSpacing / spacing[2]);

    if (m_Skips[0] == 0)
      m_Skips[0] = 1;
    if (m_Skips[1] == 0)
      m_Skips[1] = 1;
    if (m_Skips[2] == 0)
      m_Skips[2] = 1;
  }
}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::SetMovingImage(
  const typename NegativeMIImageMatchMetric<TFixedImage,TMovingImage>
  ::MovingImageType* img)
{

  itkDebugMacro(<< "SetMovingImage");

  if (img->GetImageDimension() != 3)
    itkExceptionMacro(<< "Moving image dimension invalid: only supports 3D");

  if (this->m_MovingImage != img)
  {
    this->m_MovingImage = img;
    this->Modified();
  }

  MovingImageSizeType size =
    this->m_MovingImage->GetLargestPossibleRegion().GetSize();
  MovingImageSpacingType spacing = this->m_MovingImage->GetSpacing();

  this->MapMovingImage();

}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::MapFixedImage()
{
  itkDebugMacro(<< "MapFixedImage");

  if (this->m_FixedImage.IsNull())
    return;

  if (m_UseKMeans)
  {
    m_FixedIndexImage =
      _kMeansMapIntensityToHistogramIndex<FixedImageType, IndexImageType>(
        this->m_FixedImage, m_NumberOfBins, m_KMeansSampleSpacing);
  }
  else
  {
    m_FixedIndexImage =
      _linearMapIntensityToHistogramIndex<FixedImageType, IndexImageType>(
        this->m_FixedImage, m_NumberOfBins, m_KMeansSampleSpacing);
  }
}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::MapMovingImage()
{
  itkDebugMacro(<< "MapMovingImage");

  if (this->m_MovingImage.IsNull())
    return;

  if (m_UseKMeans)
  {
    m_MovingIndexImage =
      _kMeansMapIntensityToHistogramIndex<MovingImageType, IndexImageType>(
        this->m_MovingImage, m_NumberOfBins, m_KMeansSampleSpacing);
  }
  else
  {
    m_MovingIndexImage =
      _linearMapIntensityToHistogramIndex<MovingImageType, IndexImageType>(
        this->m_MovingImage, m_NumberOfBins, m_KMeansSampleSpacing);
  }
}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::SetNumberOfBins(unsigned int n)
{

  // Clamp to minimum of 2
  if (n < 2)
  {
    itkWarningMacro(<< "Clamping number of bins to " << 2);
    n = 2;
  }

  unsigned int maxbins = itk::NumericTraits<unsigned int>::max() - 1;

  if (n > maxbins)
  {
    itkWarningMacro(<< "Clamping number of bins to " << maxbins);
    n = maxbins;
  }

  m_NumberOfBins = n;

  delete m_HistogramPointer;
  m_HistogramPointer = new HistogramType(m_NumberOfBins, m_NumberOfBins);
  m_HistogramPointer->fill(0);

  this->MapFixedImage();
  this->MapMovingImage();

  this->Modified();
}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::ComputeHistogram() const
{

  itkDebugMacro(<< "ComputeHistogram");

  HistogramType& H = *m_HistogramPointer;
  H.fill(0);

  FixedImageSizeType fixedSize =
    m_FixedIndexImage->GetLargestPossibleRegion().GetSize();

  MovingImageSizeType movingSize =
    m_MovingIndexImage->GetLargestPossibleRegion().GetSize();

  FixedImageIndexType ind;

  for (ind[2] = 0; ind[2] < fixedSize[2]; ind[2] += m_Skips[2])
    for (ind[1] = 0; ind[1] < fixedSize[1]; ind[1] += m_Skips[1])
      for (ind[0] = 0; ind[0] < fixedSize[0]; ind[0] += m_Skips[0])
      {
        // Get sampled fixed image histogram index
        unsigned int r = m_FixedIndexImage->GetPixel(ind);

        // Skip if fixed image histogram index is invalid
        if (r >= m_NumberOfBins)
        {
          continue;
        }

        FixedImagePointType fixedPoint;
        this->m_FixedImage->TransformIndexToPhysicalPoint(ind, fixedPoint);

        MovingImagePointType mappedPoint =
          this->m_Transform->TransformPoint(fixedPoint);

        // Use Partial Volume interpolation
    
        // Get continuous moving image coordinates (in voxels)
        typedef itk::ContinuousIndex<double, 3> ContinuousIndexType;
        ContinuousIndexType movingInd;
        this->m_MovingImage->TransformPhysicalPointToContinuousIndex(
          mappedPoint, movingInd);

        // Get image neighborhood
        int x0 = (int)movingInd[0];
        int y0 = (int)movingInd[1];
        int z0 = (int)movingInd[2];

        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;

        // Get distances to the image grid
        double fx = movingInd[0] - (double)x0;
        double fy = movingInd[1] - (double)y0;
        double fz = movingInd[2] - (double)z0;

        double gx = 1 - fx;
        double gy = 1 - fy;
        double gz = 1 - fz;

        // Moving image histogram index (column)
        unsigned int c = 0;

// Macro for adding trilinear weights
// Only add if inside moving image and moving index is valid
#define partialVolumeWeightMacro(x, y, z, w) \
  if ((0 <= (x)) && ((x) < movingSize[0]) && \
    (0 <= (y)) && ((y) < movingSize[1]) && \
    (0 <= (z)) && ((z) < movingSize[2])) \
  { \
    MovingImageIndexType pvind = {{(x), (y), (z)}}; \
    c = m_MovingIndexImage->GetPixel(pvind); \
    if (c < m_NumberOfBins) \
      H(r, c) += (w); \
  }

        // Fill histogram with trilinear weights
        partialVolumeWeightMacro(x0, y0, z0, gx*gy*gz);
        partialVolumeWeightMacro(x0, y0, z1, gx*gy*fz);
        partialVolumeWeightMacro(x0, y1, z0, gx*fy*gz);
        partialVolumeWeightMacro(x0, y1, z1, gx*fy*fz);
        partialVolumeWeightMacro(x1, y0, z0, fx*gy*gz);
        partialVolumeWeightMacro(x1, y0, z1, fx*gy*fz);
        partialVolumeWeightMacro(x1, y1, z0, fx*fy*gz);
        partialVolumeWeightMacro(x1, y1, z1, fx*fy*fz);

#undef partialVolumeWeightMacro

      }

}

template <class TFixedImage, class TMovingImage>
double
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::ComputeMI() const
{
  // Compute histogram
  this->ComputeHistogram();

  itkDebugMacro(<< "Start MI");

  HistogramType& H = *m_HistogramPointer;

  double totalf = 0.0;
  for (unsigned c = 0; c < m_NumberOfBins; c++)
  {
    for (unsigned r = 0; r < m_NumberOfBins; r++)
    {
      totalf += H(r, c);
    }
  }

  // All probabilities near zero, -E[log(p)] assumed to be zero
  if (totalf < 1e-20)
    return 0;

  double logtotalf = log(totalf);

  double entropyA = 0.0;
  for (unsigned r = 0; r < m_NumberOfBins; r++)
  {
    double f = 0.0;
    for (unsigned c = 0; c < m_NumberOfBins; c++)
    {
      f += H(r, c);
    }
    if (f > 0.0)
      entropyA += f*log(f);
  }
  // Negate sum and normalize histogram values
  entropyA = -entropyA / totalf + logtotalf;

  double entropyB = 0.0;
  for (unsigned c = 0; c < m_NumberOfBins; c++)
  {
    double f = 0.0;
    for (unsigned r = 0; r < m_NumberOfBins; r++)
    {
      f += H(r, c);
    }
    if (f > 0.0)
      entropyB += f*log(f);
  }
  entropyB = -entropyB / totalf + logtotalf;

  double jointEntropy = 0.0;
  for (unsigned c = 0; c < m_NumberOfBins; c++)
  {
    for (unsigned r = 0; r < m_NumberOfBins; r++)
    {
      double f = H(r, c);
      if (f > 0.0)
        jointEntropy += f*log(f);
    }
  }
  jointEntropy = -jointEntropy / totalf + logtotalf;

  if (m_Normalized)
    return (entropyA + entropyB) / jointEntropy;

  return (entropyA + entropyB) - jointEntropy;

}

template <class TFixedImage, class TMovingImage>
typename NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::MeasureType
NegativeMIImageMatchMetric<TFixedImage,TMovingImage>
::GetValue(const ParametersType& parameters) const
{

  // Make sure the transform has the current parameters
  this->m_Transform->SetParameters(parameters);

  return -1.0*this->ComputeMI();

}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::GetValueAndDerivative(const ParametersType& parameters, MeasureType& value,
  DerivativeType& derivative) const
{
  value = this->GetValue(parameters);
  this->GetDerivative(parameters, derivative);
  //this->GetStochasticDerivative(parameters, derivative);
}

template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::GetDerivative(const ParametersType& parameters, DerivativeType & derivative) const
{
  unsigned int numParams = this->m_Transform->GetNumberOfParameters();

  if (m_DerivativeStepLengths.GetSize() != numParams)
  {
    itkExceptionMacro(<< "Derivative step lengths not set");
  }

  derivative = DerivativeType(numParams);
  derivative.Fill(0);

  for (unsigned int i = 0; i < numParams; i++)
  {
    ParametersType p1 = parameters;
    p1[i] -= m_DerivativeStepLengths[i];
    double v1 = this->GetValue(p1);

    ParametersType p2 = parameters;
    p2[i] += m_DerivativeStepLengths[i];
    double v2 = this->GetValue(p2);

    derivative[i] = (v2 - v1) / (2.0*m_DerivativeStepLengths[i]);
  }
}

// Compute derivative following SPSA (Spall)
template <class TFixedImage, class TMovingImage>
void
NegativeMIImageMatchMetric<TFixedImage, TMovingImage>
::GetStochasticDerivative(const ParametersType& parameters, DerivativeType & derivative) const
{

  unsigned int numParams = this->m_Transform->GetNumberOfParameters();

  if (m_DerivativeStepLengths.GetSize() != numParams)
  {
    itkExceptionMacro(<< "Derivative step lengths not set");
  }

  derivative = DerivativeType(numParams);
  derivative.Fill(0);

  ParametersType dp = ParametersType(numParams);

  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  for (unsigned int i = 0; i < numParams; i++)
  {
    double r = rng->GenerateUniformRealClosedInterval();

    if (r >= 0.5)
      dp[i] = m_DerivativeStepLengths[i];
    else
      dp[i] = -m_DerivativeStepLengths[i];

  }

  ParametersType p1(numParams);
  ParametersType p2(numParams);
  for (unsigned int i = 0; i < numParams; i++)
  {
    p1[i] = parameters[i] + dp[i];
    p2[i] = parameters[i] - dp[i];
  }

  double v1 = this->GetValue(p1);
  double v2 = this->GetValue(p2);

  double v_diff = v1 - v2;

  for (unsigned int i = 0; i < numParams; i++)
  {
    derivative[i] = v_diff / (2.0*dp[i]);
  }
}

#endif
