
#ifndef _EMSegmentationFilter_txx
#define _EMSegmentationFilter_txx

#include "EMSegmentationFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkNumericTraits.h"

#include "vnl/algo/vnl_qr.h"
#include "vnl/vnl_math.h"

#include "ConnectedComponentsFilter.h"
#include "LLSBiasCorrector.h"
#include "Log.h"
#include "MersenneTwisterRNG.h"
#include "MSTClusteringProcess.h"

#include <iostream>

#include <math.h>
#include <stdlib.h>


template <class TInputImage, class TProbabilityImage>
EMSegmentationFilter <TInputImage, TProbabilityImage>
::EMSegmentationFilter()
{

  m_Mask = 0;

  m_Labels = 0;

  m_SampleSpacing = 4.0;

  // Bias
  m_MaxBiasDegree = 4;
  m_BiasLikelihoodTolerance = 1e-3;

  // EM convergence parameters
  m_LikelihoodTolerance = 1e-5;
  m_MaximumIterations = 40;

  m_InputModified = false;

  m_PriorWeights = VectorType(0);
  m_NumberOfGaussians = 0;

  m_PriorLookupTable = 0;

  m_FOVMask = 0;

  m_DoMSTSplit = false;

}

template <class TInputImage, class TProbabilityImage>
EMSegmentationFilter <TInputImage, TProbabilityImage>
::~EMSegmentationFilter()
{

  delete [] m_NumberOfGaussians;
  delete [] m_PriorLookupTable;

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::CheckInput()
{

  if (m_MaximumIterations == 0)
    itkExceptionMacro(<< "Maximum iterations set to zero");

  if (m_InputImages.GetSize() == 0)
    itkExceptionMacro(<< "No input images");

  if (m_Priors.GetSize() < 1)
    itkExceptionMacro(<< "Must have one or more class probabilities");

  InputImageSizeType size =
    m_InputImages[0]->GetLargestPossibleRegion().GetSize();

  for (unsigned i = 1; i < m_InputImages.GetSize(); i++) {
    if (m_InputImages[i]->GetImageDimension() != 3)
      itkExceptionMacro(<< "InputImage [" << i << "] has invalid dimension: only supports 3D images");
    InputImageSizeType isize =
      m_InputImages[i]->GetLargestPossibleRegion().GetSize();
    if (size != isize)
      itkExceptionMacro(<< "Image data 3D size mismatch");
  }

  for (unsigned i = 0; i < m_Priors.GetSize(); i++)
  {
    if (m_Priors[i]->GetImageDimension() != 3)
      itkExceptionMacro(<< "Prior [" << i << "] has invalid dimension: only supports 3D images");
    ProbabilityImageSizeType psize =
      m_Priors[i]->GetLargestPossibleRegion().GetSize();
    if (size != psize)
      itkExceptionMacro(<< "Image data and prior 3D size mismatch");
  }

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::SetInputImages(DynArray<InputImagePointer> data)
{

  itkDebugMacro(<< "SetInputImages");

  if (data.GetSize() == 0)
    itkExceptionMacro(<< "No input images");

  m_InputImages.Clear();
  for (unsigned i = 0; i < data.GetSize(); i++)
  {
    m_InputImages.Append(data[i]);
  }

  // Allocate space for bias corrected input images
  m_CorrectedImages.Clear();
  InputImageRegionType region = m_InputImages[0]->GetLargestPossibleRegion();
  for (unsigned i = 0; i < data.GetSize(); i++)
  {
    InputImagePointer img = InputImageType::New();
    img->SetRegions(region);
    img->Allocate();
    img->SetOrigin(m_InputImages[0]->GetOrigin());
    img->SetSpacing(m_InputImages[0]->GetSpacing());

    m_CorrectedImages.Append(img);
  }

  m_InputModified = true;

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::SetPriors(DynArray<ProbabilityImagePointer> priors)
{

  itkDebugMacro(<< "SetPriors");

  unsigned int numPriors = priors.GetSize();

  m_Priors.Clear();

  for (unsigned i = 0; i < numPriors; i++)
  {
    m_Priors.Append(priors[i]);
  }

  ProbabilityImageIndexType ind;
  ProbabilityImageSizeType size =
    m_Priors[0]->GetLargestPossibleRegion().GetSize();

  // Clamp minimum to zero
  for (unsigned iprior = 0; iprior < numPriors; iprior++)
  {
    ProbabilityImagePointer tmp = m_Priors[iprior];
    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          if (tmp->GetPixel(ind) < 0)
            tmp->SetPixel(ind, 0);
        }
  }

  // Normalize priors
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        double sumPrior = 0;
        for (unsigned iprior = 0; iprior < numPriors; iprior++)
          sumPrior += m_Priors[iprior]->GetPixel(ind);
        if (sumPrior == 0)
          continue;
        for (unsigned iprior = 0; iprior < numPriors; iprior++)
          m_Priors[iprior]->SetPixel(ind,
            m_Priors[iprior]->GetPixel(ind) / sumPrior);
      }

  m_PriorWeights = VectorType(numPriors);
  m_PriorWeights.fill(1.0);

  delete [] m_NumberOfGaussians;
  m_NumberOfGaussians = new unsigned int[numPriors];
  for (unsigned int i = 0; i < numPriors; i++)
    m_NumberOfGaussians[i] = 1;
  m_NumberOfGaussians[numPriors-1] = 3;

  m_InputModified = true;

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::SetPriorWeights(VectorType w)
{

  itkDebugMacro(<< "SetPriorWeights");

  if (w.size() != m_Priors.GetSize())
    itkExceptionMacro(<< "Number of prior weights invalid");

  for (unsigned i = 0; i < w.size() ; i++)
    if (w[i] == 0.0)
      itkExceptionMacro(<< "Prior weight " << i << " is zero");

  m_PriorWeights = w;

  m_InputModified = true;

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::SetNumberOfGaussians(unsigned int* ng)
{

  itkDebugMacro(<< "SetNumberOfGaussians");

  if (ng == NULL)
    itkExceptionMacro(<< "Number of cluster info invalid");

  unsigned int numPriors = m_Priors.GetSize();

  for (unsigned i = 0; i < numPriors; i++)
    if (ng[i] == 0)
      itkExceptionMacro(<< "Number of Gaussians for prior " << i << " is zero");

  for (unsigned int i = 0; i < numPriors; i++)
    m_NumberOfGaussians[i] = ng[i];

  m_InputModified = true;

}

template <class TInputImage, class TProbabilityImage>
typename EMSegmentationFilter <TInputImage, TProbabilityImage>::
ByteImagePointer
EMSegmentationFilter <TInputImage, TProbabilityImage>
::GetOutput()
{

  this->Update();

  return m_Labels;

}

template <class TInputImage, class TProbabilityImage>
DynArray<
  typename EMSegmentationFilter <TInputImage, TProbabilityImage>::
  ByteImagePointer>
EMSegmentationFilter <TInputImage, TProbabilityImage>
::GetBytePosteriors()
{

  this->Update();

  unsigned int numClasses = m_Posteriors.GetSize();

  DynArray<ByteImagePointer> bytePosts;
  bytePosts.Allocate(numClasses);

  double max = itk::NumericTraits<unsigned char>::max();

  for (unsigned int iclass = 0; iclass < numClasses; iclass++)
  { 
    ProbabilityImagePointer post = m_Posteriors[iclass];
    
    ByteImagePointer tmp = ByteImageType::New();
    tmp->SetRegions(post->GetLargestPossibleRegion());
    tmp->Allocate();
    tmp->SetOrigin(post->GetOrigin());
    tmp->SetSpacing(post->GetSpacing());
    
    ProbabilityImageIndexType ind;

    ProbabilityImageSizeType size = post->GetLargestPossibleRegion().GetSize();

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          double p = floor(max*post->GetPixel(ind) + 0.5);
          if (p < 0)
            p = 0;
          if (p > max)
            p = max;
          tmp->SetPixel(ind, (unsigned char)p);
        }

    bytePosts.Append(tmp);
  }

  return bytePosts;

}

template <class TInputImage, class TProbabilityImage>
DynArray<
  typename EMSegmentationFilter <TInputImage, TProbabilityImage>::
  ShortImagePointer>
EMSegmentationFilter <TInputImage, TProbabilityImage>
::GetShortPosteriors()
{

  this->Update();

  unsigned int numClasses = m_Posteriors.GetSize();

  DynArray<ShortImagePointer> shortPosts;
  shortPosts.Allocate(numClasses);

  double max = itk::NumericTraits<short>::max();

  for (unsigned int iclass = 0; iclass < numClasses; iclass++)
  { 
    ProbabilityImagePointer post = m_Posteriors[iclass];
    
    ShortImagePointer tmp = ShortImageType::New();
    tmp->SetRegions(post->GetLargestPossibleRegion());
    tmp->Allocate();
    tmp->SetOrigin(post->GetOrigin());
    tmp->SetSpacing(post->GetSpacing());
    
    ProbabilityImageIndexType ind;

    ProbabilityImageSizeType size = post->GetLargestPossibleRegion().GetSize();

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          double p = floor(max*post->GetPixel(ind) + 0.5);
          if (p < 0)
            p = 0;
          if (p > max)
            p = max;
          tmp->SetPixel(ind, (short)p);
        }

    shortPosts.Append(tmp);
  }

  return shortPosts;

}

template <class TInputImage, class TProbabilityImage>
DynArray<
  typename EMSegmentationFilter <TInputImage, TProbabilityImage>::
  ProbabilityImagePointer>
EMSegmentationFilter <TInputImage, TProbabilityImage>
::GetPosteriors()
{

  this->Update();

  return m_Posteriors;

}

template <class TInputImage, class TProbabilityImage>
DynArray<
  typename EMSegmentationFilter <TInputImage, TProbabilityImage>::
  InputImagePointer>
EMSegmentationFilter <TInputImage, TProbabilityImage>
::GetCorrected()
{

  this->Update();

  return m_CorrectedImages;

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::Update()
{
  if (!m_InputModified)
    return;

  this->CheckInput();

  this->ComputeMask();

  this->ComputePriorLookupTable();

  this->EMLoop();

  this->ComputeLabels();

  this->CleanUp();

  m_InputModified = false;

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::ComputePriorLookupTable()
{

  itkDebugMacro(<< "ComputePriorLookupTable");

  unsigned int numPriors = m_Priors.GetSize();

  unsigned int numClasses = 0;
  for (unsigned int i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  delete [] m_PriorLookupTable;
  m_PriorLookupTable = new unsigned int[numClasses];

  unsigned int itab = 0;
  for (unsigned int iprior = 0; iprior < numPriors; iprior++)
  {
    unsigned int n = m_NumberOfGaussians[iprior];
    for (unsigned int j = 0; j < n; j++)
    {
      m_PriorLookupTable[itab+j] = iprior;
    }
    itab += n;
  }

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::ComputeMask()
{

  itkDebugMacro(<< "ComputeMask");

  unsigned int numPriors = m_Priors.GetSize();

  m_Mask = ByteImageType::New();
  m_Mask->SetRegions(m_Priors[0]->GetLargestPossibleRegion());
  m_Mask->Allocate();
  m_Mask->SetOrigin(m_Priors[0]->GetOrigin());
  m_Mask->SetSpacing(m_Priors[0]->GetSpacing());

  ProbabilityImageIndexType ind;

  ProbabilityImageSizeType size =
    m_Priors[0]->GetLargestPossibleRegion().GetSize();

  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        double tmp = 0;
        for (unsigned int iprior = 0; iprior < (numPriors-1); iprior++)
          tmp += m_Priors[iprior]->GetPixel(ind);
        if (tmp > 0)
          m_Mask->SetPixel(ind, 1);
        else
          m_Mask->SetPixel(ind, 0);
      }

  if (!m_FOVMask.IsNull())
  {
    if (size != m_FOVMask->GetLargestPossibleRegion().GetSize())
      itkExceptionMacro(<< "FOV mask size mismatch");

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          if (m_FOVMask->GetPixel(ind) == 0)
            m_Mask->SetPixel(ind, 0);
        }
  }

  // Dilate mask
  typedef itk::BinaryBallStructuringElement<unsigned char, 3> StructElementType;
  typedef
    itk::BinaryDilateImageFilter<ByteImageType, ByteImageType,
      StructElementType> DilateType;

  StructElementType structel;
  structel.SetRadius(2);
  structel.CreateStructuringElement();

  typename DilateType::Pointer dil = DilateType::New();
  dil->SetDilateValue(1);
  dil->SetKernel(structel);
  dil->SetInput(m_Mask);

  dil->Update();

  ByteImagePointer dilmask = dil->GetOutput();

  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        m_Mask->SetPixel(ind, dilmask->GetPixel(ind));
      }

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::ComputeDistributions()
{

  unsigned numChannels = m_InputImages.GetSize();
  unsigned numPriors = m_Priors.GetSize();

  unsigned int numClasses = 0;
  for (unsigned int i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  InputImageIndexType ind;

  InputImageSizeType size =
    m_InputImages[0]->GetLargestPossibleRegion().GetSize();

  InputImageSpacingType spacing = m_InputImages[0]->GetSpacing();

  InputImageOffsetType skips;
  skips[0] = (unsigned int)(m_SampleSpacing / spacing[0]);
  skips[1] = (unsigned int)(m_SampleSpacing / spacing[1]);
  skips[2] = (unsigned int)(m_SampleSpacing / spacing[2]);

  if (skips[0] == 0)
    skips[0] = 1;
  if (skips[1] == 0)
    skips[1] = 1;
  if (skips[2] == 0)
    skips[2] = 1;

  // Compute sum of posteriors for each class
  VectorType sumClassProb(numClasses);
  for (unsigned iclass = 0; iclass < numClasses; iclass++)
  {
    double tmp = vnl_math::eps;
    for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
      for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
        for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
        {
          if (m_Mask->GetPixel(ind) == 0)
            continue;
          tmp += (double)m_Posteriors[iclass]->GetPixel(ind);
        }
     sumClassProb[iclass] = tmp;
  }

  // Compute means
  m_Means = MatrixType(numChannels, numClasses);
  for (unsigned int iclass = 0; iclass < numClasses; iclass++)
  {
    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    {

      InputImagePointer img = m_CorrectedImages[ichan];

      double mu = 0;

      // Compute the class mean
      for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
        for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
          for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
          {
            if (m_Mask->GetPixel(ind) == 0)
              continue;
            mu += (double)
              (m_Posteriors[iclass]->GetPixel(ind) * img->GetPixel(ind));
          }

      mu /= sumClassProb[iclass];

      m_Means(ichan, iclass) = mu;

    }
  } // end means loop

  // Fix mean of last class to zero vector (background)
  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    m_Means(ichan, numClasses-1) = 0;

  // Compute covariances
  DynArray<MatrixType> oldCovariances = m_Covariances;
  m_Covariances.Clear();
  for (unsigned int iclass = 0; iclass < numClasses; iclass++)
  {

    MatrixType covtmp(numChannels, numChannels);

    for (unsigned int r = 0; r < numChannels; r++)
    {

      double mu1 = m_Means(r, iclass);
      InputImagePointer img1 = m_CorrectedImages[r];

      for (unsigned int c = r; c < numChannels; c++)
      {

        double mu2 = m_Means(c, iclass);
        InputImagePointer img2 = m_CorrectedImages[c];

        double v = 0;
        double diff1 = 0;
        double diff2 = 0;
        for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
          for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
            for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
            {
              if (m_Mask->GetPixel(ind) == 0)
                continue;
              diff1 = (double)(img1->GetPixel(ind)) - mu1;
              diff2 = (double)(img2->GetPixel(ind)) - mu2;
              v += (double)(m_Posteriors[iclass]->GetPixel(ind)) *
                (diff1*diff2);
            }

        v /= sumClassProb[iclass];

        // Adjust diagonal, to make sure covariance is pos-def
        if (r == c)
          v += 1e-20;

        // Assign value to the covariance matrix (symmetric)
        covtmp(r, c) = v;
        covtmp(c, r) = v;

      }

    }

    vnl_qr<double> qr(covtmp);
    double detcov = qr.determinant();

    // Hack: use old covariance if the new one is singular
    if (detcov <= 0.0)
      m_Covariances.Append(oldCovariances[iclass]);
    else
      m_Covariances.Append(covtmp);

  } // end covariance loop

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::SplitPriorMST(unsigned int iprior)
{

  unsigned int numChannels = m_InputImages.GetSize();
  unsigned int numPriors = m_Priors.GetSize();

  if (iprior >= numPriors)
    itkExceptionMacro(<< "Invalid prior index");

  if (m_NumberOfGaussians[iprior] == 1)
    return; // No split necessary

  unsigned int numClasses = 0;
  for (unsigned int i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  muLogMacro(<< "Splitting distributions for prior " << iprior+1 << "\n");

  InputImageIndexType ind;

  InputImageSizeType size =
    m_InputImages[0]->GetLargestPossibleRegion().GetSize();

  InputImageSpacingType spacing = m_InputImages[0]->GetSpacing();

  InputImageOffsetType skips;
  skips[0] = (unsigned)floor(m_SampleSpacing / spacing[0]);
  skips[1] = (unsigned)floor(m_SampleSpacing / spacing[1]);
  skips[2] = (unsigned)floor(m_SampleSpacing / spacing[2]);

  if (skips[0] == 0)
    skips[0] = 1;
  if (skips[1] == 0)
    skips[1] = 1;
  if (skips[2] == 0)
    skips[2] = 1;

  // Find beginning class index for this prior
  unsigned int istart = 0;
  for (unsigned int j = 0; j < iprior; j++)
    istart += m_NumberOfGaussians[j];

  ProbabilityImagePointer probImg = m_Priors[iprior];

  // Find max prob
  double maxP = 0.0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++) 
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++) 
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++) 
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;
        if (probImg->GetPixel(ind) > maxP)
          maxP = probImg->GetPixel(ind);
      }

  // Select samples by thresholding prior with value above tau
  double tau = 0.85 * maxP;

  muLogMacro(<< "Sampling with threshold tau = " << tau << "\n");

  unsigned int numPossibleSamples = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;
        if (probImg->GetPixel(ind) >= tau)
          numPossibleSamples++;
     }

  // Sample selection mask
  unsigned char* selectMask = new unsigned char[numPossibleSamples];
  for (unsigned int i = 0; i < numPossibleSamples; i++)
    selectMask[i] = 0;

  unsigned int numSamples = numPossibleSamples;
  if (numSamples > 20000)
    numSamples = 20000;

  muLogMacro(<< "  Selecting " << numSamples << " / " << numPossibleSamples << "\n");

  MersenneTwisterRNG rng;

  if (numSamples < numPossibleSamples)
  {
    unsigned int c = 0;
    while (c < numSamples)
    {
      //double r = (double)rand() / (double)RAND_MAX;
      //unsigned int which = (unsigned int)((numPossibleSamples-1) * r);
      unsigned int which = (unsigned int)
        rng.GenerateUniformIntegerUpToK(numPossibleSamples - 1);
      if (selectMask[which] != 0)
        continue;
      selectMask[which] = 1;
      c++;
    }
  }
  else
  {
    for (unsigned int i = 0; i < numSamples; i++)
      selectMask[i] = 1;
  }

  DynArray<MSTClusteringProcess::VertexType> samples;
  samples.Allocate(numSamples);

  muLogMacro(<< "  Finding samples...\n");

  unsigned int r = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;
        if (probImg->GetPixel(ind) < tau)
          continue;

        if (selectMask[r] != 0)
        {
          MSTClusteringProcess::VertexType x(numChannels);
          for (unsigned int ichan = 0; ichan < numChannels; ichan++)
            x[ichan] = m_CorrectedImages[ichan]->GetPixel(ind);
          samples.Append(x);
        }

        ++r;
      }

  delete [] selectMask;

  muLogMacro(<< "  Create MST...\n");

  MSTClusteringProcess mstProc;
  mstProc.SetInputVertices(samples);
  mstProc.SortOn();

  muLogMacro(<< "  Allocate maps...\n");
  unsigned int* clusterMap = new unsigned int[numSamples];

  double T;
  for (T = 2.0; T >= 1.0; T -= 0.01)
  {
    muLogMacro(<< "  MST clustering, T = " << T << "\n");
    unsigned int numClusters = mstProc.GetClusters(clusterMap, T);

    if (numClusters < m_NumberOfGaussians[iprior])
      continue;

    // Check cluster sizes
    bool sizeOK = true;
    for (unsigned int m = 0; m < m_NumberOfGaussians[iprior]; m++)
    {
      unsigned int count = 0;
      for (unsigned int i = 0; i < numSamples; i++)
        if (clusterMap[i] == m)
          count++;
      if (count < 10)
        sizeOK = false;
    }
    
    if (sizeOK)
      break;
  }

  if (T < 1.0)
    itkExceptionMacro(<< "Failed clustering prior " << iprior+1);

  // Use the largest clusters to estimate the Gaussian parameters
  for (unsigned int m = 0; m < m_NumberOfGaussians[iprior]; m++)
  {
    unsigned int iclass = istart+m;

    unsigned int count = 0;
    for (unsigned int i = 0; i < numSamples; i++)
      if (clusterMap[i] == m)
        count++;

    muLogMacro(<< "  Estimating Gaussian parameters for class " << iclass);
    muLogMacro(<< " with " << count << " samples\n");

    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    { 
      double mu = 0.0; 
      for (unsigned int i = 0; i < numSamples; i++)
        if (clusterMap[i] == m)
        {
          MSTClusteringProcess::VertexType x = samples[i];
          mu += x[ichan];
        }
      mu /= count;
      m_Means(ichan, iclass) = mu;
    }

    MatrixType covtmp(numChannels, numChannels);
    for (unsigned int r = 0; r < numChannels; r++)
    { 
      double mu1 = m_Means(r, iclass);
      for (unsigned int c = r; c < numChannels; c++)
      { 
        double mu2 = m_Means(c, iclass);

        double v = 0.0;
        double diff1 = 0;
        double diff2 = 0;
        for (unsigned int i = 0; i < numSamples; i++)
          if (clusterMap[i] == m)
          {
            MSTClusteringProcess::VertexType x = samples[i];
            diff1 = x[r] - mu1;
            diff2 = x[c] - mu2;
            v += (diff1*diff2);
          }
        v /= (count - 1 + vnl_math::eps);

        if (r == c)
          v += 1e-10;

        covtmp(r, c) = v;
        covtmp(c, r) = v;
      }
    }

    // Scale the covariance up a bit for softer initialization
    covtmp *= 1.2;

    m_Covariances[iclass] = covtmp;

  }

  samples.Clear();

  delete [] clusterMap;

  // Special case for background, always set the "darkest" mean
  // to be the last and set it to zero
  if (iprior == (numPriors-1))
  {
    unsigned int imin = istart;
    double minm = m_Means.get_column(istart).squared_magnitude();
    for (unsigned int m = 1; m < m_NumberOfGaussians[iprior]; m++)
    {
      unsigned int iclass = istart+m;
      double mag = m_Means.get_column(iclass).squared_magnitude();
      if (mag < minm)
      {
        imin = iclass;
        minm = mag;
      }
    }

    if (imin != (numClasses-1))
    {
      muLogMacro(
        << "  Replacing " << m_Means.get_column(imin) << " with zero\n");
      VectorType v = m_Means.get_column(numClasses-1);
      m_Means.set_column(imin, v);
    }
    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
      m_Means(ichan, numClasses-1) = 0;
  }

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::SplitPrior(unsigned int iprior)
{

  unsigned int numChannels = m_InputImages.GetSize();
  unsigned int numPriors = m_Priors.GetSize();

  if (iprior >= numPriors)
    itkExceptionMacro(<< "Invalid prior index");

  if (m_NumberOfGaussians[iprior] == 1)
    return; // No split necessary

  unsigned int numClasses = 0;
  for (unsigned int i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  muLogMacro(<< "Splitting distributions for prior " << iprior+1 << "\n");

  // Get the start class index
  unsigned int istart = 0;
  for (unsigned int i = 0; i < iprior; i++)
    istart += m_NumberOfGaussians[i];

  // Scale the means
  for (unsigned int k = 0; k < m_NumberOfGaussians[iprior]; k++)
  {
    double s = 1.0 / pow(1.25, (double)k);
    unsigned int iclass = istart+k;
    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
      m_Means(ichan, iclass) = s * m_Means(ichan, iclass);
  }

  // Set mean of last class to zero
  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    m_Means(ichan, numClasses-1) = 0;

}


template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::ComputePosteriors(bool fullRes)
{

  itkDebugMacro(<< "ComputePosteriors");

  unsigned numChannels = m_InputImages.GetSize();
  unsigned numPriors = m_Priors.GetSize();

  unsigned int numClasses = 0;
  for (unsigned i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  if (fullRes)
    muLogMacro(<< "Computing posteriors at full resolution\n");

  ProbabilityImageIndexType ind;

  ProbabilityImageSizeType size =
    m_Priors[0]->GetLargestPossibleRegion().GetSize();

  ProbabilityImageSpacingType spacing = m_Priors[0]->GetSpacing();

  ProbabilityImageOffsetType skips;
  skips[0] = (unsigned int)(m_SampleSpacing / spacing[0]);
  skips[1] = (unsigned int)(m_SampleSpacing / spacing[1]);
  skips[2] = (unsigned int)(m_SampleSpacing / spacing[2]);

  if (skips[0] == 0 || fullRes)
    skips[0] = 1;
  if (skips[1] == 0 || fullRes)
    skips[1] = 1;
  if (skips[2] == 0 || fullRes)
    skips[2] = 1;

  for (unsigned int iclass = 0; iclass < numClasses; iclass++)
  {

    ProbabilityImagePointer post = m_Posteriors[iclass];
    ProbabilityImagePointer prior = m_Priors[m_PriorLookupTable[iclass]];

    double priorScale =
      m_PriorWeights[m_PriorLookupTable[iclass]]
      /
      m_NumberOfGaussians[m_PriorLookupTable[iclass]];

    vnl_qr<double> qr(m_Covariances[iclass]);
    double detcov = qr.determinant();

    if (detcov <= 0.0)
      itkExceptionMacro(<< "Determinant of covariance for class " << iclass
        << " is <= 0.0 (" << detcov << "), covariance matrix:\n"
        << m_Covariances[iclass]);

    // Normalizing constant for the Gaussian
    double denom =
      pow(2*vnl_math::pi, m_InputImages.GetSize()/2.0) * sqrt(detcov)
      + vnl_math::eps;
      
    //MatrixType invcov = MatrixInverseType(m_Covariances[iclass]);
    MatrixType id(numChannels, numChannels);
    id.set_identity();
    MatrixType invcov = qr.solve(id);

    for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
      for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
        for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
        {
          if (m_Mask->GetPixel(ind) == 0)
          {
            post->SetPixel(ind, 0);
            continue;
          }

          MatrixType X(numChannels, 1);
          for (unsigned int ichan = 0; ichan < numChannels; ichan++)
            X(ichan, 0) =
              m_CorrectedImages[ichan]->GetPixel(ind) - m_Means(ichan, iclass);

          MatrixType Y = invcov * X;

          double mahalo = 0.0;
          for (unsigned int ichan = 0; ichan < numChannels; ichan++)
            mahalo += X(ichan, 0) * Y(ichan, 0);

          double likelihood = exp(-0.5 * mahalo) / denom;
          
          post->SetPixel(ind, (ProbabilityImagePixelType)
            (priorScale * prior->GetPixel(ind) * likelihood));
        } // for 0


  } // end class loop

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::CorrectBias(unsigned int degree, bool fullRes)
{

  muLogMacro(<< "Bias correction, polynomial degree = " << degree << "\n");

  unsigned int numPriors = m_Priors.GetSize();

  unsigned int numClasses = 0;
  for (unsigned int i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  unsigned int numFGClasses = numClasses - m_NumberOfGaussians[numPriors-1];

  // Perform bias correction
  DynArray<ProbabilityImagePointer> biasPosteriors;
  for (unsigned j = 0; j < (numClasses-1); j++)
    biasPosteriors.Append(m_Posteriors[j]);

  typedef LLSBiasCorrector<InputImageType, ProbabilityImageType>
    BiasCorrectorType;
  typedef typename BiasCorrectorType::Pointer BiasCorrectorPointer;

  BiasCorrectorPointer biascorr = BiasCorrectorType::New();

  biascorr->SetClampBias(false);
  biascorr->SetMaximumBiasMagnitude(5.0);
  biascorr->SetMaxDegree(degree);
  biascorr->SetSampleSpacing(m_SampleSpacing);

  biascorr->SetMask(m_Mask);
  biascorr->SetProbabilities(biasPosteriors);

  if (this->GetDebug())
    biascorr->DebugOn();

  for (unsigned int ichan = 0; ichan < m_InputImages.GetSize(); ichan++)
  {
    itkDebugMacro(<< "Correcting image " << ichan+1 << "...");

// Let the bias corrector compute the distribution
// (uses log scaled intensities)
#if 0
    BiasCorrectorType::VectorType means(numClasses-1);
    for (unsigned int k = 0; k < (numClasses-1); k++)
      means[k] = m_Means(ichan, k);
    biascorr->SetMeans(means);

    BiasCorrectorType::VectorType vars(numClasses-1);
    for (unsigned int k = 0; k < (numClasses-1); k++)
      vars[k] = (m_Covariances[k])(ichan, ichan);
    biascorr->SetVariances(vars);
#endif

    biascorr->Correct(m_InputImages[ichan], m_CorrectedImages[ichan],
      fullRes);
  }

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::EMLoop()
{

  unsigned int numChannels = m_InputImages.GetSize();
  unsigned int numPriors = m_Priors.GetSize();

  unsigned int numClasses = 0;
  for (unsigned int i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  InputImageIndexType ind;

  InputImageSizeType size =
    m_InputImages[0]->GetLargestPossibleRegion().GetSize();

  InputImageSpacingType spacing = m_InputImages[0]->GetSpacing();

  InputImageOffsetType skips;
  skips[0] = (unsigned)floor(m_SampleSpacing / spacing[0]);
  skips[1] = (unsigned)floor(m_SampleSpacing / spacing[1]);
  skips[2] = (unsigned)floor(m_SampleSpacing / spacing[2]);

  if (skips[0] == 0)
    skips[0] = 1;
  if (skips[1] == 0)
    skips[1] = 1;
  if (skips[2] == 0)
    skips[2] = 1;

  // Allocate space for posteriors
  for (unsigned j = 0; j < m_Posteriors.GetSize(); j++)
  {
    m_Posteriors[j] = 0;
  }
  m_Posteriors.Clear();
  for (unsigned iclass = 0; iclass < numClasses; iclass++)
  {
    ProbabilityImagePointer post = ProbabilityImageType::New();
    post->SetRegions(m_InputImages[0]->GetLargestPossibleRegion());
    post->Allocate();
    post->SetOrigin(m_InputImages[0]->GetOrigin());
    post->SetSpacing(m_InputImages[0]->GetSpacing());

    m_Posteriors.Append(post);
  }

  // Initialize first iteration posteriors with priors
  for (unsigned iclass = 0; iclass < numClasses; iclass++)
  {
     ProbabilityImagePointer post = m_Posteriors[iclass];
     ProbabilityImagePointer prior = m_Priors[m_PriorLookupTable[iclass]];
     double scale =
       1.0 / m_NumberOfGaussians[m_PriorLookupTable[iclass]];
     for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
       for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
         for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
         {
           if (m_Mask->GetPixel(ind) == 0)
             post->SetPixel(ind, 0);
           else
             post->SetPixel(ind, (ProbabilityImagePixelType)
               (prior->GetPixel(ind) * scale));
         }
  }

  // Compute the reference mean (first class)
  VectorType refMean(numChannels);
  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
  {
    InputImagePointer img = m_InputImages[ichan];
    double mu = 0;
    double sumfirst = vnl_math::eps;
    for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
      for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
        for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
        {
          mu += (double)
            (m_Posteriors[0]->GetPixel(ind) * img->GetPixel(ind));
          sumfirst += (double)m_Posteriors[0]->GetPixel(ind);
        }
    mu /= sumfirst;
    refMean(ichan) = mu;
  }

  // Initial copy, since the computation of the distribution
  // parameters is done using the corrected images
  for (unsigned int i = 0; i < m_InputImages.GetSize(); i++)
  {
    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          double v = m_InputImages[i]->GetPixel(ind);
          m_CorrectedImages[i]->SetPixel(ind, v);
        }
  }

  // Compute initial distribution parameters
  this->ComputeDistributions();

  // Split distributions with the same prior
  muLogMacro(<< "Splitting distributions with same prior\n");
  for (unsigned int iprior = 0; iprior < numPriors; iprior++)
  {
    if (m_NumberOfGaussians[iprior] > 1)
    {
      if (m_DoMSTSplit)
        this->SplitPriorMST(iprior);
      else
        this->SplitPrior(iprior);
    }
  }

  // Update posteriors using the initial distribution parameters
  this->ComputePosteriors(false);

  //double logLikelihood = vnl_huge_val(1.0);
  double logLikelihood = 1.0 / vnl_math::eps;
  double deltaLogLikelihood = 1.0;

  unsigned int biasdegree = 0;

  // EM loop
  bool converged = false;
  unsigned int iter = 0;
  while (!converged)
  {

    iter++;

    muLogMacro(<< "\n\nEM iteration " << iter << "\n");
    muLogMacro(<< "---------------------\n");

    // Fix the mean for the first class (speeds up convergence)
    // Makes sure that image intensity in the same range after bias correction
    MatrixType offsetMeans(numChannels, numClasses, 0);
    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    {
      for (unsigned int iclass = 0; iclass < (numClasses-1); iclass++)
        offsetMeans(ichan, iclass) =
          -1.0 * (m_Means(ichan, 0) - refMean(ichan));
    }
    m_Means += offsetMeans;

    for (unsigned int iclass = 0; iclass < numClasses; iclass++)
    {
      muLogMacro(<< "Class " << (iclass+1) << " mean: ");
      for (unsigned int ichan = 0; ichan < numChannels; ichan++)
        muLogMacro(<< m_Means(ichan, iclass) << "\t");
      muLogMacro(<< "\n");
    }
/*
    muLogMacro(<< "\n");
    for (unsigned int iclass = 0; iclass < numClasses; iclass++)
    {
      muLogMacro(<< "Class " << (iclass+1) << " covariance:\n");
      muLogMacro(<< m_Covariances[iclass]);
    }
*/

    // Bias correction
    if (m_MaxBiasDegree > 0)
    {
      if ((deltaLogLikelihood < m_BiasLikelihoodTolerance)
          &&
          (biasdegree < m_MaxBiasDegree))
      {
        biasdegree++;
      }
      this->CorrectBias(biasdegree, false);
    }

    // Update distribution parameters after bias correction
    this->ComputeDistributions();

    // Recompute posteriors, not at full resolution
    this->ComputePosteriors(false);

    double prevLogLikelihood = logLikelihood;
    if (prevLogLikelihood == 0)
      prevLogLikelihood = vnl_math::eps;

    // Compute log-likelihood and normalize posteriors
    logLikelihood = 0;
    for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
      for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
        for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
        {
          if (m_Mask->GetPixel(ind) == 0)
            continue;

          double tmp = vnl_math::eps;
          for (unsigned int iclass = 0; iclass < numClasses; iclass++)
            tmp += m_Posteriors[iclass]->GetPixel(ind);
          logLikelihood += log(tmp);

          // Normalize posteriors
          // TODO: normalize to min-max range of prob pixel type?
          // Assume that probability pixel type is float / double for now
          for (unsigned int iclass = 0; iclass < numClasses; iclass++)
            m_Posteriors[iclass]->SetPixel(ind, (ProbabilityImagePixelType)
              (m_Posteriors[iclass]->GetPixel(ind) / tmp));
        }
    muLogMacro(<< "log(likelihood) = " << logLikelihood << "\n");

// TODO: move to before prevL update
    deltaLogLikelihood =
      vnl_math_abs((logLikelihood - prevLogLikelihood) / prevLogLikelihood);
      //(logLikelihood - prevLogLikelihood) / vnl_math_abs(prevLogLikelihood);

    muLogMacro(<< "delta log(likelihood) = " << deltaLogLikelihood << "\n");

    // Convergence check
    converged =
      (iter > m_MaximumIterations)
      ||
      (deltaLogLikelihood < 0)
      ||
      ((deltaLogLikelihood < m_LikelihoodTolerance)
        &&
        (biasdegree == m_MaxBiasDegree));

  } // end EM loop

  muLogMacro(<< "Done computing posteriors with " << iter << " iterations\n");

  // Bias correction at full resolution, still using downsampled images
  // for computing the bias field coeficients
  if (m_MaxBiasDegree > 0)
  {
    this->CorrectBias(biasdegree, true);
    this->ComputeDistributions();
  }

  // Recompute posteriors at full resolution
  this->ComputePosteriors(true);

  // Normalize posteriors at full resolution
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;

        double tmp = vnl_math::eps;
        for (unsigned int iclass = 0; iclass < numClasses; iclass++)
          tmp += m_Posteriors[iclass]->GetPixel(ind);

        for (unsigned int iclass = 0; iclass < numClasses; iclass++)
          m_Posteriors[iclass]->SetPixel(ind, (ProbabilityImagePixelType)
            (m_Posteriors[iclass]->GetPixel(ind) / tmp));
      }

}

// Labeling using maximum a posteriori, also do brain stripping using
// mathematical morphology and connected component
template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::ComputeLabels()
{

  itkDebugMacro(<< "ComputeLabels");

  unsigned int numPriors = m_Priors.GetSize();

  unsigned int numClasses = 0;
  for (unsigned int i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  unsigned int numFGClasses = numClasses - m_NumberOfGaussians[numPriors-1];

  InputImageRegionType region = m_InputImages[0]->GetLargestPossibleRegion();

  m_Labels = ByteImageType::New();
  m_Labels->SetRegions(region);
  m_Labels->Allocate();
  m_Labels->SetOrigin(m_InputImages[0]->GetOrigin());
  m_Labels->SetSpacing(m_InputImages[0]->GetSpacing());

  ByteImagePointer mask = ByteImageType::New();
  mask->SetRegions(region);
  mask->Allocate();
  mask->SetOrigin(m_InputImages[0]->GetOrigin());
  mask->SetSpacing(m_InputImages[0]->GetSpacing());

  MersenneTwisterRNG rng;

  ByteImageIndexType ind;
  ByteImageSizeType size = m_Labels->GetLargestPossibleRegion().GetSize();

  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;

        double maxv = m_Posteriors[0]->GetPixel(ind);
        unsigned int imax = 0;

        for (unsigned int iclass = 1; iclass < numClasses; iclass++)
        {
          double v = m_Posteriors[iclass]->GetPixel(ind);
          if (v > maxv)
          {
            maxv = v;
            imax = iclass;
          }
        }

        // Choose label at random if there's a tie
        // Only for brain classes
        DynArray<unsigned int> tieIndices;
        tieIndices.Allocate(numFGClasses);
        tieIndices.Append(imax);

        for (unsigned int iclass = 0; iclass < numFGClasses; iclass++)
        {
          if (imax == iclass)
            continue;
          double v = m_Posteriors[iclass]->GetPixel(ind);
          if (vnl_math_abs(v - maxv) < 1e-20)
            tieIndices.Append(iclass);
        }

        if (tieIndices.GetSize() > 1)
        {
          unsigned int whichIndex =
            (unsigned int)
              rng.GenerateUniformIntegerUpToK(tieIndices.GetSize() - 1);
          imax = tieIndices[whichIndex];
          maxv = m_Posteriors[imax]->GetPixel(ind);
        }

        unsigned char label = 0;
        unsigned char fgflag = 0;

        // Only use non-zero probabilities and foreground classes
        if (maxv > 0 && imax < numFGClasses)
        {
          label = (unsigned char)(imax+1);
          fgflag = 1;
        }

        m_Labels->SetPixel(ind, label);
        mask->SetPixel(ind, fgflag);

      }

  // Binary opening
  typedef itk::BinaryBallStructuringElement<unsigned char, 3> StructElementType;
  typedef
    itk::BinaryDilateImageFilter<ByteImageType, ByteImageType,
      StructElementType> DilateType;
  typedef
    itk::BinaryErodeImageFilter<ByteImageType, ByteImageType,
      StructElementType> ErodeType;

  StructElementType structel;
  structel.SetRadius(2);
  structel.CreateStructuringElement();

  typename ErodeType::Pointer erode = ErodeType::New();
  erode->SetErodeValue(1);
  erode->SetInput(mask);
  erode->SetKernel(structel);

  erode->Update();

  typename DilateType::Pointer dil = DilateType::New();
  dil->SetDilateValue(1);
  dil->SetInput(erode->GetOutput());
  dil->SetKernel(structel);

  dil->Update();

  // Find largest cluster
  typedef itk::Image<unsigned short, 3> ShortImageType;
  typedef ShortImageType::Pointer ShortImagePointer;

  typedef ConnectedComponentsFilter<ByteImageType, ShortImageType> CCType;
  typename CCType::Pointer ccfilter = CCType::New();

  ccfilter->SetInput(dil->GetOutput());
  ccfilter->Update();

  ShortImagePointer components = ccfilter->GetOutput();

  // Remove floating clusters
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        if (components->GetPixel(ind) != 1)
          m_Labels->SetPixel(ind, 0);
      }

}

template <class TInputImage, class TProbabilityImage>
void
EMSegmentationFilter <TInputImage, TProbabilityImage>
::CleanUp()
{

  itkDebugMacro(<< "CleanUp");

  unsigned int numChannels = m_InputImages.GetSize();

  unsigned int numPriors = m_Priors.GetSize();

  unsigned int numClasses = 0;
  for (unsigned int i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  unsigned int numFGClasses = numClasses - m_NumberOfGaussians[numPriors-1];

  InputImageIndexType ind;

  InputImageSizeType size =
    m_InputImages[0]->GetLargestPossibleRegion().GetSize();

  // Make sure all output corrected image intensities are within a "nice" range
  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
  {
    InputImagePointer img = m_CorrectedImages[ichan];

    double minv = vnl_huge_val(1.0f);
    double maxv = -vnl_huge_val(1.0f);
    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          if (m_Labels->GetPixel(ind) == 0)
            continue;

          double v = img->GetPixel(ind);
          if (v < minv)
            minv = v;
          if (v > maxv)
            maxv = v;
        }

    double range = maxv - minv;

    double outmax = itk::NumericTraits<short>::max();

    // Renormalize to [0, outmax]
    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          double v = img->GetPixel(ind);
          if (v < minv)
            v = minv;
          if (v > maxv)
            v = maxv;
          v = (v - minv) / range * outmax;
          v = floor(v + 0.5);
          img->SetPixel(ind, v);
        }
  }

  // Zero the foreground posteriors with label 0
  for (unsigned int iclass = 0; iclass < numFGClasses; iclass++)
  {
    ProbabilityImagePointer post = m_Posteriors[iclass];

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          if (m_Labels->GetPixel(ind) == 0)
            post->SetPixel(ind, 0);
        }
  }

  // Normalize posteriors
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        double tmp = vnl_math::eps;
        for (unsigned int iclass = 0; iclass < numClasses; iclass++)
          tmp += m_Posteriors[iclass]->GetPixel(ind);

        for (unsigned int iclass = 0; iclass < numClasses; iclass++)
          m_Posteriors[iclass]->SetPixel(ind, (ProbabilityImagePixelType)
            (m_Posteriors[iclass]->GetPixel(ind) / tmp));
      }

}

#endif
