
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"

#include "vnl/algo/vnl_qr.h"
#include "vnl/vnl_math.h"

// Image processing
#include "ConnectedComponentsFilter.h"
#include "LLSBiasCorrector.h"
#include "NeonateSegmentationFilter.h"
//#include "UnsharpMaskingImageFilter.h"

// Robust statistics / clustering
#include "FastMCDSampleFilter.h"
#include "MSTClusteringProcess.h"

// SPR
#include "KernelWidthEstimator.h"
#include "KNNClassifier.h"
#include "ReducedSetDensityEstimator.h"

#include "Log.h"
#include "MersenneTwisterRNG.h"

#include "muMacros.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.1415926535898
#endif

#define MCD_REFONLY 0
#define DROP_EXTREME 1

#define MIN_CLASS_SAMPLES 10

////////////////////////////////////////////////////////////////////////////////
// Class definitions
////////////////////////////////////////////////////////////////////////////////

NeonateSegmentationFilter
::NeonateSegmentationFilter()
{

  m_Mask = 0;

  m_LabelImage = 0;

  m_GImage = 0;

  m_SampleSpacing = 2.0;

  // Bias
  m_MaxBiasDegree = 4;
  m_BiasLikelihoodTolerance = 1e-3;

  // EM convergence parameters
  m_LikelihoodTolerance = 1e-5;
  m_MaximumIterations = 30;

  m_InputModified = false;

  m_PriorWeights = VectorType(0);
  m_NumberOfGaussians = 0;

  m_PriorLookupTable = 0;

  m_MahalanobisThreshold = 3.0;
  m_PriorThresholdFraction = 0.9;
  m_KernelWidthFraction = 0.05;

  m_ReferenceImageIndex = 0;
  m_ReferenceModality = T2;

  m_MyelinPriorWeight = 1.0;

  m_Prototypes.Clear();
  m_PrototypeLabels.Clear();

  m_RescaledMax = itk::NumericTraits<short>::max();

  m_MaxSamples = 20000;

  m_DoMSTSplit = true;

  m_FOVMask = 0;

  m_DoneEM = false;
  m_DoneNPPrototypes = false;

}

NeonateSegmentationFilter
::~NeonateSegmentationFilter()
{

  if (m_PriorLookupTable != 0)
    delete [] m_PriorLookupTable;

}

void
NeonateSegmentationFilter
::CheckInput()
{

  if (m_MaximumIterations == 0)
  {
    itkExceptionMacro(<< "Maximum iterations set to zero");
  }

  if (m_InputImages.GetSize() == 0)
  {
    itkExceptionMacro(<< "No input images");
  }

  if (m_Priors.GetSize() < 1)
  {
    itkExceptionMacro(<< "Must have one or more class probabilities");
  }


  itkDebugMacro(
    << "CheckInput: " << m_InputImages.GetSize() << " images, "
    << m_Priors.GetSize() << " priors");

  InputImageSizeType size =
    m_InputImages[0]->GetLargestPossibleRegion().GetSize();

  for (unsigned i = 1; i < m_InputImages.GetSize(); i++)
  {
    if (m_InputImages[i]->GetImageDimension() != 3)
      itkExceptionMacro(
        << "InputImage [" << i << "] has invalid dimension: only supports 3D images");
    InputImageSizeType isize =
      m_InputImages[i]->GetLargestPossibleRegion().GetSize();
    if (size != isize)
      itkExceptionMacro(
        << "Image data 3D size inconsistent: " << size << " vs " << isize);
  }

  for (unsigned i = 0; i < m_Priors.GetSize(); i++)
  {
    if (m_Priors[i]->GetImageDimension() != 3)
      itkExceptionMacro(
        << "Prior [" << i << "] has invalid dimension: only supports 3D images");
    ProbabilityImageSizeType psize =
      m_Priors[i]->GetLargestPossibleRegion().GetSize();
    if (size != psize)
      itkExceptionMacro(
        << "Image data and prior 3D size inconsistent: " << size << " vs " << psize);
  }

}

void
NeonateSegmentationFilter
::SetInputImages(DynArray<InputImagePointer> data) {

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

  // Initial fill of corrected images using raw input
  InputImageIndexType ind;
  InputImageSizeType size =
    m_InputImages[0]->GetLargestPossibleRegion().GetSize();
  for (unsigned int i = 0; i < m_InputImages.GetSize(); i++)
  {
    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          m_CorrectedImages[i]->SetPixel(ind, m_InputImages[i]->GetPixel(ind));
        }
  }

  m_InputModified = true;

}

void
NeonateSegmentationFilter
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
        double sumPrior = DBL_EPSILON;
        for (unsigned iprior = 0; iprior < numPriors; iprior++)
          sumPrior += m_Priors[iprior]->GetPixel(ind);
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
  m_NumberOfGaussians[0] = 2;
  m_NumberOfGaussians[numPriors-1] = 3;

  m_InputModified = true;

}

void
NeonateSegmentationFilter
::SetPriorWeights(const VectorType& w)
{

  itkDebugMacro(<< "SetPriorWeights");

  if (w.size() != m_Priors.GetSize())
    itkExceptionMacro(<< "Number of prior weights invalid");

  for (unsigned i = 0; i < w.size() ; i++)
  {
    if (w[i] == 0.0)
      itkExceptionMacro(<< "Prior weight " << i << " is zero");
  }

  m_PriorWeights = w;

  m_InputModified = true;

}

void
NeonateSegmentationFilter
::SetNumberOfGaussians(const unsigned int* c)
{

  itkDebugMacro(<< "SetNumberOfGaussians");

  if (c == NULL)
    itkExceptionMacro(<< "Number of cluster info invalid");

  for (unsigned i = 0; i < m_Priors.GetSize(); i++)
    m_NumberOfGaussians[i] = c[i];

  m_InputModified = true;

}

void
NeonateSegmentationFilter
::SplitPriorMST(unsigned int iprior)
{

  unsigned int numChannels = m_InputImages.GetSize();
  unsigned int numPriors = m_Priors.GetSize();

  if (iprior >= numPriors)
    itkExceptionMacro(<< "Invalid prior index: " << iprior+1);

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
  skips[0] = (unsigned)floor(4.0 / spacing[0]);
  skips[1] = (unsigned)floor(4.0 / spacing[1]);
  skips[2] = (unsigned)floor(4.0 / spacing[2]);
  //skips[0] = (unsigned)floor(m_SampleSpacing / spacing[0]);
  //skips[1] = (unsigned)floor(m_SampleSpacing / spacing[1]);
  //skips[2] = (unsigned)floor(m_SampleSpacing / spacing[2]);

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
        if (probImg->GetPixel(ind) > maxP)
          maxP = probImg->GetPixel(ind);
      }

  // Select samples by thresholding prior with value above tau
  double tau = m_PriorThresholdFraction * maxP;

  unsigned int numSamples = 0;
  DynArray<MSTClusteringProcess::VertexType> samples;

  {
    MatrixType tmpmat =
      this->ProbabilitySampling(probImg, tau, vnl_huge_val(1.0), 4.0);
    numSamples = tmpmat.rows();

    samples.Allocate(numSamples);
    for (unsigned int r = 0; r < numSamples; r++)
    {
      MSTClusteringProcess::VertexType x(numChannels);
      for (unsigned int ichan = 0; ichan < numChannels; ichan++)
        x[ichan] = tmpmat(r, ichan);
      samples.Append(x);
    }
  }

  MSTClusteringProcess mstProc;
  mstProc.SetInputVertices(samples);
  mstProc.SortOn();

  DynArray<unsigned int> clusterMap;
  clusterMap.Initialize(numSamples, 0);

  unsigned int lastNumClusters = 0;

  double T;
  for (T = 2.5; T >= 0.999; T -= 0.01)
  {
    muLogMacro(<< "  MST clustering, T = " << T << "\n");

    unsigned int numClusters = mstProc.GetClusters(clusterMap.GetRawArray(), T);

    if (numClusters < m_NumberOfGaussians[iprior])
      continue;

    if (numClusters == lastNumClusters)
      continue;

    lastNumClusters = numClusters;

    // Check cluster sizes
    bool sizeOK = true;
    for (unsigned int m = 0; m < m_NumberOfGaussians[iprior]; m++)
    {
      unsigned int count = 0;
      for (unsigned int i = 0; i < numSamples; i++)
        if (clusterMap[i] == m)
          count++;
      if (count < MIN_CLASS_SAMPLES)
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
        v /= (count-1);

        if (r == c)
          v += 1e-10;

        covtmp(r, c) = v;
        covtmp(c, r) = v;
      }
    }

    //covtmp *= 1.1;

    m_Covariances[iclass] = covtmp;

  }

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
        << "  Replacing background mean " << m_Means.get_column(imin)
        << " with zeros\n");
      VectorType v = m_Means.get_column(numClasses-1);
      m_Means.set_column(imin, v);
    }
    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
      m_Means(ichan, numClasses-1) = 0;
  }

}

void
NeonateSegmentationFilter
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

  // Scale the parameters
  for (unsigned int k = 0; k < m_NumberOfGaussians[iprior]; k++)
  {
    double s = 1.0 / pow(1.5, (double)k);
    unsigned int iclass = istart+k;
    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
      m_Means(ichan, iclass) = s * m_Means(ichan, iclass);
    //MatrixType cov = m_Covariances[iclass];
    //cov *= s;
    //m_Covariances[iclass] = cov;
  }

  // Set mean of last class to zero
  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    m_Means(ichan, numClasses-1) = 0;

}

NeonateSegmentationFilter::ByteImagePointer
NeonateSegmentationFilter
::GetOutput()
{

  this->Update();

  return m_LabelImage;

}

DynArray<NeonateSegmentationFilter::ByteImagePointer>
NeonateSegmentationFilter
::GetBytePosteriors()
{
  this->Update();

  unsigned int numClasses = m_Posteriors.GetSize();

  DynArray<ByteImagePointer> bytePosts;
  bytePosts.Allocate(numClasses);

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
          double p = floor(255.0*post->GetPixel(ind) + 0.5);
          if (p < 0)
            p = 0;
          if (p > 255)
            p = 255;
          tmp->SetPixel(ind, (unsigned char)p);
        }

    bytePosts.Append(tmp);
  }

  return bytePosts;

}

DynArray<NeonateSegmentationFilter::ProbabilityImagePointer>
NeonateSegmentationFilter
::GetPosteriors()
{

  this->Update();

  return m_Posteriors;

}

DynArray<NeonateSegmentationFilter::ProbabilityImagePointer>
NeonateSegmentationFilter
::GetCorrected()
{

  this->Update();

  return m_CorrectedImages;

}

void
NeonateSegmentationFilter
::Update()
{

  if (!m_InputModified)
    return;

  this->CheckInput();

  this->ComputeInitialMask();

  this->ComputePriorLookupTable();

  this->ComputeGImage();

  this->EMLoop();

  //this->RescaleCorrected();

  // Recompute G image after EM bias corr
  this->ComputeGImage();

  m_DoneEM = true;

  m_InputModified = false;

}

void
NeonateSegmentationFilter
::DoKNearestNeighbors()
{
  // Make sure we've done an EM loop
  if (!m_DoneEM)
    this->Update();

  itkDebugMacro(<< "DoKNearest");

  // Make sure we have the non-parametric prototypes
  if (!m_DoneNPPrototypes)
  {
    this->ComputePrototypesMSTFull();
    m_DoneNPPrototypes = true;
  }

  unsigned int numChannels = m_InputImages.GetSize();

  unsigned int numPrototypes = m_Prototypes.GetSize();

  MatrixType tmp(numPrototypes, numChannels);
  for (unsigned int i = 0; i < numPrototypes; i++)
  {
    tmp.set_row(i, m_Prototypes[i]);
  }

  KNNClassifier knn;
  knn.SetDimension(numChannels);
  knn.SetKNeighbors(5);

  knn.SetTrainingData(tmp, m_PrototypeLabels);
  tmp = MatrixType(0, numChannels);

  muLogMacro(<< "Condensing knn prototypes...\n");
  knn.Update();

  ByteImageIndexType ind;
  ByteImageSizeType size =
    m_Priors[0]->GetLargestPossibleRegion().GetSize();

  muLogMacro(<< "K-nearest neighbors loop...\n");
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;

        KNNClassifier::VectorType x(numChannels);
        for (unsigned int ichan = 0; ichan < numChannels; ichan++)
        {
          x[ichan] = m_CorrectedImages[ichan]->GetPixel(ind);
        }

        m_LabelImage->SetPixel(ind, knn.Classify(x));
      }

  // TODO: clean up KNN labels?
  // Later, KNN not used right now
}

void
NeonateSegmentationFilter
::DoKernelDensityEstimation()
{
  if (!m_DoneEM)
    this->Update();

  if (!m_DoneNPPrototypes)
  {
    this->ComputePrototypesMSTFull();
    m_DoneNPPrototypes = true;
  }

  unsigned int numChannels = m_InputImages.GetSize();
  unsigned int numPriors = m_Priors.GetSize();

  unsigned int numPrototypes = m_Prototypes.GetSize();

  // Kernel width
  double lambda0 = m_KernelWidthFraction * m_RescaledMax;
  double d_lambda = 1e-4 * m_RescaledMax;

  muLogMacro(<< "Initial guess for kernel width = " << lambda0 << "\n");

  unsigned int numFGClasses = 0;
  for (unsigned int i = 0; i < (numPriors-1); i++)
    numFGClasses += m_NumberOfGaussians[i];

  unsigned int numClasses = numFGClasses + m_NumberOfGaussians[numPriors-1];

  // Do background? Assume removed by previous EM seg phase
  //numFGClasses += 1;

  muLogMacro(
    << "Kernel density estimation with " << numFGClasses << " classes\n");

  // Build list of density estimators
  ReducedSetDensityEstimator* rsdeList = 
    new ReducedSetDensityEstimator[numFGClasses];

  for (unsigned int iclass = 0; iclass < numFGClasses; iclass++)
  {
    unsigned char label = iclass+1;

    unsigned int n = 0;
    for (unsigned int i = 0; i < numPrototypes; i++)
      if (m_PrototypeLabels[i] == label)
        n++;

    if (n == 0)
      itkExceptionMacro(<< "No voxels labeled " << label << "\n");

    MatrixType tmp(n, numChannels);

    unsigned int r = 0;
    for (unsigned int i = 0; i < numPrototypes; i++)
      if (m_PrototypeLabels[i] == label)
      {
        tmp.set_row(r, m_Prototypes[i]);
        r++;
      }

/*
    MatrixType cov = m_Covariances[iclass];
    double avgStd = 0.0;
    for (int i = 0; i < numChannels; i++)
      avgStd += cov(i, i);
    avgStd /= numChannels;

    muLogMacro(<< "Average eigvalues of covariance = " << avgStd << "\n");
*/

    KernelWidthEstimator bandwidthEstimator;
    bandwidthEstimator.SetMaximumIterations(200);
    bandwidthEstimator.SetNumberOfTestPoints(n/5 + 10);
    bandwidthEstimator.SetLearningRate(0.1);
    bandwidthEstimator.SetStepLength(d_lambda);
    bandwidthEstimator.SetGaussianParameters(
      m_Means.get_column(iclass), m_Covariances[iclass]);

    double classBandwidth = bandwidthEstimator.GetKernelWidth(tmp, lambda0);

    muLogMacro(
      << "Class " << iclass+1 << " estimated kernel width "
      << classBandwidth << "\n");

    rsdeList[iclass].SetDimension(numChannels);
    rsdeList[iclass].SetInputSet(tmp);
//PP
    //rsdeList[iclass].SetKernelWidth(lambda0);
    rsdeList[iclass].SetKernelWidth(classBandwidth);

    if (n > 10)
    {
      muLogMacro(<< "RSDE update for class " << (unsigned int)label << "...\n");
      rsdeList[iclass].Update();
    }
    else
    {
      muLogMacro(
        << "Too few samples for class " << (unsigned int)label
        << ": " << n << ", not using RSDE\n");
    }
  }

  ByteImageIndexType ind;
  ByteImageSizeType size =
    m_Priors[0]->GetLargestPossibleRegion().GetSize();

  muLogMacro(<< "Kernel density estimation for each voxel...\n");
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {

        if (m_Mask->GetPixel(ind) == 0)
        {
          for (unsigned int i = 0; i < numFGClasses; i++)
            m_Posteriors[i]->SetPixel(ind, 0);
          double sump = DBL_EPSILON;
          for (unsigned int i = numFGClasses; i < numClasses; i++)
            sump += m_Posteriors[i]->GetPixel(ind);
          for (unsigned int i = numFGClasses; i < numClasses; i++)
            m_Posteriors[i]->SetPixel(ind, m_Posteriors[i]->GetPixel(ind)/sump);
          continue;
        }

        ReducedSetDensityEstimator::VectorType x(numChannels);
        for (unsigned int ichan = 0; ichan < numChannels; ichan++)
          x[ichan] = m_CorrectedImages[ichan]->GetPixel(ind);

        DynArray<double> probs;
        probs.Initialize(numFGClasses, 0);

        for (unsigned int c = 0; c < numFGClasses; c++)
        {
          if (rsdeList[c].GetNumberOfTrainingPoints() > 10)
            probs[c] = rsdeList[c].Evaluate(x);
          else
            probs[c] = rsdeList[c].EvaluateWithoutReduce(x);
          // Brain pdf should be non-zero
          probs[c] += 1e-10;
        }

        // Modulate with prior probabilities
        for (unsigned int i = 0; i < numFGClasses; i++)
        {
          unsigned int iprior = m_PriorLookupTable[i];
          probs[i] = probs[i] *
            m_Priors[iprior]->GetPixel(ind)
            /
            m_NumberOfGaussians[iprior];
          if (i > 0)
            probs[i] *= m_PriorWeights[iprior];
          else
            probs[i] *= m_MyelinPriorWeight;
        }

        // Rescale and update posterior probability values
        double sump = DBL_EPSILON;
        for (unsigned int i = 0; i < numFGClasses; i++)
          sump += probs[i];
        for (unsigned int i = 0; i < numFGClasses; i++)
          m_Posteriors[i]->SetPixel(ind, probs[i] / sump);

        // Zero the background probabilities, assume it's already handled
        // through the mask from EM
        for (unsigned int i = numFGClasses; i < numClasses; i++)
          m_Posteriors[i]->SetPixel(ind, 0);

      }

  delete [] rsdeList;

  this->ComputeLabelImage();

}

void
NeonateSegmentationFilter
::ComputePriorLookupTable()
{

  itkDebugMacro(<< "ComputePriorLookupTable");

  if (m_PriorLookupTable != 0)
    delete [] m_PriorLookupTable;

  unsigned int numPriors = m_Priors.GetSize();

  unsigned numClasses = 0;
  for (unsigned i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

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

void
NeonateSegmentationFilter
::ComputeInitialMask()
{

  itkDebugMacro(<< "ComputeInitialMask");

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
        for (unsigned iprior = 0; iprior < (numPriors-1); iprior++)
          tmp += m_Priors[iprior]->GetPixel(ind);
        if (tmp > 0)
          m_Mask->SetPixel(ind, 1);
        else
          m_Mask->SetPixel(ind, 0);
      }

  // Mask out voxels outside FOV, if specified
  if (!m_FOVMask.IsNull())
  {
    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          if (m_FOVMask->GetPixel(ind) == 0)
            m_Mask->SetPixel(ind, 0);
        }
  }

#if 0
  // Set a minimum floor probability for brain pixels inside mask
  for (unsigned int iprior = 0; iprior < (numPriors-1); iprior++)
  {
    ProbabilityImagePointer prob = m_Priors[iprior];

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          if (m_Mask->GetPixel(ind) == 0)
            continue;
          prob->SetPixel(ind, prob->GetPixel(ind) + 1e-10);
        }
  }
#endif

}

void
NeonateSegmentationFilter
::ComputeGImage()
{

  itkDebugMacro(<< "ComputeGImage");

  m_GImage = InputImageType::New();
  m_GImage->SetRegions(m_CorrectedImages[0]->GetLargestPossibleRegion());
  m_GImage->Allocate();
  m_GImage->SetOrigin(m_CorrectedImages[0]->GetOrigin());
  m_GImage->SetSpacing(m_CorrectedImages[0]->GetSpacing());

  InputImageIndexType ind;

  InputImageSizeType size =
    m_CorrectedImages[0]->GetLargestPossibleRegion().GetSize();

  // Clear G
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        m_GImage->SetPixel(ind, 0);

  //typedef UnsharpMaskingImageFilter<InputImageType> UnsharpMaskingType;
  //UnsharpMaskingType::Pointer sharpFilter = UnsharpMaskingType::New();

  //sharpFilter->SetBlurVariance(4.0);

  typedef itk::DiscreteGaussianImageFilter<InputImageType, InputImageType>
    BlurType;
  BlurType::Pointer blurFilter = BlurType::New();

  blurFilter->SetVariance(1.0);

  typedef itk::GradientMagnitudeImageFilter<InputImageType, InputImageType>
    GradMagFilterType;
  GradMagFilterType::Pointer gradMagFilter = GradMagFilterType::New();

  // Sum of squares
  for (unsigned int k = 0; k < m_CorrectedImages.GetSize(); k++)
  {
    blurFilter->SetInput(m_CorrectedImages[k]);
    blurFilter->Update();

    //sharpFilter->SetInput(m_CorrectedImages[k]);
    //sharpFilter->SetInput(blurFilter->GetOutput());
    //sharpFilter->Update();

    gradMagFilter->SetInput(blurFilter->GetOutput());
    //gradMagFilter->SetInput(sharpFilter->GetOutput());
    gradMagFilter->Update();

    InputImageType::Pointer tmp = gradMagFilter->GetOutput();

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          double g = tmp->GetPixel(ind);
          m_GImage->SetPixel(ind, m_GImage->GetPixel(ind) + g*g);
        }
  }

  // Square root for 2-norm
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        m_GImage->SetPixel(ind, sqrt(m_GImage->GetPixel(ind)));
      }

  // TODO: grayscale opening
  // no need for trimmed /clamped g values

}

void
NeonateSegmentationFilter
::ComputeGaussians()
{

  itkDebugMacro(<< "ComputeGaussians");

  unsigned numChannels = m_InputImages.GetSize();
  unsigned numPriors = m_Priors.GetSize();

  unsigned numClasses = 0;
  for (unsigned i = 0; i < numPriors; i++)
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

  // Compute sum of posteriors for each class
  muLogMacro(<< "Computing sum of posteriors for each class...\n");
  VectorType sumClassProb(numClasses);
  for (unsigned iclass = 0; iclass < numClasses; iclass++)
  {
    double tmp = DBL_EPSILON;
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
  muLogMacro(<< "Computing means...\n");
  m_Means = MatrixType(numChannels, numClasses, 0.0);
  for (unsigned int iclass = 0; iclass < numClasses; iclass++)
  {
    if (sumClassProb[iclass] < 1e-20)
      continue;

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

  // Fix mean of last class to zero (background)
  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    m_Means(ichan, numClasses-1) = 0;

  // Compute covariances
  muLogMacro(<< "Computing covariances...\n");
  DynArray<MatrixType> oldCovariances = m_Covariances;
  m_Covariances.Clear();
  m_Covariances.Allocate(numClasses);
  for (unsigned int iclass = 0; iclass < numClasses; iclass++)
  {

    if (sumClassProb[iclass] < 1e-20)
      continue;

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
              diff1 = img1->GetPixel(ind) - mu1;
              diff2 = img2->GetPixel(ind) - mu2;
              v += m_Posteriors[iclass]->GetPixel(ind) * (diff1*diff2);
            }

        v /= sumClassProb[iclass];

        if (r == c)
          v += 1e-10;

        // Assign value to the covariance matrix (symmetric)
        covtmp(r, c) = v;
        covtmp(c, r) = v;
      }
    }

    vnl_qr<float> qr(covtmp);
    double detcov = qr.determinant();

    if (detcov <= 0)
      m_Covariances.Append(oldCovariances[iclass]);
    else
      m_Covariances.Append(covtmp);

  } // end covariance loop
}

void
NeonateSegmentationFilter
::ComputePosteriors(bool fullRes)
{

  itkDebugMacro(<< "ComputePosteriors");

  unsigned numChannels = m_InputImages.GetSize();
  unsigned numPriors = m_Priors.GetSize();

  unsigned numClasses = 0;
  for (unsigned i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  if (fullRes)
    itkDebugMacro(<< "Computing posteriors at full resolution");

  ProbabilityImageIndexType ind;

  ProbabilityImageSizeType size =
    m_Priors[0]->GetLargestPossibleRegion().GetSize();

  ProbabilityImageSpacingType spacing = m_Priors[0]->GetSpacing();

  ProbabilityImageOffsetType skips;
  skips[0] = (unsigned)floor(m_SampleSpacing / spacing[0]);
  skips[1] = (unsigned)floor(m_SampleSpacing / spacing[1]);
  skips[2] = (unsigned)floor(m_SampleSpacing / spacing[2]);

  if (skips[0] == 0 || fullRes)
    skips[0] = 1;
  if (skips[1] == 0 || fullRes)
    skips[1] = 1;
  if (skips[2] == 0 || fullRes)
    skips[2] = 1;

  for (unsigned iclass = 0; iclass < numClasses; iclass++)
  {

    ProbabilityImagePointer post = m_Posteriors[iclass];
    ProbabilityImagePointer prior = m_Priors[m_PriorLookupTable[iclass]];

    double priorScale = 1.0;

    if (iclass > 0)
    {
      priorScale = 
        m_PriorWeights[m_PriorLookupTable[iclass]]
        /
        m_NumberOfGaussians[m_PriorLookupTable[iclass]];
    }
    else
    {
      priorScale = 
        m_MyelinPriorWeight
        /
        m_NumberOfGaussians[m_PriorLookupTable[iclass]];

    }

    vnl_qr<float> qr(m_Covariances[iclass]);
    double detcov = qr.determinant();

    if (detcov < 0)
      itkExceptionMacro(<< "Detcov class " << iclass << " <= 0 : " << detcov << "\n" << m_Covariances[iclass]);

    double denom =
      pow(2*M_PI, m_InputImages.GetSize()/2.0) * sqrt(detcov) + DBL_EPSILON;

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

void
NeonateSegmentationFilter
::RescaleCorrected()
{

  unsigned int numChannels = m_CorrectedImages.GetSize();
  InputImageIndexType ind;
  InputImageSizeType size =
    m_CorrectedImages[0]->GetLargestPossibleRegion().GetSize();

  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
  {
    InputImagePointer img = m_CorrectedImages[ichan];

    double maxI = -vnl_huge_val(1.0);
    double minI = vnl_huge_val(1.0);

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          if (m_LabelImage->GetPixel(ind) == 0)
            continue;
          double v = img->GetPixel(ind);
          if (v > maxI)
            maxI = v;
          if (v < minI)
            minI = v;
        }

    double rangeI = maxI - minI;

    if (rangeI < 1e-10)
      rangeI = 1e-10;

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          double v = img->GetPixel(ind);

          // Clamp to brain region intensity range
          if (v < minI)
            v = minI;
          if (v > maxI)
            v = maxI;

          // Rescale
          v = (v - minI) / rangeI * m_RescaledMax + 1;

          // Round
          v = floor(v + 0.5);

          img->SetPixel(ind, v);
        }
  }

}

void
NeonateSegmentationFilter
::CorrectBias(unsigned int degree, bool fullRes)
{

  muLogMacro(<< "Bias correction, max degree = " << degree << "\n");

  unsigned int numPriors = m_Priors.GetSize();

  unsigned int numClasses = 0;
  for (unsigned i = 0; i < numPriors; i++)
    numClasses += m_NumberOfGaussians[i];

  // Perform bias correction
  typedef LLSBiasCorrector<InputImageType, ProbabilityImageType>
    BiasCorrectorType;
  typedef BiasCorrectorType::Pointer BiasCorrectorPointer;

  BiasCorrectorPointer biascorr = BiasCorrectorType::New();

  //biascorr->DebugOn();

  DynArray<ProbabilityImagePointer> biasPosteriors;
  for (unsigned int j = 0; j < (numClasses-1); j++)
    biasPosteriors.Append(m_Posteriors[j]);

  biascorr->SetSampleSpacing(2.0*m_SampleSpacing);
  biascorr->SetWorkingSpacing(m_SampleSpacing);
  biascorr->SetMaxDegree(degree);

  biascorr->SetMask(m_Mask);
  biascorr->SetProbabilities(biasPosteriors);

  biascorr->SetClampBias(false);
  biascorr->SetMaximumBiasMagnitude(5.0);

  if (this->GetDebug())
    biascorr->DebugOn();

  for (unsigned int ichan = 0; ichan < m_InputImages.GetSize(); ichan++)
  {
    itkDebugMacro(<< "Correcting image " << ichan+1);

// Use the distribution of log-intensities instead
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

    biascorr->Correct(m_InputImages[ichan], m_CorrectedImages[ichan], fullRes);
  }

}

void
NeonateSegmentationFilter
::NormalizePosteriors(bool fullRes)
{
  itkDebugMacro(<< "NormalizePosteriors");

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

  if (skips[0] == 0 || fullRes)
    skips[0] = 1;
  if (skips[1] == 0 || fullRes)
    skips[1] = 1;
  if (skips[2] == 0 || fullRes)
    skips[2] = 1;

  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;

        for (unsigned int iclass = 0; iclass < numClasses; iclass++)
          if (m_Posteriors[iclass]->GetPixel(ind) < 0.0)
            m_Posteriors[iclass]->SetPixel(ind, 0.0);

        double tmp = 0.0;
        for (unsigned int iclass = 0; iclass < numClasses; iclass++)
          tmp += m_Posteriors[iclass]->GetPixel(ind);

        if (tmp < 1e-20)
          continue;

        // TODO: normalize to min-max range of prob pixel type?
        // Always use float for now
        for (unsigned int iclass = 0; iclass < numClasses; iclass++)
          m_Posteriors[iclass]->SetPixel(ind, (ProbabilityImagePixelType)
            (m_Posteriors[iclass]->GetPixel(ind) / tmp));
      }
}

void
NeonateSegmentationFilter
::EMLoop()
{

  unsigned numChannels = m_InputImages.GetSize();
  unsigned numPriors = m_Priors.GetSize();

  unsigned numClasses = 0;
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

  // Initialize posterior probabilities with atlas priors
  for (unsigned iclass = 0; iclass < numClasses; iclass++)
  {
    unsigned int iprior = m_PriorLookupTable[iclass];

    ProbabilityImagePointer post = m_Posteriors[iclass];
    ProbabilityImagePointer prior = m_Priors[iprior];

    double scale = 1.0;
    if (iclass > 0)
    {
      scale =
        m_PriorWeights[iprior] * m_NumberOfGaussians[iprior];
    }
    else
    {
      scale =
        m_MyelinPriorWeight * m_NumberOfGaussians[iprior];
    }

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


  // Hack, for allocating parameters
  this->ComputeGaussians();

  for (unsigned int iprior = 1; iprior < numPriors; iprior++)
    if (m_NumberOfGaussians[iprior] > 1)
    {
      if (m_DoMSTSplit)
        this->SplitPriorMST(iprior);
      else
        this->SplitPrior(iprior);
    }

  // Compute initial distribution parameters from the atlas using robust
  // clustering (MST) and robust Gaussian estimation (FastMCD)
  this->ComputeInitialDistributions();

  muLogMacro(<< "\nInitial Gaussian parameters\n");
  muLogMacro(<< "---------------------------\n");
  muLogMacro(<< "Means =\n" << m_Means << "\n");
  for (unsigned int iclass = 0; iclass < numClasses; iclass++)
  {
    muLogMacro(
      << "Covar " << iclass+1 << ":\n" << m_Covariances[iclass] << "\n");
  }
  muLogMacro(<< "\n");

  // Use the first white matter mean as reference mean
  VectorType refMean(numChannels);
  for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    refMean[ichan] = m_Means(ichan, 3);

  // Compute initial posteriors with updated distributions
  this->ComputePosteriors(false);

  double logLikelihood = vnl_huge_val(1.0);
  double deltaLogLikelihood = 1.0;

  unsigned int biasdegree = 0;

  // EM loop
  bool converged = false;
  unsigned int iter = 0;
  while (!converged)
  {

    iter++;

    // Fix the mean for one class to avoid drift and to make
    // sure bias correction gives intensities in same range
    itkDebugMacro(<< "Fixing means of the first class");
    MatrixType offsetMeans(numChannels, numClasses, 0.0);
    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    {
      for (unsigned int iclass = 0; iclass < (numClasses-1); iclass++)
        offsetMeans(ichan, iclass) =
          -1.0 * (m_Means(ichan, 3) - refMean(ichan));
    }
    m_Means += offsetMeans;

    muLogMacro(<< "\n");
    muLogMacro(<< "EM iteration " << iter << "\n");
    muLogMacro(<< "---------------------" << "\n");
    for (unsigned int iclass = 0; iclass < numClasses; iclass++)
    {
      muLogMacro(<< "Class " << (iclass+1) << " means: ");
      for (unsigned int ichan = 0; ichan < numChannels; ichan++)
         muLogMacro(<< m_Means(ichan, iclass) << "\t");
      muLogMacro(<< "\n");
    }

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

    // Update distribution info after bias correction
    this->ComputeGaussians();

    // Recompute posteriors, not at full resolution
    this->ComputePosteriors(false);

    double prevLogLikelihood = logLikelihood;
    if (prevLogLikelihood == 0)
      prevLogLikelihood = DBL_EPSILON;

    // Compute log-likelihood and normalize posteriors
    muLogMacro(<< "Computing log likelihood...\n");
    logLikelihood = 0;
    for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
      for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
        for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
        {
          if (m_Mask->GetPixel(ind) == 0)
            continue;

          double tmp = DBL_EPSILON;
          for (unsigned int iclass = 0; iclass < numClasses; iclass++)
            tmp += m_Posteriors[iclass]->GetPixel(ind);
          logLikelihood += log(tmp);
        }
    muLogMacro(<< "log(likelihood) = " << logLikelihood << "\n");

    this->NormalizePosteriors(false);

    deltaLogLikelihood =
      fabs((logLikelihood - prevLogLikelihood) / prevLogLikelihood);
      //(logLikelihood - prevLogLikelihood) / fabs(prevLogLikelihood);

    muLogMacro(<< "delta log(likelihood) = " << deltaLogLikelihood << "\n");

    // Convergence check
    converged =
      (iter >= m_MaximumIterations)
      ||
      (deltaLogLikelihood < 0)
      ||
      ((deltaLogLikelihood < m_LikelihoodTolerance)
        &&
        (biasdegree == m_MaxBiasDegree));

  } // end EM loop

  muLogMacro(<< "Done computing posteriors with " << iter << " iterations\n");

  // Bias correction at full resolution, still using downsampled data
  // for computing the coeficients though
  if (m_MaxBiasDegree > 0)
  {
    this->CorrectBias(biasdegree, true);
    this->ComputeGaussians();
  }

  // Compute posteriors at full resolution
  this->ComputePosteriors(true);

  // Normalize posteriors at full resolution
  this->NormalizePosteriors(true);

  // Compute labels from EM posteriors
  this->ComputeLabelImage();

}

void
NeonateSegmentationFilter
::ComputeLabelImage()
{

  itkDebugMacro(<< "ComputeLabelImage");

  unsigned int numPriors = m_Priors.GetSize();

  unsigned numClasses = m_Posteriors.GetSize();

  unsigned int numFGClasses =
    numClasses - m_NumberOfGaussians[numPriors-1];

  InputImageRegionType region = m_InputImages[0]->GetLargestPossibleRegion();

  m_LabelImage = ByteImageType::New();
  m_LabelImage->SetRegions(region);
  m_LabelImage->Allocate();
  m_LabelImage->SetOrigin(m_InputImages[0]->GetOrigin());
  m_LabelImage->SetSpacing(m_InputImages[0]->GetSpacing());

  ByteImagePointer mask = ByteImageType::New();
  mask->SetRegions(region);
  mask->Allocate();
  mask->SetOrigin(m_InputImages[0]->GetOrigin());
  mask->SetSpacing(m_InputImages[0]->GetSpacing());

  ByteImageIndexType ind;
  ByteImageSizeType size = m_LabelImage->GetLargestPossibleRegion().GetSize();

  // Mersenne twister
  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {

        if (m_Mask->GetPixel(ind) == 0)
        {
          m_LabelImage->SetPixel(ind, 0);
          mask->SetPixel(ind, 0);
          continue;
        }

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
          if (fabs(v - maxv) < 1e-20)
            tieIndices.Append(iclass);
        }

        if (tieIndices.GetSize() > 1)
        {
          unsigned int whichIndex =
            (unsigned int)
              rng->GenerateUniformIntegerUpToK(tieIndices.GetSize() - 1);
          imax = tieIndices[whichIndex];
          maxv = m_Posteriors[imax]->GetPixel(ind);
        }

        unsigned char label = 0;
        unsigned char isFG = 0;

        // Only use non-zero probabities and foreground classes
        if (maxv > 0 && imax < numFGClasses)
        {
          label = (unsigned char)(imax+1);
          isFG = 1;
        }

        m_LabelImage->SetPixel(ind, label);
        mask->SetPixel(ind, isFG);

      } // for ind


  // Binary opening
  typedef itk::BinaryBallStructuringElement<unsigned char, 3> StructElementType;
  typedef
    itk::BinaryDilateImageFilter<ByteImageType, ByteImageType,
      StructElementType> DilateType;
  typedef
    itk::BinaryErodeImageFilter<ByteImageType, ByteImageType,
      StructElementType> ErodeType;

  StructElementType structel;

  StructElementType::RadiusType radius;
  for (unsigned int i = 0; i < 3; i++)
    radius[i] = 2;
  structel.SetRadius(radius);

  ErodeType::Pointer erode = ErodeType::New();
  erode->SetErodeValue(1);
  erode->SetKernel(structel);
  erode->SetInput(mask);

  erode->Update();

  DilateType::Pointer dil = DilateType::New();
  dil->SetDilateValue(1);
  dil->SetKernel(structel);
  dil->SetInput(erode->GetOutput());

  dil->Update();

  // Find largest cluster
  typedef itk::Image<unsigned short, 3> UShortImageType;
  typedef UShortImageType::Pointer UShortImagePointer;

  typedef ConnectedComponentsFilter<ByteImageType, UShortImageType> CCType;

  CCType::Pointer ccfilter = CCType::New();
  ccfilter->SetInput(dil->GetOutput());
  ccfilter->Update();

  UShortImagePointer components = ccfilter->GetOutput();

  // Strip floating clusters
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        // Update brain mask
        if (m_LabelImage->GetPixel(ind) == 0)
          m_Mask->SetPixel(ind, 0);
        else
          m_Mask->SetPixel(ind, 1);

        // Strip floating voxels
        if (components->GetPixel(ind) != 1)
        {
          m_LabelImage->SetPixel(ind, 0);
          m_Mask->SetPixel(ind, 0);
        }

        // Relabel non-brain
        if (m_LabelImage->GetPixel(ind) > numFGClasses)
        {
          m_LabelImage->SetPixel(ind, 0);
          m_Mask->SetPixel(ind, 0);
        }

      }

  // Set brain class posteriors to zero if mask is 0
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        if (m_Mask->GetPixel(ind) == 0)
        {
          for (unsigned int i = 0; i < numFGClasses; i++)
            m_Posteriors[i]->SetPixel(ind, 0);
          double sump = DBL_EPSILON;
          for (unsigned int i = numFGClasses; i < numClasses; i++)
            sump += m_Posteriors[i]->GetPixel(ind);
          for (unsigned int i = numFGClasses; i < numClasses; i++)
            m_Posteriors[i]->SetPixel(ind, m_Posteriors[i]->GetPixel(ind)/sump);
        }
      }


#if 0
  // Clean up white matter probability
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        //if (m_Posteriors[0]->GetPixel(ind) > 0.05)
        if (m_LabelImage->GetPixel(ind) == 1)
          mask->SetPixel(ind, 1);
        else
          mask->SetPixel(ind, 0);
      }

  CCType::Pointer wmccfilter = CCType::New();
  wmccfilter->SetInput(mask);
  wmccfilter->Update();

  components = wmccfilter->GetOutput();

  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        if (components->GetPixel(ind) != 1)
          m_Posteriors[0]->SetPixel(ind, 0);
      }

  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        //if (m_Posteriors[1]->GetPixel(ind) > 0.05)
        if (m_LabelImage->GetPixel(ind) == 2)
          mask->SetPixel(ind, 1);
        else
          mask->SetPixel(ind, 0);
      }

  wmccfilter->SetInput(mask);
  wmccfilter->Update();

  components = wmccfilter->GetOutput();

  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        if (components->GetPixel(ind) != 1)
          m_Posteriors[1]->SetPixel(ind, 0);
      }

  // Normalize class posteriors after wm cleaning
  // Recompute labels
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
      {
        if (m_Mask->GetPixel(ind) != 0)
        {
          unsigned int imax = 0;
          double maxp = 0;
          double sump = DBL_EPSILON;
          for (unsigned int i = 0; i < numFGClasses; i++)
          {
            double p = m_Posteriors[i]->GetPixel(ind);
            sump += p;
            if (p > maxp)
            {
              maxp = p;
              imax = i;
            }
          }

          for (unsigned int i = 0; i < numFGClasses; i++)
            m_Posteriors[i]->SetPixel(ind, m_Posteriors[i]->GetPixel(ind)/sump);

          if (maxp > 0)
            m_LabelImage->SetPixel(ind, imax+1);
        }
      }
#endif

}

void
NeonateSegmentationFilter
::ComputePrototypesMSTFull()
{

  itkDebugMacro(<< "ComputePrototypesMSTFull");

  unsigned int numChannels = m_InputImages.GetSize();

  unsigned int numPriors = m_Priors.GetSize();

  unsigned int numFGClasses = 0;
  for (unsigned int i = 0; i < (numPriors-1); i++)
    numFGClasses += m_NumberOfGaussians[i];

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

// Don't subsample?
#if 1
  skips[0] = 1;
  skips[1] = 1;
  skips[2] = 1;
#endif

  // Find range of G image
  double maxG = -vnl_huge_val(0.0);
  double minG = vnl_huge_val(0.0);
  for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;
        double g = m_GImage->GetPixel(ind);
        if (g < minG)
          minG = g;
        if (g > maxG)
          maxG = g;
      }
  double rangeG = maxG - minG;

  muLogMacro(<< "Max gradient magnitude = " << maxG << "\n");

  // Clamp threshold
  double upperG = maxG - 0.01*rangeG;

  // Compute gradient threshold
  double sumG = 0;
  double sumWM = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;
        double g = m_GImage->GetPixel(ind);
        if (g > upperG)
          g = upperG;
        double p =
          m_Posteriors[0]->GetPixel(ind) + m_Posteriors[1]->GetPixel(ind);
        sumG += p * g;
        sumWM += p;
      }

  if (sumWM < 1e-20)
    itkExceptionMacro(<< "Zero sum of white matter priors");

  double wmAvgGrad = sumG / sumWM;

  // Create list of possible samples
  DynArray<MSTClusteringProcess::VertexType> sampleList;

  // Create list of initial labels
  DynArray<unsigned char> initLabels;

  double csfReference = 0.0;

  // Get samples
  for (int iclass = (numFGClasses-1); iclass >= 0; iclass--)
  {
    muLogMacro(<< "Obtaining samples for class " << iclass+1 << "\n");

    ProbabilityImagePointer probImg = m_Priors[m_PriorLookupTable[iclass]];

    ByteImagePointer labelImg = m_LabelImage;

    double gradThres = wmAvgGrad;
    //if (iclass > 1)
    if (iclass == 3)
      gradThres = vnl_huge_val(1.0);

    MatrixType samples = this->LabelSampling(labelImg, iclass+1, gradThres, 0);

#if 0
    // Pre-prune "unimodal" samples using FastMCD
    {
      muLogMacro(<< "  Running MCD...\n");
      FastMCDSampleFilter mcdf;
      mcdf.SetMahalanobisThreshold(m_MahalanobisThreshold);
      mcdf.SetMaxCStepIterations(100);
      samples = mcdf.GetInliers(samples);
      muLogMacro(<< "  MCD gives " << samples.rows() << " inliers\n");
    }
#endif

    if (iclass == 3) 
    {
      csfReference = 0.0;
      for (unsigned int r = 0; r < samples.rows(); r++)
        csfReference += samples(r, m_ReferenceImageIndex);
      csfReference /= samples.rows();
      muLogMacro(
        << "  csf mean reference modality value = " << csfReference << "\n");
    }

    unsigned char mark = iclass+1;

    unsigned int count = 0;
    for (unsigned int r = 0; r < samples.rows(); r++)
    {
#if DROP_EXTREME
      // T1w: Pre-prune wm/gm samples darker than csf mean
      if ((iclass < 3) && (m_ReferenceModality == T1)
          &&
          (samples(r, m_ReferenceImageIndex) < csfReference))
        continue;
      // T2w: Pre-prune wm/gm samples brighter than csf mean
      if ((iclass < 3) && (m_ReferenceModality == T2)
          &&
          (samples(r, m_ReferenceImageIndex) > csfReference))
        continue;
#endif

      MSTClusteringProcess::VertexType x(numChannels);
      for (unsigned int k = 0; k < numChannels; k++)
        x[k] = samples(r, k);

      sampleList.Append(x);
      initLabels.Append(mark);

      count++;
    }

    if (iclass < 3)
    {
      muLogMacro(
        << "  Class " << iclass+1 << ": " << count << " samples"
        << " after removing extreme ref image intensities\n");
    }

  } // for iclass

  unsigned int totalNumSamples = sampleList.GetSize();
  muLogMacro(<< "Total #samples = " << totalNumSamples << "\n");

  // Build MST
  muLogMacro(<< "Building MST...\n");
  MSTClusteringProcess mstProc;
  mstProc.SetInputVertices(sampleList);
  mstProc.SortOn();

  // Generate clusters from MST
  muLogMacro(<< "Clustering from MST...\n");

  DynArray<unsigned char> clusterLabels;
  clusterLabels.Allocate(totalNumSamples);
  for (unsigned int k = 0; k < totalNumSamples; k++)
    clusterLabels.Append(0);

  DynArray<unsigned int> treeMap;
  treeMap.Initialize(totalNumSamples, 0);

  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  // Break edges
  double T;
  for(T = 2.5; T >= 0.999; T -= 0.01)
  {

    muLogMacro(<< "--------------\n");
    muLogMacro(<< "Breaking edges, T = " << T << "\n");

    unsigned int numClusters = mstProc.GetClusters(treeMap.GetRawArray(), T);

    if (numClusters < numFGClasses)
    {
      muLogMacro(<< "Less than " << numFGClasses << " clusters, skipping...\n");
      continue;
    }

    // Find the label majority for each non-empty cluster
    DynArray<unsigned int> majorLabels;
    majorLabels.Initialize(numClusters, 0);

    for(unsigned int k = 0; k < numClusters; k++)
    {
      DynArray<unsigned int> labelCount;
      labelCount.Initialize(numFGClasses, 0);

      for (unsigned int i = 0; i < totalNumSamples; i++)
      {
        if (treeMap[i] == k)
          labelCount[initLabels[i]-1]++;
      }

      unsigned int imax = 0;
      unsigned int cmax = labelCount[0];
      for (unsigned int j = 1; j < numFGClasses; j++)
        if (labelCount[j] > cmax)
        {
          imax = j;
          cmax = labelCount[j];
        }

      DynArray<unsigned int> tieLabels;
      tieLabels.Allocate(numFGClasses);
      tieLabels.Append(imax);

      for (unsigned int j = 0; j < numFGClasses; j++)
        if ((j != imax) && (labelCount[j] == cmax))
        {
          tieLabels.Append(j);
        }

      if (tieLabels.GetSize() > 1)
      {
        unsigned int whichLabel =
          (unsigned int)
            rng->GenerateUniformIntegerUpToK(tieLabels.GetSize() - 1);
        imax = tieLabels[whichLabel];
      }

      majorLabels[k] = imax+1;

      //muLogMacro(<< "Found majority " << imax+1 << "\n");

    } // for k

    // Search for the IDs of the major clusters
    muLogMacro(<< "Search for major clusters...\n");

    DynArray<unsigned int> majorIDs;
    majorIDs.Initialize(numFGClasses, 0);

    DynArray<bool> foundCluster;
    foundCluster.Initialize(numFGClasses, true);

    // Find major clusters for WM1, WM2, GM, CSF
    for (unsigned int m = 0; m < numFGClasses; m++)
    {
      unsigned int label = m+1;

      unsigned int ic = 0;
      for (ic = 0; ic < numClusters; ic++)
      {
        if (majorLabels[ic] == label)
          break;
      }
      if (ic < numClusters)
      {
        majorIDs[m] = ic;
      }
      else
      {
        foundCluster[m] = false; // Search failed
      }
    }

    bool allOK = true;
    for (unsigned int m = 0; m < numFGClasses; m++)
    {
      if (!foundCluster[m])
      {
        muLogMacro(<< "X Failed search for class " << m+1 << "\n");
        allOK = false;
        //break;
      }
    }

    if (!allOK)
    {
      continue;
    }

    // Check if the sequence of samples along first feature is valid
    // Use means or medians for location estimate
    DynArray<double> locations;
    locations.Initialize(numFGClasses, 0);

    DynArray<unsigned int> majorSizes;
    majorSizes.Initialize(numFGClasses, 0);

    for (unsigned int m = 0; m < numFGClasses; m++)
    {
      unsigned int id = majorIDs[m];

      unsigned int label = m+1;

      unsigned int n = 0;
      for (unsigned int k = 0; k < totalNumSamples; k++)
      {
        if (treeMap[k] == id && initLabels[k] == label)
          n++;
      }

      majorSizes[m] = n;

      muLogMacro(<< "[" << m << "] n = " << n << " ");
      muLogMacro(<< "id = " << id << "\n");

// MCD on reference intensity only?
// or use full multi-dimensional MCD?
#if MCD_REFONLY
      MatrixType mm(n, 1);
      unsigned int r = 0;
      for (unsigned int k = 0; k < totalNumSamples; k++)
      {
        if (treeMap[k] == id && initLabels[k] == label)
        {
          SampleType x = sampleList[k];
          mm(r, 0) = x[m_ReferenceImageIndex];
          r++;
        }
      }

      // Use MCD for location estimate
      MatrixType mu(1, 1);
      MatrixType cov(numChannels, numChannels);

      FastMCDSampleFilter mcdf;
      mcdf.SetMaxCStepIterations(100);
      mcdf.GetRobustEstimate(mu, cov, mm);

      locations[m] = mu(0, 0);
#else
      MatrixType mm(n, numChannels);
      unsigned int r = 0;
      for (unsigned int k = 0; k < totalNumSamples; k++)
      {
        if (treeMap[k] == id && initLabels[k] == label)
        {
          SampleType x = sampleList[k];
          for (unsigned int c = 0; c < numChannels; c++)
            mm(r, c) = x[c];
          r++;
        }
      }

      // Use MCD for location estimate
      MatrixType mu(1, numChannels);
      MatrixType cov(numChannels, numChannels);

      FastMCDSampleFilter mcdf;
      mcdf.SetMaxCStepIterations(100);
      mcdf.GetRobustEstimate(mu, cov, mm);

      locations[m] = mu(0, m_ReferenceImageIndex);
#endif


    }

    muLogMacro(
      << "WM locations = " << locations[0] << ", " << locations[1] << "\n");
    muLogMacro(
     << "WM major IDs = " << majorIDs[0] << ", " << majorIDs[1] << "\n");

    muLogMacro(<< "Locations = \n  ");
    for (unsigned int i = 0; i < numFGClasses; i++)
      muLogMacro(<< locations[i] << " ");
    muLogMacro(<< "\n");

    bool goodseq = false;
    if (m_ReferenceModality == T1)
      goodseq = (locations[0] > locations[2]) && (locations[2] > locations[1])
        && (locations[1] > locations[3]);
    if (m_ReferenceModality == T2)
      goodseq = (locations[0] < locations[2]) && (locations[2] < locations[1])
        && (locations[1] < locations[3]);

    // Check size
    // TODO: Check proportion instead?
    for (unsigned int m = 0; m < numFGClasses; m++)
      if (majorSizes[m] < MIN_CLASS_SAMPLES)
        goodseq = false;

    if (goodseq)
    {
      // Update final label list, drop false positives
      muLogMacro(<< "Clusters found, updating final labels...\n");
      for (unsigned int m = 0; m < numFGClasses; m++)
      {
        unsigned int id = majorIDs[m];
        unsigned int label = m+1;

        for (unsigned int k = 0; k < totalNumSamples; k++)
          if (treeMap[k] == id)
          {
            if (initLabels[k] == label)
              clusterLabels[k] = m+1;
            else
              clusterLabels[k] = 0;
          }
      }
      break;
    }
    else
    {
      muLogMacro(<< "X Bad sequence\n");
    }


  } // for T

  if (T < 1.0)
    itkExceptionMacro(<< "Whole brain MST clustering: T < 1");

  //
  // Drop samples with zero labels, condense list of samples to a matrix form
  //

  unsigned int numProto = 0;
  for (unsigned int k = 0; k < totalNumSamples; k++)
    if (clusterLabels[k] != 0)
      numProto++;

  muLogMacro(<< "Allocate prototypes...\n");
  m_Prototypes.Clear();
  m_Prototypes.Allocate(numProto);
  m_PrototypeLabels.Clear();
  m_PrototypeLabels.Allocate(numProto);

  muLogMacro(<< "Prototype selection from samples...\n");
  for (unsigned int k = 0; k < totalNumSamples; k++)
  {
    if (clusterLabels[k] != 0)
    {
      m_Prototypes.Append(sampleList[k]);
      m_PrototypeLabels.Append(clusterLabels[k]);
    }
  }

}


void
NeonateSegmentationFilter
::ComputeInitialDistributions()
{

  itkDebugMacro(<< "ComputeInitialDistributions");

// Handled by ComputeGaussians()
#if 0
  // Allocate means and covariances
  m_Means = MatrixType(numChannels, numClasses);
  m_Covariances.Clear();
  for (unsigned int i = 0; i < numClasses; i++)
  {
    MatrixType c(numChannels, numChannels, 0.0);
    m_Covariances.Append(c);
  }
#endif

  unsigned int numChannels = m_InputImages.GetSize();

  unsigned int numPriors = m_Priors.GetSize();

  unsigned int numFGClasses = 0;
  for (unsigned int i = 0; i < (numPriors-1); i++)
    numFGClasses += m_NumberOfGaussians[i];

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

// Don't subsample?
#if 1
  skips[0] = 1;
  skips[1] = 1;
  skips[2] = 1;
#endif

  // Find range of G image
  double maxG = -vnl_huge_val(0.0);
  double minG = vnl_huge_val(0.0);
  for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;
        double g = m_GImage->GetPixel(ind);
        if (g < minG)
          minG = g;
        if (g > maxG)
          maxG = g;
      }
  double rangeG = maxG - minG;

  muLogMacro(<< "Max gradient magnitude = " << maxG << "\n");

  // Clamp threshold
  double upperG = maxG - 0.01*rangeG;

  // Compute gradient threshold
  double sumG = 0;
  double sumWM = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
      {
        if (m_Mask->GetPixel(ind) == 0)
          continue;
        double g = m_GImage->GetPixel(ind);
        if (g > upperG)
          g = upperG;

        double p = m_Priors[0]->GetPixel(ind);
        sumG += p * g;
        sumWM += p;
      }

  if (sumWM < 1e-20)
    itkExceptionMacro(<< "Zero sum of white matter priors");

  double wmAvgGrad = sumG / sumWM;

  double csfRef = 0.0;
  double gmRef = 0.0;

  // Get samples
  for (int iclass = (numFGClasses-1); iclass > 1; iclass--)
  {
    muLogMacro(<< "Obtaining samples for class " << iclass+1 << "\n");

    ProbabilityImagePointer probImg = m_Priors[m_PriorLookupTable[iclass]];

    ByteImagePointer labelImg = m_LabelImage;

    // Compute threshold
    double maxProb = 0.0;
    for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
      for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
        for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
        {
          if (m_Mask->GetPixel(ind) == 0)
            continue;
          if (probImg->GetPixel(ind) > maxProb)
            maxProb = probImg->GetPixel(ind);
        }

    double probThres = m_PriorThresholdFraction * maxProb;
    //double probThres = m_PriorThresholdFraction;

    double gradThres = vnl_huge_val(1.0);

    MatrixType samples =
      this->ProbabilitySampling(probImg, probThres, gradThres, 0);

    MatrixType mu(1, numChannels);
    MatrixType cov(numChannels, numChannels);

    muLogMacro(<< "  Running MCD...\n");
    FastMCDSampleFilter mcdf;
    mcdf.SetMaxCStepIterations(100);
    mcdf.GetRobustEstimate(mu, cov, samples);

    if (iclass == 3)
      csfRef = mu(0, m_ReferenceImageIndex);
    if (iclass == 2)
      gmRef = mu(0, m_ReferenceImageIndex);

    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
      m_Means(ichan, iclass) = mu(0, ichan);

    for (unsigned int i = 0; i < numChannels; i++)
      cov(i, i) = cov(i, i) + 1e-10;

    // Stretch covariance a bit, softer EM initialization
    //cov *= 1.1;

    m_Covariances[iclass] = cov;

  } // for iclass

  if ((m_ReferenceModality == T1) && (gmRef < csfRef))
    itkExceptionMacro(
      << "Invalid T1 gm and csf ordering: " << gmRef << ", " << csfRef
      << " (bad prior threshold?)");
  if ((m_ReferenceModality == T2) && (gmRef > csfRef))
    itkExceptionMacro(
      << "Invalid T2 gm and csf ordering: " << gmRef << ", " << csfRef
      << " (bad prior threshold?)");

  // Sample wm
  DynArray<MSTClusteringProcess::VertexType> wmSamples;
  {
    muLogMacro(<< "Obtaining samples for class " << 1 << "\n");

    ProbabilityImagePointer probImg = m_Priors[m_PriorLookupTable[0]];

    // Compute probability threshold
    double maxProb = 0.0;
    for (ind[2] = 0; ind[2] < size[2]; ind[2]+=skips[2])
      for (ind[1] = 0; ind[1] < size[1]; ind[1]+=skips[1])
        for (ind[0] = 0; ind[0] < size[0]; ind[0]+=skips[0])
        {
          if (m_Mask->GetPixel(ind) == 0)
            continue;
          if (probImg->GetPixel(ind) > maxProb)
            maxProb = probImg->GetPixel(ind);
        }

    double probThres = m_PriorThresholdFraction * maxProb;
    //double probThres = m_PriorThresholdFraction;

    double gradThres = wmAvgGrad;

    MatrixType samples =
      this->ProbabilitySampling(probImg, probThres, gradThres, 0);

    unsigned int count = 0;
    for (unsigned int r = 0; r < samples.rows(); r++)
    {
#if DROP_EXTREME
      // Pre-prune
      if ((m_ReferenceModality == T1)
          &&
          (samples(r, m_ReferenceImageIndex) < csfRef))
        continue;
      if ((m_ReferenceModality == T2)
          &&
          (samples(r, m_ReferenceImageIndex) > csfRef))
        continue;
#endif

      MSTClusteringProcess::VertexType x(numChannels);
      for (unsigned int k = 0; k < numChannels; k++)
        x[k] = samples(r, k);

      wmSamples.Append(x);

      count++;
    }

    muLogMacro(<< "  Class " << 1 << ": " << count << " samples\n");
  }

  unsigned int numWMSamples = wmSamples.GetSize();

  DynArray<unsigned char> wmMarks;
  wmMarks.Initialize(numWMSamples, 0);

  // Build MST
  muLogMacro(<< "Building MST...\n");
  MSTClusteringProcess mstProc;
  mstProc.SetInputVertices(wmSamples);
  mstProc.SortOn();

  // Generate clusters from MST
  muLogMacro(<< "Clustering from MST...\n");

  DynArray<unsigned int> treeMap;
  treeMap.Initialize(numWMSamples, 0);

  unsigned int lastNumClusters = 0;

  // Break edges
  double T;
  for(T = 2.5; T >= 0.999; T -= 0.01)
  {

    muLogMacro(<< "--------------\n");
    muLogMacro(<< "Breaking edges, T = " << T << "\n");

    unsigned int numClusters = mstProc.GetClusters(treeMap.GetRawArray(), T);

    if (numClusters < 2)
    {
      muLogMacro(<< "Less than " << 2 << " clusters, skipping...\n");
      continue;
    }

    if (numClusters == lastNumClusters)
    {
      muLogMacro(<< "No change from last T, skipping...\n");
      continue;
    }

    lastNumClusters = numClusters;

    muLogMacro(<< "Edges broken, #clusters = " << numClusters << "\n");

    DynArray<double> locArray;
    locArray.Initialize(numClusters, 0);
    for (unsigned int id = 0; id < numClusters; id++)
    {
      unsigned int n = 0;
      for (unsigned int k = 0; k < numWMSamples; k++)
        if (treeMap[k] == id)
          n++;

// MCD location on reference modality only?
// or on whole multi-dimensional input?
#if MCD_REFONLY
      MatrixType mm(n, 1);
      unsigned int r = 0;
      for (unsigned int k = 0; k < numWMSamples; k++)
      {
        if (treeMap[k] == id)
        {
          SampleType x = wmSamples[k];
          mm(r, 0) = x[m_ReferenceImageIndex];
          r++;
        }
      }

      // Use MCD for location estimate
      MatrixType mu(1, 1);
      MatrixType cov(numChannels, numChannels);

      FastMCDSampleFilter mcdf;
      mcdf.SetMaxCStepIterations(100);
      mcdf.GetRobustEstimate(mu, cov, mm);

      locArray[id] = mu(0, 0);
#else 
      MatrixType mm(n, numChannels);
      unsigned int r = 0;
      for (unsigned int k = 0; k < numWMSamples; k++)
      {
        if (treeMap[k] == id)
        {
          SampleType x = wmSamples[k];
          for (unsigned int c = 0; c < numChannels; c++)
            mm(r, c) = x[c];
          r++;
        }
      }

      MatrixType mu(1, numChannels);
      MatrixType cov(numChannels, numChannels);

      FastMCDSampleFilter mcdf;
      mcdf.SetMaxCStepIterations(100);
      mcdf.GetRobustEstimate(mu, cov, mm);

      locArray[id] = mu(0, m_ReferenceImageIndex);
#endif
    }

    unsigned int clusterIDs[2];

    // Find first WM cluster (myelinated)
    // T1: mu > gm
    // T2: mu < gm
    clusterIDs[0] = numClusters;
    for (unsigned int id = 0; id < numClusters; id++)
    {
      bool validOrder = false;
      if (m_ReferenceModality == T1)
        validOrder = (locArray[id] > gmRef);
      if (m_ReferenceModality == T2)
        validOrder = (locArray[id] < gmRef);

      if (validOrder)
      {
        clusterIDs[0] = id;
        break;
      }
    }

    // Find second wm cluster (non-myelinated)
    // T1: gm > mu > csf
    // T2: gm < mu < csf
    clusterIDs[1] = numClusters;
    for (unsigned int id = 0; id < numClusters; id++)
    {
      bool validOrder = false;
      if (m_ReferenceModality == T1)
        validOrder = (gmRef > locArray[id]) && (locArray[id] > csfRef);
      if (m_ReferenceModality == T2)
        validOrder = (gmRef < locArray[id]) && (locArray[id] < csfRef);

      if (validOrder)
      {
        clusterIDs[1] = id;
        break;
      }
    }

    if (clusterIDs[0] == numClusters || clusterIDs[1] == numClusters)
    {
      muLogMacro(<< "X Failed search for well-ordered WM clusters\n");
      continue;
    }

    muLogMacro(
      << "WM1, WM2, GM, CSF reference intensity locations: "
      << locArray[clusterIDs[0]]
      << ", " << locArray[clusterIDs[1]]
      << ", " << gmRef
      << ", " << csfRef << "\n");

    locArray.Clear();

    // Check size
    // TODO: Check proportion instead?
    muLogMacro(<< "WM cluster sizes: ");
    bool sizeOK = true;
    for (unsigned int m = 0; m < 2; m++)
    {
      unsigned int id = clusterIDs[m];

      unsigned int csize = 0;
      for (unsigned int k = 0; k < numWMSamples; k++)
        if (treeMap[k] == id)
          csize++;
      muLogMacro(<< csize << " ");

      if (csize < MIN_CLASS_SAMPLES)
        sizeOK = false;
    }
    muLogMacro(<< "\n");

    if (!sizeOK)
    {
      muLogMacro(<< "X Cluster sizes too small\n");
      continue;
    }

    // Update final label list, drop false positives
    if (sizeOK)
    {
      muLogMacro(<< "Clusters found, updating final labels...\n");
      for (unsigned int m = 0; m < 2; m++)
      {
        unsigned int id = clusterIDs[m];

        for (unsigned int k = 0; k < numWMSamples; k++)
          if (treeMap[k] == id)
          {
            wmMarks[k] = m+1;
          }
      }
      break;
    }

  } // for T

  if (T < 1.0)
    itkExceptionMacro(<< "WM MST clustering: T < 1");

  // Compute white matter Gaussian parameters
  for (unsigned int iclass = 0; iclass < 2; iclass++)
  {
    unsigned char mark = iclass+1;

    unsigned int count = 0;
    for (unsigned int i = 0; i < numWMSamples; i++)
      if (wmMarks[i] == mark)
        count++;

    muLogMacro(<< "WM " << iclass+1 << " #samples in cluster = " << count << "\n");
    if (count == 0)
      itkExceptionMacro(<< "Zero samples for WM " << iclass+1);

#if 1
    // Use MCD
    MatrixType mm(count, numChannels);
    unsigned int r = 0;
    for (unsigned int k = 0; k < numWMSamples; k++)
    {
      if (wmMarks[k] == mark)
      {
        SampleType x = wmSamples[k];
        for (unsigned int c = 0; c < numChannels; c++)
          mm(r, c) = x[c];
        r++;
      }
    }

    MatrixType mu(1, numChannels);
    MatrixType cov(numChannels, numChannels);

    FastMCDSampleFilter mcdf;
    mcdf.SetMaxCStepIterations(100);
    mcdf.GetRobustEstimate(mu, cov, mm);

    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
      m_Means(ichan, iclass) = mu(0, ichan);

    //cov *= 1.1;
    m_Covariances[iclass] = cov;
#else
    // Use ordinary mean - cov, already pruned / robustified
    for (unsigned int ichan = 0; ichan < numChannels; ichan++)
    {
      double mu = 0.0;
      for (unsigned int i = 0; i < numWMSamples; i++)
      {
        if (wmMarks[i] == mark)
        {
          SampleType x = wmSamples[i];
          mu += x[ichan];
        }
      }
      mu /= count;
      m_Means(ichan, iclass) = mu;
    }

    MatrixType covtmp(numChannels, numChannels, 0.0);

    if (count == 1)
    {
      m_Covariances[iclass] = covtmp;
      continue;
    }

    for (unsigned int r = 0; r < numChannels; r++)
    {
      double mu1 = m_Means(r, iclass);
      for (unsigned int c = r; c < numChannels; c++)
      {
        double mu2 = m_Means(c, iclass);

        double v = 0.0;
        double diff1 = 0;
        double diff2 = 0;
        for (unsigned int i = 0; i < numWMSamples; i++)
        {
          if (wmMarks[i] == mark)
          {
            SampleType x = wmSamples[i];
            diff1 = x[r] - mu1;
            diff2 = x[c] - mu2;
            v += (diff1*diff2);
          }
        }
        v /= (count-1);

        if (r == c)
          v += 1e-10;

        covtmp(r, c) = v;
        covtmp(c, r) = v;
      }
    }

    //covtmp *= 1.1;

    m_Covariances[iclass] = covtmp;
#endif

  }

}

NeonateSegmentationFilter::MatrixType
NeonateSegmentationFilter
::ProbabilitySampling(
  ProbabilityImagePointer probImg, double tau, double gamma,
  double sampleSpacing)
{

#if 0
  muLogMacro(
    << "  Obtaining samples with probability > " << tau << ", gradient < "
    << gamma << "\n");
#else
  muLogMacro(
    << "  Obtaining samples with gradient < " << gamma << "\n");
#endif

  unsigned numChannels = m_InputImages.GetSize();

  ProbabilityImageIndexType ind;

  ProbabilityImageSizeType size = probImg->GetLargestPossibleRegion().GetSize();

  InputImageSpacingType spacing = m_CorrectedImages[0]->GetSpacing();

  InputImageOffsetType skips;
  skips[0] = (unsigned)floor(sampleSpacing / spacing[0]);
  skips[1] = (unsigned)floor(sampleSpacing / spacing[1]);
  skips[2] = (unsigned)floor(sampleSpacing / spacing[2]);

  if (skips[0] == 0)
    skips[0] = 1;
  if (skips[1] == 0)
    skips[1] = 1;
  if (skips[2] == 0)
    skips[2] = 1;

  // Mersenne twister
  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  // Allocate mask
  ByteImagePointer tmpmask = ByteImageType::New();
  tmpmask->SetRegions(probImg->GetLargestPossibleRegion());
  tmpmask->Allocate();
  tmpmask->SetOrigin(probImg->GetOrigin());
  tmpmask->SetSpacing(probImg->GetSpacing());

  // Zero the mask
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        tmpmask->SetPixel(ind, 0);

  // Threshold criteria
  unsigned int numPossibleSamples = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        double p = probImg->GetPixel(ind);
        double g = m_GImage->GetPixel(ind);
//TODO:
        //double r = rng->GenerateUniformRealClosedOpenInterval();

        //if (p > r && g < gamma)
        if (p > tau && g < gamma)
        {
          tmpmask->SetPixel(ind, 1);
          numPossibleSamples++;
        }
      }

  // Sample selection mask
  DynArray<unsigned char> selectMask;
  selectMask.Initialize(numPossibleSamples, 0);

  unsigned int numSamples = numPossibleSamples;
  if (numSamples > m_MaxSamples)
    numSamples = m_MaxSamples;

  muLogMacro(<< "  Selecting " << numSamples << " / " << numPossibleSamples << "\n");

  if (numSamples < numPossibleSamples)
  {
    unsigned int* seq =
      rng->GenerateIntegerSequence(numSamples, numPossibleSamples-1);

    for (unsigned int i = 0; i < numSamples; i++)
    {
      unsigned int which = seq[i];
      selectMask[which] = 1;
    }

    delete [] seq;
  }
  else
  {
    for (unsigned int i = 0; i < numPossibleSamples; i++)
      selectMask[i] = 1;
  }

  MatrixType samples(numSamples, numChannels);

  unsigned int pix = 0;
  unsigned int r = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        if (tmpmask->GetPixel(ind) == 0)
          continue;
        if (selectMask[pix] != 0)
        {
          for (unsigned int ichan = 0; ichan < numChannels; ichan++)
            samples(r, ichan) = m_CorrectedImages[ichan]->GetPixel(ind);
          ++r;
        }
        ++pix;
      }

  return samples;

}

NeonateSegmentationFilter::MatrixType
NeonateSegmentationFilter
::LabelSampling(
  ByteImagePointer labelImg, unsigned char label, double gamma,
  double sampleSpacing)
{

  muLogMacro(
    << "  Obtaining samples with label = " << (unsigned int)label
    << ", gradient < "     << gamma << "\n");

  unsigned numChannels = m_InputImages.GetSize();

  ProbabilityImageIndexType ind;

  ProbabilityImageSizeType size =
    labelImg->GetLargestPossibleRegion().GetSize();

  InputImageSpacingType spacing = m_CorrectedImages[0]->GetSpacing();

  InputImageOffsetType skips;
  skips[0] = (unsigned)floor(sampleSpacing / spacing[0]);
  skips[1] = (unsigned)floor(sampleSpacing / spacing[1]);
  skips[2] = (unsigned)floor(sampleSpacing / spacing[2]);

  if (skips[0] == 0)
    skips[0] = 1;
  if (skips[1] == 0)
    skips[1] = 1;
  if (skips[2] == 0)
    skips[2] = 1;

  // Allocate mask
  ByteImagePointer tmpmask = ByteImageType::New();
  tmpmask->SetRegions(labelImg->GetLargestPossibleRegion());
  tmpmask->Allocate();
  tmpmask->SetOrigin(labelImg->GetOrigin());
  tmpmask->SetSpacing(labelImg->GetSpacing());

  // Zero the mask
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        tmpmask->SetPixel(ind, 0);

  // Threshold criteria
  unsigned int numPossibleSamples = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        unsigned char c = labelImg->GetPixel(ind);
        double g = m_GImage->GetPixel(ind);

        if (c == label && g < gamma)
        {
          tmpmask->SetPixel(ind, 1);
          numPossibleSamples++;
        }
      }

  // Sample selection mask
  DynArray<unsigned char> selectMask;
  selectMask.Initialize(numPossibleSamples, 0);

  unsigned int numSamples = numPossibleSamples;
  if (numSamples > m_MaxSamples)
    numSamples = m_MaxSamples;

  muLogMacro(<< "  Selecting " << numSamples << " / " << numPossibleSamples << "\n");

  // Mersenne twister
  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  if (numSamples < numPossibleSamples)
  {
    unsigned int* seq =
      rng->GenerateIntegerSequence(numSamples, numPossibleSamples-1);

    for (unsigned int i = 0; i < numSamples; i++)
    {
      unsigned int which = seq[i];
      selectMask[which] = 1;
    }

    delete [] seq;
  }
  else
  {
    for (unsigned int i = 0; i < numPossibleSamples; i++)
      selectMask[i] = 1;
  }

  MatrixType samples(numSamples, numChannels);

  unsigned int pix = 0;
  unsigned int r = 0;
  for (ind[2] = 0; ind[2] < size[2]; ind[2] += skips[2])
    for (ind[1] = 0; ind[1] < size[1]; ind[1] += skips[1])
      for (ind[0] = 0; ind[0] < size[0]; ind[0] += skips[0])
      {
        if (tmpmask->GetPixel(ind) == 0)
          continue;
        if (selectMask[pix] != 0)
        {
          for (unsigned int ichan = 0; ichan < numChannels; ichan++)
            samples(r, ichan) = m_CorrectedImages[ichan]->GetPixel(ind);
          ++r;
        }
        ++pix;
      }

  return samples;

}
