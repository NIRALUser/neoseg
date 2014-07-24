
/*******************************************************************************

Blah...

*******************************************************************************/

#ifndef _NeonateSegmentationFilter_h
#define _NeonateSegmentationFilter_h

#include "itkArray.h"
#include "itkImage.h"
#include "itkObject.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"

#include "DynArray.h"

class NeonateSegmentationFilter: public itk::Object
{

public:

  // Standard ITK class typedefs
  typedef NeonateSegmentationFilter Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  // Method for creation through the object factory
  itkNewMacro(Self);

  // Image types
  typedef itk::Image<float, 3> InputImageType;
  typedef InputImageType::Pointer InputImagePointer;
  typedef InputImageType::IndexType InputImageIndexType;
  typedef InputImageType::OffsetType InputImageOffsetType;
  typedef InputImageType::PixelType InputImagePixelType;
  typedef InputImageType::SizeType InputImageSizeType;
  typedef InputImageType::RegionType InputImageRegionType;
  typedef InputImageType::SpacingType InputImageSpacingType;

  typedef itk::Image<unsigned char, 3> ByteImageType;
  typedef ByteImageType::Pointer ByteImagePointer;
  typedef ByteImageType::IndexType ByteImageIndexType;
  typedef ByteImageType::OffsetType ByteImageOffsetType;
  typedef ByteImageType::PixelType ByteImagePixelType;
  typedef ByteImageType::RegionType ByteImageRegionType;
  typedef ByteImageType::SizeType ByteImageSizeType;
  typedef ByteImageType::SpacingType ByteImageSpacingType;

  typedef itk::Image<float, 3> ProbabilityImageType;
  typedef ProbabilityImageType::Pointer ProbabilityImagePointer;
  typedef ProbabilityImageType::IndexType ProbabilityImageIndexType;
  typedef ProbabilityImageType::OffsetType ProbabilityImageOffsetType;
  typedef ProbabilityImageType::PixelType ProbabilityImagePixelType;
  typedef ProbabilityImageType::RegionType ProbabilityImageRegionType;
  typedef ProbabilityImageType::SizeType ProbabilityImageSizeType;
  typedef ProbabilityImageType::SpacingType ProbabilityImageSpacingType;

  typedef vnl_vector<float> VectorType;
  typedef vnl_matrix<float> MatrixType;

  typedef vnl_matrix_inverse<float> MatrixInverseType;
  typedef vnl_symmetric_eigensystem<float> MatrixEigenType;

  typedef vnl_vector<float> SampleType;

  typedef enum{T1, T2} ImageModalityEnumType;

  // Set/Get the maximum polynomial degree of the bias field estimate
  itkSetMacro(MaxBiasDegree, unsigned int);
  itkGetMacro(MaxBiasDegree, unsigned int);

  itkSetMacro(BiasLikelihoodTolerance, double);
  itkGetMacro(BiasLikelihoodTolerance, double);

  itkSetMacro(LikelihoodTolerance, double);
  itkGetMacro(LikelihoodTolerance, double);

  itkSetMacro(MaximumIterations, unsigned int);
  itkGetMacro(MaximumIterations, unsigned int);

  itkSetMacro(SampleSpacing, double);
  itkGetMacro(SampleSpacing, double);

  itkSetMacro(MahalanobisThreshold, double);
  itkGetMacro(MahalanobisThreshold, double);

  itkSetMacro(KernelWidthFraction, double);
  itkGetMacro(KernelWidthFraction, double);

  itkSetMacro(PriorThresholdFraction, double);
  itkGetMacro(PriorThresholdFraction, double);

  itkSetMacro(MyelinPriorWeight, double);
  itkGetMacro(MyelinPriorWeight, double);

  itkSetMacro(ReferenceModality, ImageModalityEnumType);
  itkGetMacro(ReferenceModality, ImageModalityEnumType);

  itkSetMacro(ReferenceImageIndex, unsigned int);
  itkGetMacro(ReferenceImageIndex, unsigned int);

  itkSetMacro(DoMSTSplit, bool);
  itkGetMacro(DoMSTSplit, bool);

  void SetInputImages(DynArray<InputImagePointer> data);

  void SetPriors(DynArray<ProbabilityImagePointer> probs);

  void SetPriorWeights(const VectorType& w);
  void SetNumberOfGaussians(const unsigned int* c);

  ByteImagePointer GetOutput();

  DynArray<ByteImagePointer> GetBytePosteriors();
  DynArray<ProbabilityImagePointer> GetPosteriors();

  //TODO: ???
  //DynArray<ShortImagePointer> GetShortPosteriors();

  DynArray<InputImagePointer> GetCorrected();

  void Update();

  // Refine the current segmentation with KNN or kernel density estimation
  void DoKNearestNeighbors();
  void DoKernelDensityEstimation();

  itkSetMacro(FOVMask, ByteImagePointer);
  itkGetMacro(FOVMask, ByteImagePointer);

protected:

  NeonateSegmentationFilter();
  ~NeonateSegmentationFilter();

  void CheckInput();

  void ComputePriorLookupTable();
  void ComputeInitialMask();

  // Compute image of 2-norm of gradient magnitudes
  void ComputeGImage();

  void SplitPriorMST(unsigned int i);
  void SplitPrior(unsigned int i);

  // Initial Gaussian params with MCD and MST
  void ComputeInitialDistributions();

  // Gaussian params from posterior probs
  void ComputeGaussians();

  void ComputePosteriors(bool fullRes);

  void CorrectBias(unsigned int degree, bool fullRes);

  void NormalizePosteriors(bool fullRes);

  // EM segmentation - bias correction loop
  void EMLoop();

  // Compute the labels from class posteriors
  void ComputeLabelImage();

  void RescaleCorrected();

  // Compute prototypes from label or priors
  //void ComputePrototypesMSTFull(bool useLabels=false);
  void ComputePrototypesMSTFull();

  MatrixType
    ProbabilitySampling(ProbabilityImagePointer p, double tau, double gamma,
      double sampleSpacing);
  MatrixType
    LabelSampling(ByteImagePointer p, unsigned char label, double gamma,
      double sampleSpacing);

private:

  DynArray<InputImagePointer> m_InputImages;
  DynArray<InputImagePointer> m_CorrectedImages;

  DynArray<ProbabilityImagePointer> m_Priors;
  DynArray<ProbabilityImagePointer> m_Posteriors;

  ByteImagePointer m_Mask;

  ByteImagePointer m_LabelImage;

  InputImagePointer m_GImage;

  bool m_InputModified;

  // Desired grid sample spacing in mm
  double m_SampleSpacing;

  unsigned int m_MaxBiasDegree;

  // EM iteration parameters
  double m_BiasLikelihoodTolerance;
  double m_LikelihoodTolerance;
  unsigned m_MaximumIterations;

  VectorType m_PriorWeights;

  // Number of Gaussian distributions for each associated prior
  unsigned int* m_NumberOfGaussians;

  // Map class index to prior index
  unsigned int* m_PriorLookupTable;

  MatrixType m_Means;
  DynArray<MatrixType> m_Covariances;

  DynArray<SampleType> m_Prototypes;
  DynArray<unsigned char> m_PrototypeLabels;

  ImageModalityEnumType m_ReferenceModality;

  // Which image is the reference
  unsigned int m_ReferenceImageIndex;

  unsigned int m_RescaledMax;

  unsigned int m_MaxSamples;

  ByteImagePointer m_FOVMask;

  bool m_DoMSTSplit;

  double m_MahalanobisThreshold;

  double m_KernelWidthFraction;

  double m_PriorThresholdFraction;

  double m_MyelinPriorWeight;

  bool m_DoneEM;
  bool m_DoneNPPrototypes;

};

#endif
