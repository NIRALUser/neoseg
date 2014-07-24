
////////////////////////////////////////////////////////////////////////////////
//
// Atlas based segmentation using the Expectation Maximization algorithm
//
// Designed for 3D MRI
//
// Van Leemput K, Maes F, Vandermeulen D, Suetens P. Automated model based
// tissue classification of MR images of the brain. IEEE TMI 1999; 18:897-908.
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 3/2004

#ifndef _EMSegmentationFilter_h
#define _EMSegmentationFilter_h

#include "itkArray.h"
#include "itkImage.h"
#include "itkObject.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

#include "vnl/algo/vnl_matrix_inverse.h"

#include "DynArray.h"

template <class TInputImage, class TProbabilityImage>
class EMSegmentationFilter: public itk::Object
{

public:

  // Standard class typedefs
  typedef EMSegmentationFilter Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  // Method for creation through the object factory
  itkNewMacro(Self);

  // The dimension of the image we're working with
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  // Image types
  typedef TInputImage InputImageType;
  typedef typename TInputImage::Pointer InputImagePointer;
  typedef typename TInputImage::IndexType InputImageIndexType;
  typedef typename TInputImage::OffsetType InputImageOffsetType;
  typedef typename TInputImage::PixelType InputImagePixelType;
  typedef typename TInputImage::RegionType InputImageRegionType;
  typedef typename TInputImage::SizeType InputImageSizeType;
  typedef typename TInputImage::SpacingType InputImageSpacingType;

  typedef itk::Image<unsigned char, itkGetStaticConstMacro(ImageDimension)> ByteImageType;
  typedef typename ByteImageType::Pointer ByteImagePointer;
  typedef typename ByteImageType::IndexType ByteImageIndexType;
  typedef typename ByteImageType::OffsetType ByteImageOffsetType;
  typedef typename ByteImageType::PixelType ByteImagePixelType;
  typedef typename ByteImageType::RegionType ByteImageRegionType;
  typedef typename ByteImageType::SizeType ByteImageSizeType;

  typedef itk::Image<short, itkGetStaticConstMacro(ImageDimension)> ShortImageType;
  typedef typename ShortImageType::Pointer ShortImagePointer;
  typedef typename ShortImageType::IndexType ShortImageIndexType;
  typedef typename ShortImageType::OffsetType ShortImageOffsetType;
  typedef typename ShortImageType::PixelType ShortImagePixelType;
  typedef typename ShortImageType::RegionType ShortImageRegionType;
  typedef typename ShortImageType::SizeType ShortImageSizeType;

  typedef TProbabilityImage ProbabilityImageType;
  typedef typename ProbabilityImageType::Pointer ProbabilityImagePointer;
  typedef typename ProbabilityImageType::IndexType ProbabilityImageIndexType;
  typedef typename ProbabilityImageType::OffsetType ProbabilityImageOffsetType;
  typedef typename ProbabilityImageType::PixelType ProbabilityImagePixelType;
  typedef typename ProbabilityImageType::RegionType ProbabilityImageRegionType;
  typedef typename ProbabilityImageType::SizeType ProbabilityImageSizeType;
  typedef typename ProbabilityImageType::SpacingType ProbabilityImageSpacingType;

  typedef vnl_vector<double> VectorType;
  typedef vnl_matrix<double> MatrixType;
  typedef vnl_matrix_inverse<double> MatrixInverseType;

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

  void SetInputImages(DynArray<InputImagePointer> data);

  void SetPriors(DynArray<ProbabilityImagePointer> probs);

  void SetPriorWeights(VectorType w);

  unsigned int* GetNumberOfGaussians() { return m_NumberOfGaussians; }
  void SetNumberOfGaussians(unsigned int* n);

  ByteImagePointer GetOutput();

  DynArray<ByteImagePointer> GetBytePosteriors();
  DynArray<ShortImagePointer> GetShortPosteriors();
  DynArray<ProbabilityImagePointer> GetPosteriors();

  DynArray<InputImagePointer> GetCorrected();

  void Update();

  itkGetMacro(FOVMask, ByteImagePointer);
  itkSetMacro(FOVMask, ByteImagePointer);

  itkGetMacro(DoMSTSplit, bool);
  itkSetMacro(DoMSTSplit, bool);

protected:

  EMSegmentationFilter();
  ~EMSegmentationFilter();

  void CheckInput();

  void ComputeMask();
  void ComputePriorLookupTable();

  // Split the initial parameters of classes with the same prior using
  // MST clustering
  void SplitPriorMST(unsigned int iprior);

  void SplitPrior(unsigned int iprior);

  void ComputeDistributions();
  void ComputePosteriors(bool fullRes);

  void CorrectBias(unsigned int degree, bool fullRes);

  void EMLoop();

  void ComputeLabels();

  void CleanUp();

private:

  DynArray<InputImagePointer> m_InputImages;
  DynArray<InputImagePointer> m_CorrectedImages;

  DynArray<ProbabilityImagePointer> m_Priors;
  DynArray<ProbabilityImagePointer> m_Posteriors;

  ByteImagePointer m_Mask;

  ByteImagePointer m_Labels;

  bool m_InputModified;

  double m_SampleSpacing;

  unsigned int m_MaxBiasDegree;
  double m_BiasLikelihoodTolerance;
  double m_LikelihoodTolerance;
  unsigned m_MaximumIterations;

  VectorType m_PriorWeights;

  // Number of Gaussian distributions for each associated prior
  unsigned int* m_NumberOfGaussians;

  unsigned int* m_PriorLookupTable;

  MatrixType m_Means;
  DynArray<MatrixType> m_Covariances;

  ByteImagePointer m_FOVMask;

  bool m_DoMSTSplit;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "EMSegmentationFilter.txx"
#endif

#endif
