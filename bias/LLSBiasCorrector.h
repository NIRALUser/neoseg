
////////////////////////////////////////////////////////////////////////////////
//
// Bias field estimation using linear least squares polynomial fitting.
// Requires input image intensities and probability estimates for each voxel.
//
// Corrections are done either on subsampled grid or full resolution.
// This is designed to be run iteratively. Corrections should not be
// accumulated. Always use original image as the input, otherwise may get
// strange results.
//
// Van Leemput K, Maes F, Vandermeulen D, Suetens P. Automated model based
// bias field correction of MR images of the brain. IEEE TMI 1999; 18:885-896.
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 3/2004

#ifndef _LLSBiasCorrector_h
#define _LLSBiasCorrector_h

#include "itkImage.h"
#include "itkObject.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_qr.h"

#include "DynArray.h"

template <class TInputImage, class TProbabilityImage>
class LLSBiasCorrector : public itk::Object
{

public:

  /** Standard class typedefs. */
  typedef LLSBiasCorrector Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** The dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  // Image types
  typedef TInputImage InputImageType;
  typedef typename TInputImage::Pointer InputImagePointer;
  typedef typename TInputImage::IndexType InputImageIndexType;
  typedef typename TInputImage::PixelType InputImagePixelType;
  typedef typename TInputImage::RegionType InputImageRegionType;
  typedef typename TInputImage::SizeType InputImageSizeType;
  typedef typename TInputImage::SpacingType InputImageSpacingType;

  typedef itk::Image<unsigned char, 3> MaskImageType;
  typedef MaskImageType::Pointer MaskImagePointer;
  typedef MaskImageType::IndexType MaskImageIndexType;
  typedef MaskImageType::PixelType MaskImagePixelType;
  typedef MaskImageType::RegionType MaskImageRegionType;
  typedef MaskImageType::SizeType MaskImageSizeType;

  typedef TProbabilityImage ProbabilityImageType;
  typedef typename ProbabilityImageType::Pointer ProbabilityImagePointer;
  typedef typename ProbabilityImageType::IndexType ProbabilityImageIndexType;
  typedef typename ProbabilityImageType::PixelType ProbabilityImagePixelType;
  typedef typename ProbabilityImageType::RegionType ProbabilityImageRegionType;
  typedef typename ProbabilityImageType::SizeType ProbabilityImageSizeType;

  typedef itk::Image<float, 3> InternalImageType;
  typedef InternalImageType::Pointer InternalImagePointer;
  typedef InternalImageType::IndexType InternalImageIndexType;
  typedef InternalImageType::PixelType InternalImagePixelType;
  typedef InternalImageType::RegionType InternalImageRegionType;
  typedef InternalImageType::SizeType InternalImageSizeType;

  typedef vnl_matrix<double> MatrixType;
  typedef vnl_vector<double> VectorType;

  typedef vnl_matrix_inverse<double> MatrixInverseType;
  typedef vnl_qr<double> MatrixQRType;

// TODO: update of max degree and sample spacing should trigger mask update
// to generate new LHS matrix

  // The maximum polynomial degree of the bias field estimate
  void SetMaxDegree(unsigned int);
  itkGetMacro(MaxDegree, unsigned int);

  itkSetMacro(SampleSpacing, double);
  itkGetMacro(SampleSpacing, double);

  itkSetMacro(WorkingSpacing, double);
  itkGetMacro(WorkingSpacing, double);

  itkSetMacro(ClampBias, bool);
  itkGetMacro(ClampBias, bool);

  void AdditiveOn() { m_DoLog = false; }
  void MultiplicativeOn() { m_DoLog = true; }

  // Bias field max magnitude
  itkSetMacro(MaximumBiasMagnitude, double);
  itkGetMacro(MaximumBiasMagnitude, double);

  void SetMask(MaskImageType* mask);
  void SetProbabilities(DynArray<ProbabilityImagePointer> probs);

  void SetMeans(const VectorType& v);
  void SetVariances(const VectorType& v);

  void Correct(InputImagePointer input, InputImagePointer output, bool fullRes);

protected:

  LLSBiasCorrector();
  ~LLSBiasCorrector();

  void CheckInput();
  void ComputeDistributions();

private:

  InputImagePointer m_InputData;

  DynArray<ProbabilityImagePointer> m_Probabilities;

  bool m_DoLog;

  unsigned int m_MaxDegree;

  double m_SampleSpacing;
  double m_WorkingSpacing;

  bool m_ClampBias;
  double m_MaximumBiasMagnitude;

  VectorType m_Means;
  VectorType m_Variances;

  MaskImagePointer m_Mask;

  MatrixType m_LHS;
  MatrixType m_LHS_Qt; // Q' given LHS = Q*R

  // Coordinate scaling and offset, computed from input probabilities
  double m_XMu[3];
  double m_XStd[3];

};

#ifndef MU_MANUAL_INSTANTIATION
#include "LLSBiasCorrector.txx"
#endif

#endif
