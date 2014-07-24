
////////////////////////////////////////////////////////////////////////////////
//
// (Normalized) Mutual information metric for 3D images
//
// Based on MIRIT, uses downsampling and partial volume interpolation
//
// Maes, F., Collignon, A., Vandermeulen, D., Marchal, G., Suetens, P., April
// 1997. Multimodality image registration by maximization of mutual information.
// IEEE TMI 16 (2), 187-198.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef _NegativeMIImageMatchMetric_h
#define _NegativeMIImageMatchMetric_h

#include "itkImage.h"
#include "itkImageToImageMetric.h"
#include "itkIndex.h"
#include "itkPoint.h"
#include "itkSingleValuedCostFunction.h"

#include "vnl/vnl_matrix.h"

template <class TFixedImage, class TMovingImage>
class NegativeMIImageMatchMetric:
  public itk::ImageToImageMetric<TFixedImage, TMovingImage>
{
public:

  /** Standard class typedefs. */
  typedef NegativeMIImageMatchMetric  Self;
  typedef itk::ImageToImageMetric<TFixedImage, TMovingImage> Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NegativeMIImageMatchMetric, itk::ImageToImageMetric);

  /** Types inherited from Superclass. */
  typedef typename Superclass::TransformType            TransformType;
  typedef typename Superclass::TransformPointer         TransformPointer;
  typedef typename Superclass::TransformJacobianType    TransformJacobianType;
  typedef typename Superclass::InterpolatorType         InterpolatorType;
  typedef typename Superclass::MeasureType              MeasureType;
  typedef typename Superclass::DerivativeType           DerivativeType;
  typedef typename Superclass::ParametersType           ParametersType;
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::FixedImageConstPointer   FixedImageConstPointer;
  typedef typename Superclass::MovingImageConstPointer  MovingImageCosntPointer;

  // Image parameters type
  typedef typename FixedImageType::IndexType FixedImageIndexType;
  typedef typename FixedImageType::SizeType FixedImageSizeType;
  typedef typename FixedImageType::SpacingType FixedImageSpacingType;

  typedef typename MovingImageType::IndexType MovingImageIndexType;
  typedef typename MovingImageType::SizeType MovingImageSizeType;
  typedef typename MovingImageType::SpacingType MovingImageSpacingType;

  typedef typename TransformType::InputPointType FixedImagePointType;
  typedef typename TransformType::OutputPointType MovingImagePointType;

  // Image type containing histogram indices
  typedef itk::Image<unsigned int, 3> IndexImageType;
  typedef IndexImageType::IndexType IndexImageIndexType;
  typedef IndexImageType::Pointer IndexImagePointer;
  typedef IndexImageType::RegionType IndexImageRegionType;
  typedef IndexImageType::SizeType IndexImageSizeType;
  typedef IndexImageType::SpacingType IndexImageSpacingType;

  typedef vnl_matrix<double> HistogramType;

  /** Enum of the moving image dimension. */
  itkStaticConstMacro(MovingImageDimension, unsigned int,
                      MovingImageType::ImageDimension);

  /** Get the derivatives of the match measure. */
  void GetDerivative(
    const ParametersType& parameters,
    DerivativeType & Derivative ) const;
  void GetStochasticDerivative(
    const ParametersType& parameters,
    DerivativeType & Derivative ) const;

  /**  Get the value. */
  MeasureType GetValue( const ParametersType& parameters ) const;

  /**  Get the value and derivatives for single valued optimizers. */
  void GetValueAndDerivative( const ParametersType& parameters,
    MeasureType& Value, DerivativeType& Derivative ) const;

  void SetFixedImage(const FixedImageType* img);
  void SetMovingImage(const MovingImageType* img);

  itkGetConstMacro(SampleSpacing, double);
  void SetSampleSpacing(double s);

  itkGetConstMacro(Normalized, bool);
  itkSetMacro(Normalized, bool);

  itkGetConstMacro(UseKMeans, bool);
  itkSetMacro(UseKMeans, bool);

  itkGetConstMacro(DerivativeStepLengths, ParametersType);
  itkSetMacro(DerivativeStepLengths, ParametersType);

  itkGetConstMacro(NumberOfBins, unsigned int);
  void SetNumberOfBins(unsigned int numbins);

  virtual unsigned int GetNumberOfParameters() const
  {
    return this->m_Transform->GetNumberOfParameters();
  }

protected:
  NegativeMIImageMatchMetric();
  virtual ~NegativeMIImageMatchMetric() { delete m_HistogramPointer; }
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

  void MapFixedImage();
  void MapMovingImage();

  void ComputeHistogram() const;
  double ComputeMI() const;

private:
  NegativeMIImageMatchMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int m_NumberOfBins;

  IndexImagePointer m_FixedIndexImage;
  IndexImagePointer m_MovingIndexImage;

  bool m_UseKMeans;

  double m_KMeansSampleSpacing;
  double m_SampleSpacing;

  unsigned int m_Skips[3];

  HistogramType* m_HistogramPointer;

  bool m_Normalized;

  ParametersType m_DerivativeStepLengths;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "NegativeMIImageMatchMetric.txx"
#endif

#endif

