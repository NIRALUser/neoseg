
////////////////////////////////////////////////////////////////////////////////
//
// Rescale an image (in-place) to a specified range
// Some fraction of the extreme values are clamped to output min max
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 8/2004

#ifndef _TrimmedRescaleImageModifier_h
#define _TrimmedRescaleImageModifier_h

template <class TImage>
class TrimmedRescaleImageModifier: public itk::Object
{

public:

  /** Standard class typedefs. */
  typedef TrimmedRescaleImageModifier Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** The dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TImage::ImageDimension);

  // Image types
  typedef TImage ImageType;
  typedef typename TImage::Pointer ImagePointer;
  typedef typename TImage::IndexType ImageIndexType;
  typedef typename TImage::OffsetType ImageOffsetType;
  typedef typename TImage::PixelType ImagePixelType;
  typedef typename TImage::RegionType ImageRegionType;
  typedef typename TImage::SizeType ImageSizeType;
  typedef typename TImage::SpacingType ImageSpacingType;

  void TrimAboveOn() { m_TrimAboveFlag = true; }
  void TrimAboveOff() { m_TrimAboveFlag = false; }

  void TrimBelowOn() { m_TrimBelowFlag = true; }
  void TrimBelowOff() { m_TrimBelowFlag = false; }

  itkGetConstMacro(TrimFraction, double);
  itkSetMacro(TrimFraction, double);

  itkGetConstMacro(OutputMinimum, double);
  itkSetMacro(OutputMinimum, double);

  itkGetConstMacro(OutputMaximum, double);
  itkSetMacro(OutputMaximum, double);

  void Rescale(ImagePointer img);

protected:

  TrimmedRescaleImageModifier();
  ~TrimmedRescaleImageModifier() { }

private:

  double m_TrimFraction;

  bool m_TrimAboveFlag;
  bool m_TrimBelowFlag;

  double m_OutputMinimum;
  double m_OutputMaximum;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "TrimmedRescaleImageModifier.txx"
#endif

#endif
