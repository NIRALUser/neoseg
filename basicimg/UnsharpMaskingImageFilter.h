
////////////////////////////////////////////////////////////////////////////////
//
// Edge enhancement by unsharp masking (blurring + addition)
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 2/2004

#ifndef _UnsharpMaskingImageFilter_h
#define _UnsharpMaskingImageFilter_h

template <class TImage>
class UnsharpMaskingImageFilter:
public itk::ImageToImageFilter<TImage, TImage>
{

public:

  /** Standard class typedefs. */
  typedef UnsharpMaskingImageFilter Self;
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

  itkGetConstMacro(BlurVariance, double);
  itkSetMacro(BlurVariance, double);

  itkGetConstMacro(OutputMinimum, double);
  itkSetMacro(OutputMinimum, double);

protected:

  UnsharpMaskingImageFilter();
  ~UnsharpMaskingImageFilter() { }

  void GenerateData();

private:

  double m_BlurVariance;

  double m_OutputMinimum;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "UnsharpMaskingImageFilter.txx"
#endif

#endif
