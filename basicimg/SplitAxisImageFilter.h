
////////////////////////////////////////////////////////////////////////////////
//
// Split image along an axis (e.g. dual echo -> PD + T2)
// The even and odd indices form separate images
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 12/2005

#ifndef _SplitAxisImageFilter_h
#define _SplitAxisImageFilter_h

template <class TImage>
class SplitAxisImageFilter:
//public itk::ImageToImageFilter<TImage, TImage>
public itk::Object
{

public:

  /** Standard class typedefs. */
  typedef SplitAxisImageFilter Self;
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

  void  SetInput(ImageType* img);

  ImageType* GetEvenImage();
  ImageType* GetOddImage();

  itkGetConstMacro(SplitAxis, unsigned int);
  //itkSetMacro(SplitAxis, unsigned int);
  void SetSplitAxis(unsigned int i)
  { m_SplitAxis = i; m_Modified = true; }

  void Update();

protected:

  SplitAxisImageFilter();
  ~SplitAxisImageFilter() { }

private:

  bool m_Modified;

  unsigned int m_SplitAxis;

  ImagePointer m_Input;

  ImagePointer m_OddImage;
  ImagePointer m_EvenImage;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "SplitAxisImageFilter.txx"
#endif

#endif
