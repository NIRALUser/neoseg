
////////////////////////////////////////////////////////////////////////////////
//
// Connected component
//
// Returns component labels, sorted according to size (1 is the largest)
//
////////////////////////////////////////////////////////////////////////////////

#ifndef _ConnectedComponentsFilter_h
#define _ConnectedComponentsFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"

template <class TInputImage, class TOutputImage>
class ConnectedComponentsFilter:
public itk::ImageToImageFilter<TInputImage, TOutputImage>
{

public:

  /** Standard class typedefs. */
  typedef ConnectedComponentsFilter Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(ConnectedComponentsFilter, ImageToImageFilter);

  /** The dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  // Image types
  typedef TInputImage InputImageType;
  typedef typename TInputImage::Pointer InputImagePointer;
  typedef typename TInputImage::IndexType InputImageIndexType;
  typedef typename TInputImage::OffsetType InputImageOffsetType;
  typedef typename TInputImage::PixelType InputImagePixelType;
  typedef typename TInputImage::SizeType InputImageSizeType;
  typedef typename TInputImage::RegionType InputImageRegionType;

  typedef TOutputImage OutputImageType;
  typedef typename TOutputImage::Pointer OutputImagePointer;
  typedef typename TOutputImage::IndexType OutputImageIndexType;
  typedef typename TOutputImage::OffsetType OutputImageOffsetType;
  typedef typename TOutputImage::PixelType OutputImagePixelType;
  typedef typename TOutputImage::SizeType OutputImageSizeType;
  typedef typename TOutputImage::RegionType OutputImageRegionType;

protected:

  ConnectedComponentsFilter() { }
  ~ConnectedComponentsFilter() { }

  void GenerateData();

};

#ifndef MU_MANUAL_INSTANTIATION
#include "ConnectedComponentsFilter.txx"
#endif

#endif
