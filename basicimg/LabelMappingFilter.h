
////////////////////////////////////////////////////////////////////////////////
//
// Transform labels using user defined mappings
//
// If a label in domain has no defined mapping, it is left as is
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 5/2004

#ifndef _LabelMappingFilter_h
#define _LabelMappingFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"

#include "vnl/vnl_vector.h"

template <class TImage>
class LabelMappingFilter:
  public itk::ImageToImageFilter<TImage, TImage>
{

public:

  /** Standard class typedefs. */
  typedef LabelMappingFilter Self;
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

  typedef vnl_vector<double> VectorType;

  itkGetConstMacro(InputVector, VectorType);
  itkSetMacro(InputVector, VectorType);

  itkGetConstMacro(OutputVector, VectorType);
  itkSetMacro(OutputVector, VectorType);

protected:

  LabelMappingFilter();
  ~LabelMappingFilter() { }

  void GenerateData();

private:

  VectorType m_InputVector;
  VectorType m_OutputVector;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "LabelMappingFilter.txx"
#endif

#endif
