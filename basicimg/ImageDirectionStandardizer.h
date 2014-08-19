
#ifndef _ImageDirectionStandardizer_h
#define _ImageDirectionStandardizer_h

#include "itkImage.h"

#include "itkSpatialOrientation.h"

#include <string>

template <class TImage>
class ImageDirectionStandardizer: public itk::Object
{

public:

  typedef ImageDirectionStandardizer Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self);

  typedef TImage ImageType;
  typedef typename TImage::Pointer ImagePointer;
  typedef typename TImage::DirectionType ImageDirectionType;
  typedef typename TImage::IndexType ImageIndexType;
  typedef typename TImage::PointType ImagePointType;
  typedef typename TImage::SizeType ImageSizeType;
  typedef typename TImage::SpacingType ImageSpacingType;

  typedef itk::SpatialOrientation::ValidCoordinateOrientationFlags
    CoordinateOrientationCode;

  void SetTargetDirectionFromString(ImageType* img, std::string& s);

  // Standardize an image given string option dir
  // dir = file -> use dir info from file
  // dir = RAI/ASR/etc -> override direction using appropriate matrix from code
  ImagePointer Standardize(ImageType* img, std::string& dir);

protected:

  ImageDirectionStandardizer();
  ~ImageDirectionStandardizer();

  CoordinateOrientationCode _ParseOrientationString(std::string& s);

  ImageDirectionType _GetDirectionFromString(std::string& s);

  ImageDirectionType m_TargetImageOrientation;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "ImageDirectionStandardizer.txx"
#endif

#endif
