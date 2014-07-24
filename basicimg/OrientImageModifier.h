
#ifndef _OrientImageModifier_h
#define _OrientImageModifier_h

#include "itkImage.h"

#include <string>

template <class TImage>
class OrientImageModifier: public itk::Object
{

public:

  typedef OrientImageModifier Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self);

  typedef TImage ImageType;
  typedef typename TImage::Pointer ImagePointer;
  typedef typename TImage::IndexType ImageIndexType;
  typedef typename TImage::SizeType ImageSizeType;
  typedef typename TImage::SpacingType ImageSpacingType;

  // Orient encoding
  // R/L -> 0, A/P -> 1, I/S ->2
  // R is false, L is true (X/Y -> false/true)
  typedef struct 
  {
    unsigned char axes[3];
    bool flips[3];
  }
  OrientInfo;

  bool IsValidOrientation(std::string s);

  void SetSourceOrientation(std::string s);
  void SetTargetOrientation(std::string s);

  void Modify(ImageType* img);

protected:

  OrientImageModifier();
  ~OrientImageModifier();

  OrientInfo _GetOrientInfo(std::string& s);

  std::string m_SourceOrientString;
  std::string m_TargetOrientString;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "OrientImageModifier.txx"
#endif

#endif
