
#ifndef _OrientImageModifier_txx
#define _OrientImageModifier_txx

#include "itkExceptionObject.h"
#include "itkFixedArray.h"
#include "itkFlipImageFilter.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkImageRegionIterator.h"

#include "OrientImageModifier.h"

#include <ctype.h>

template <class TImage>
OrientImageModifier<TImage>
::OrientImageModifier()
{

  std::string rai = "RAI";

  this->SetSourceOrientation(rai);
  this->SetTargetOrientation(rai);

}

template <class TImage>
OrientImageModifier<TImage>
::~OrientImageModifier()
{

}

template <class TImage>
bool
OrientImageModifier<TImage>
::IsValidOrientation(std::string s)
{

  try
  {
    OrientInfo dummy = this->_GetOrientInfo(s);
  }
  catch (itk::ExceptionObject& err)
  {
    return false;
  }

  return true;

}

template <class TImage>
typename OrientImageModifier<TImage>::OrientInfo
OrientImageModifier<TImage>
::_GetOrientInfo(std::string& s)
{

  if (s.length() != 3)
    itkExceptionMacro(<< "Invalid orientation: " << s);

  OrientInfo info;

  for (unsigned int i = 0; i < 3; i++)
  {
    char ch = toupper(s[i]);

    switch (ch)
    {
      case 'R':
        info.axes[i] = 0;
        info.flips[i] = false;
        break;
      case 'L':
        info.axes[i] = 0;
        info.flips[i] = true;
        break;

      case 'A':
        info.axes[i] = 1;
        info.flips[i] = false;
        break;
      case 'P':
        info.axes[i] = 1;
        info.flips[i] = true;
        break;

      case 'I':
        info.axes[i] = 2;
        info.flips[i] = false;
        break;
      case 'S':
        info.axes[i] = 2;
        info.flips[i] = true;
        break;

      default:
       itkExceptionMacro(<< "Invalid orientation: " << s);
    }
  }

  // Verify that there are no duplicates in axes encoding
  for (unsigned int i = 0; i < 3; i++)
    for (unsigned int j = i+1; j < 3; j++)
      if (info.axes[i] == info.axes[j])
        itkExceptionMacro(<< "Invalid orientation: " << s); 

  return info;

}

template <class TImage>
void
OrientImageModifier<TImage>
::SetSourceOrientation(std::string s)
{

  m_SourceOrientString = s;

}

template <class TImage>
void
OrientImageModifier<TImage>
::SetTargetOrientation(std::string s)
{

  m_TargetOrientString = s;

}

template <class TImage>
void
OrientImageModifier<TImage>
::Modify(ImageType* img)
{

  if (img->GetImageDimension() != 3)
    itkExceptionMacro(<< "Only supports 3D image");

  // Convert strings to orient info
  OrientInfo sourceInfo = this->_GetOrientInfo(m_SourceOrientString);
  OrientInfo targetInfo = this->_GetOrientInfo(m_TargetOrientString);

  // Get permutation and flips
  itk::FixedArray<unsigned int, 3> permutation;
  itk::FixedArray<bool, 3> flips;

  for (unsigned int j = 0; j < 3; j++)
  {
    unsigned char target_ax = targetInfo.axes[j];
    bool target_fl = targetInfo.flips[j];

    for (unsigned int i = 0; i < 3; i++)
    {
      unsigned char source_ax = sourceInfo.axes[i];
      bool source_fl = sourceInfo.flips[i];

      if (target_ax == source_ax)
      {
        permutation[j] = i;
        if (target_fl != source_fl)
          flips[j] = true;
        else
          flips[j] = false;
      }
    }
  }

  // Get oriented image
  typedef itk::PermuteAxesImageFilter<ImageType> PermuteFilterType;
  typename PermuteFilterType::Pointer permfilter = PermuteFilterType::New();

  permfilter->SetOrder(permutation);
  permfilter->SetInput(img);
  permfilter->Update();

  typedef itk::FlipImageFilter<ImageType> FlipFilterType;
  typename FlipFilterType::Pointer flipfilter = FlipFilterType::New();

  flipfilter->SetFlipAxes(flips);
  flipfilter->SetInput(permfilter->GetOutput());
  flipfilter->Update();

  ImagePointer tmp = flipfilter->GetOutput();

  // Set origin to zero
  unsigned int dim = ImageType::GetImageDimension();

  float* v = new float[dim];
  for (unsigned int i = 0; i < dim; i++)
    v[i] = 0;

  tmp->SetOrigin(v);

  delete [] v;

  // Modify input space specification (and reallocate memory)
  img->Initialize();
  img->SetRegions(tmp->GetLargestPossibleRegion());
  img->Allocate();
  img->SetOrigin(tmp->GetOrigin());
  img->SetSpacing(tmp->GetSpacing());

  // Copy result to input
  typedef itk::ImageRegionIterator<ImageType> IteratorType;

  IteratorType inIt(tmp, tmp->GetLargestPossibleRegion());
  IteratorType outIt(img, img->GetLargestPossibleRegion());

  inIt.GoToBegin();
  outIt.GoToBegin();
  while (!outIt.IsAtEnd())
  {
    outIt.Set(inIt.Get());
    ++inIt;
    ++outIt;
  }

}

#endif
