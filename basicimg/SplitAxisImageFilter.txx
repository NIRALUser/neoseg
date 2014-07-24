
#ifndef _SplitAxisImageFilter_txx
#define _SplitAxisImageFilter_txx

#include "SplitAxisImageFilter.h"

#include "itkImageRegionConstIteratorWithIndex.h"

template <class TImage>
SplitAxisImageFilter<TImage>
::SplitAxisImageFilter()
{

  m_SplitAxis = 0;

  m_Modified = false;

}

template <class TImage>
void
SplitAxisImageFilter<TImage>
::SetInput(ImageType* img)
{
  m_Modified = true;

  m_Input = img;
}

template <class TImage>
typename SplitAxisImageFilter<TImage>::ImageType*
SplitAxisImageFilter<TImage>
::GetEvenImage()
{
  if (m_Modified)
    this->Update();

  return m_EvenImage;
}

template <class TImage>
typename SplitAxisImageFilter<TImage>::ImageType*
SplitAxisImageFilter<TImage>
::GetOddImage()
{
  if (m_Modified)
    this->Update();

  return m_OddImage;
}

template <class TImage>
void
SplitAxisImageFilter<TImage>
::Update()
{

  if (m_SplitAxis >= Self::ImageDimension)
    itkExceptionMacro(<< "Invalid axis: " << m_SplitAxis);

  ImageRegionType region = m_Input->GetLargestPossibleRegion();

  ImageSizeType size = region.GetSize();

  unsigned int axisSize_even = (size[m_SplitAxis]+1) / 2;
  unsigned int axisSize_odd = size[m_SplitAxis] / 2;

  ImageSizeType halfSize = size;

  ImageRegionType halfRegion = region;

  halfSize[m_SplitAxis] = axisSize_even;
  halfRegion.SetSize(halfSize);
  m_EvenImage = ImageType::New();
  m_EvenImage->SetRegions(halfRegion);
  m_EvenImage->Allocate();
  m_EvenImage->SetOrigin(m_Input->GetOrigin());
  m_EvenImage->SetSpacing(m_Input->GetSpacing());

  halfSize[m_SplitAxis] = axisSize_odd;
  halfRegion.SetSize(halfSize);
  m_OddImage = ImageType::New();
  m_OddImage->SetRegions(halfRegion);
  m_OddImage->Allocate();
  m_OddImage->SetOrigin(m_Input->GetOrigin());
  m_OddImage->SetSpacing(m_Input->GetSpacing());

  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ConstIteratorType;

  ConstIteratorType inputIt(m_Input, region);

  inputIt.GoToBegin();
  while(!inputIt.IsAtEnd())
  {
    ImageIndexType ind = inputIt.GetIndex();

    ImagePixelType v = inputIt.Get();

    if ((ind[m_SplitAxis] % 2) == 0)
    {
      ImageIndexType ind_even = ind;
      ind_even[m_SplitAxis] = ind[m_SplitAxis] / 2;
      m_EvenImage->SetPixel(ind_even, v);
    }
    else
    {
      ImageIndexType ind_odd = ind;
      ind_odd[m_SplitAxis] = (ind[m_SplitAxis]-1) / 2;
      m_OddImage->SetPixel(ind_odd, v);
    }

    ++inputIt;
  }

  m_Modified = false;

}

#endif
