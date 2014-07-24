
#ifndef _UnsharpMaskingImageFilter_txx
#define _UnsharpMaskingImageFilter_txx

#include "UnsharpMaskingImageFilter.h"

#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageRegionIterator.h"

template <class TImage>
UnsharpMaskingImageFilter<TImage>
::UnsharpMaskingImageFilter()
{

  m_BlurVariance = 2.0;

  m_OutputMinimum = 0.0;

}

template <class TImage>
void
UnsharpMaskingImageFilter<TImage>
::GenerateData()
{

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType>
    BlurType;
  typename BlurType::Pointer blur = BlurType::New();

  blur->SetInput(this->GetInput());
  blur->SetMaximumError(1e-3);
  blur->SetMaximumKernelWidth(48);
  blur->SetVariance(m_BlurVariance);
  blur->Update();

  ImageRegionType region = this->GetInput()->GetLargestPossibleRegion();

  this->GetOutput()->SetRegions(region);
  this->GetOutput()->Allocate();

  this->GetOutput()->SetOrigin(this->GetInput()->GetOrigin());
  this->GetOutput()->SetSpacing(this->GetInput()->GetSpacing());

  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
  typedef itk::ImageRegionIterator<ImageType> IteratorType;

  ConstIteratorType inputIt(this->GetInput(), region);
  IteratorType blurIt(blur->GetOutput(), region);
  IteratorType outputIt(this->GetOutput(), region);

  inputIt.GoToBegin();
  blurIt.GoToBegin();
  outputIt.GoToBegin();
  while(!outputIt.IsAtEnd())
  {
    double x = inputIt.Get();
    double b = blurIt.Get();
    double y = 2*x - b;

    if (y < m_OutputMinimum)
      y = m_OutputMinimum;

    outputIt.Set(y);

    ++inputIt;
    ++blurIt;
    ++outputIt;
  }

}

#endif
