
#ifndef _TrimmedRescaleImageModifier_txx
#define _TrimmedRescaleImageModifier_txx

#include "TrimmedRescaleImageModifier.h"

#include "itkImageRegionIterator.h"

template <class TImage>
TrimmedRescaleImageModifier<TImage>
::TrimmedRescaleImageModifier()
{
  m_TrimFraction = 0.01;

  m_OutputMinimum = 0.0;
  m_OutputMaximum = 1.0;

  m_TrimAboveFlag = true;
  m_TrimBelowFlag = true;
}

template <class TImage>
void
TrimmedRescaleImageModifier<TImage>
::Rescale(ImagePointer img)
{

  double outRange = m_OutputMaximum - m_OutputMinimum;
  if (outRange <= 0)
    itkExceptionMacro(<< "Invalid output range: " << outRange);

  typedef itk::ImageRegionIterator<ImageType> IteratorType;

  IteratorType it(img, img->GetLargestPossibleRegion());

  // Get image range
  it.GoToBegin();
  double min = it.Get();
  double max = min;

  it.GoToBegin();
  while(!it.IsAtEnd())
  {
    double v = it.Get();

    if (v < min)
      min = v;
    if (v > max)
      max = v;

    ++it;
  }

  double range = max - min;
  if (range <= 0)
    itkExceptionMacro(<< "Invalid range: " << range);

  // Trim range
  double trimAmount = m_TrimFraction * range;
  double trimRange = range;

  if (m_TrimBelowFlag)
  {
    min += trimAmount;
    trimRange -= trimAmount;
  }

  if (m_TrimAboveFlag)
  {
    max -= trimAmount;
    trimRange -= trimAmount;
  }

  if (trimRange <= 0)
    itkExceptionMacro(<< "Invalid trimmed range: " << trimRange);

  it.GoToBegin();
  while(!it.IsAtEnd())
  {
    double u = (it.Get() - min) / trimRange;

    if (u < 0)
      u = 0;
    if (u > 1)
      u = 1;

    it.Set((u * outRange) + m_OutputMinimum);

    ++it;
  }

}

#endif
