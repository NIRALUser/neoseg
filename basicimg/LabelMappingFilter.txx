
#ifndef _LabelMappingFilter_txx
#define _LabelMappingFilter_txx

#include "LabelMappingFilter.h"

#include "itkImageRegionIterator.h"

template <class TImage>
LabelMappingFilter<TImage>
::LabelMappingFilter()
{

  m_InputVector = VectorType(0);
  m_OutputVector = VectorType(0);

}

template <class TImage>
void
LabelMappingFilter<TImage>
::GenerateData()
{

  if (this->GetInput() == 0)
    return;
 
  if (m_InputVector.size() != m_OutputVector.size())
    itkExceptionMacro(<< "Mapping vectors size mismatch");

  ImageRegionType region = this->GetInput()->GetLargestPossibleRegion();

  this->GetOutput()->SetRegions(region);
  this->GetOutput()->Allocate();
  this->GetOutput()->SetOrigin(this->GetInput()->GetOrigin());
  this->GetOutput()->SetSpacing(this->GetInput()->GetSpacing());

  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
  typedef itk::ImageRegionIterator<ImageType> IteratorType;

  ConstIteratorType inputIt(this->GetInput(), region);
  IteratorType outputIt(this->GetOutput(), region);

  // Map labels
  inputIt.GoToBegin();
  outputIt.GoToBegin();
  while(!outputIt.IsAtEnd())
  {
    double x = inputIt.Get();

    double mapx = x;
    for (unsigned int k = 0; k < m_InputVector.size(); k++)
      if (x == m_InputVector[k])
      {
        mapx = m_OutputVector[k];
        break;
      }

    outputIt.Set(mapx);

    ++inputIt;
    ++outputIt;
  }

}

#endif
