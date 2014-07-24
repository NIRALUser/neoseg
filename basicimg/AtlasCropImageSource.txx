
#ifndef _AtlasCropImageSource_txx
#define _AtlasCropImageSource_txx

#include "itkImageRegionIteratorWithIndex.h"

#include "AtlasCropImageSource.h"

template <class TInputImage, class TProbabilityImage>
AtlasCropImageSource<TInputImage, TProbabilityImage>
::AtlasCropImageSource()
{

  m_Padding = 8.0;

  m_LowerBound.Fill(0);
  m_UpperBound.Fill(0);

  m_OriginalSize.Fill(0);

}

template <class TInputImage, class TProbabilityImage>
bool
AtlasCropImageSource<TInputImage, TProbabilityImage>
::CheckBounds()
{

  for (unsigned int i = 0; i < ImageDimension; i++)
    if (m_LowerBound[i] > m_UpperBound[i])
      return false;

  InputImageSizeType croppedSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
    croppedSize[i] = m_UpperBound[i] - m_LowerBound[i] + 1;

  for (unsigned int i = 0; i < ImageDimension; i++)
    if (croppedSize[i] > m_OriginalSize[i])
      return false;

  return true;
}

template <class TInputImage, class TProbabilityImage>
void
AtlasCropImageSource<TInputImage, TProbabilityImage>
::UseProbabilities(ProbabilityImageList probs)
{

  if (probs.GetSize() == 0)
    itkExceptionMacro(<< "Need at least one class probability image");

  ProbabilityImageSizeType size =
    probs[0]->GetLargestPossibleRegion().GetSize();

  ProbabilityImageSpacingType spacing =
    probs[0]->GetSpacing();

  // Make sure all class probabilities have the same space
  for (unsigned int i = 1; i < probs.GetSize(); i++)
  {
    ProbabilityImageSizeType othersize =
      probs[i]->GetLargestPossibleRegion().GetSize();
    if (size != othersize)
      itkExceptionMacro(<< "Probability list size mismatch");
  }

  // Transform padding to voxel counts
  InputImageOffsetType padding;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    padding[i] = (unsigned int)floor(m_Padding / spacing[i] + 0.5);
  }

  // Make sure padding is sensible
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    if (size[i] <= padding[i])
      itkExceptionMacro(
        << "Bounding box padding larger than or equal to image size");
  }

  // Save size info
  m_OriginalSize = size;

  // Initial bounds: whole image
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    m_LowerBound[i] = size[i]-1; 
    m_UpperBound[i] = 0;
  }

  // Go through image and update bounds
  typedef itk::ImageRegionIteratorWithIndex<ProbabilityImageType>
    IteratorType;

  IteratorType it(probs[0], probs[0]->GetLargestPossibleRegion());

  it.GoToBegin();
  while (!it.IsAtEnd())
  {
    ProbabilityImageIndexType ind = it.GetIndex();

    double sumProb = 0;
    for (unsigned int i = 0; i < probs.GetSize(); i++)
      sumProb += probs[i]->GetPixel(ind);

    if (sumProb > 0)
    {
      for (unsigned int i = 0; i < ImageDimension; i++)
      {
        if (ind[i] < m_LowerBound[i])
          m_LowerBound[i] = ind[i];
        if (ind[i] > m_UpperBound[i])
          m_UpperBound[i] = ind[i];
      }
    }

    ++it;
  }

  // Enlarge/pad bounding box
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    if (m_LowerBound[i] > padding[i])
      m_LowerBound[i] -= padding[i];
    else
      m_LowerBound[i] = 0;

    if (m_UpperBound[i] < (size[i]-1-padding[i]))
      m_UpperBound[i] += padding[i];
    else
      m_UpperBound[i] = size[i]-1;
  }

  //m_CropInfo.offset =
  //m_CropInfo.size =

}

template <class TInputImage, class TProbabilityImage>
typename AtlasCropImageSource<TInputImage, TProbabilityImage>
  ::InputImagePointer
AtlasCropImageSource<TInputImage, TProbabilityImage>
::Restore(InputImagePointer img)
{

  if (!this->CheckBounds())
  {
    itkExceptionMacro(<< "Invalid bounds");
  }

  InputImageSizeType size = img->GetLargestPossibleRegion().GetSize();

  InputImageSizeType croppedSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
    croppedSize[i] = m_UpperBound[i] - m_LowerBound[i] + 1;

  if (size != croppedSize)
    itkExceptionMacro(<< "Input size does not match size of cropped image");

  InputImageRegionType region;
  region.SetSize(m_OriginalSize);

  InputImagePointer output = InputImageType::New();
  output->SetRegions(region);
  output->Allocate();
  output->SetOrigin(img->GetOrigin());
  output->SetSpacing(img->GetSpacing());
  output->FillBuffer(0);

  InputImageOffsetType offt;
  for (unsigned int i = 0; i < ImageDimension; i++)
    offt[i] = m_LowerBound[i];

  typedef itk::ImageRegionIteratorWithIndex<InputImageType> IteratorType;

  IteratorType inIt(img, img->GetLargestPossibleRegion());

  inIt.GoToBegin();
  while (!inIt.IsAtEnd())
  {
    InputImageIndexType ind = inIt.GetIndex();

    output->SetPixel(ind + offt, inIt.Get());

    ++inIt;
  }

  return output;

}


template <class TInputImage, class TProbabilityImage>
typename AtlasCropImageSource<TInputImage, TProbabilityImage>
  ::InputImagePointer
AtlasCropImageSource<TInputImage, TProbabilityImage>
::Crop(InputImagePointer img)
{

  if (!this->CheckBounds())
  {
    itkExceptionMacro(<< "Invalid bounds");
  }

  InputImageSizeType size = img->GetLargestPossibleRegion().GetSize();

  if (size != m_OriginalSize)
    itkExceptionMacro(<< "Input size does not match size of probability images");

  InputImageSizeType croppedSize;
  for (unsigned int i = 0; i < ImageDimension; i++)
    croppedSize[i] = m_UpperBound[i] - m_LowerBound[i] + 1;

  InputImageRegionType cropRegion;
  cropRegion.SetSize(croppedSize);

  InputImagePointer output = InputImageType::New();

  output->SetRegions(cropRegion);
  output->Allocate();
  output->SetOrigin(img->GetOrigin());
  output->SetSpacing(img->GetSpacing());

  InputImageOffsetType offt;
  for (unsigned int i = 0; i < ImageDimension; i++)
    offt[i] = m_LowerBound[i];

  typedef itk::ImageRegionIteratorWithIndex<InputImageType> IteratorType;

  IteratorType outIt(output, cropRegion);

  outIt.GoToBegin();
  while (!outIt.IsAtEnd())
  {
    InputImageIndexType ind = outIt.GetIndex();

    outIt.Set(img->GetPixel(ind+offt));

    ++outIt;
  }

  return output;

}

#endif
