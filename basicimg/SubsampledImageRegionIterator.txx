
#ifndef _SubsampledImageRegionIterator_txx
#define _SubsampledImageRegionIterator_txx

#include "SubsampledImageRegionIterator.h"

#include <cmath>

template<typename TImage>
SubsampledImageRegionIterator<TImage>
::SubsampledImageRegionIterator(TImage *ptr, const RegionType& region, double sampleSpacing)
  : itk::ImageRegionConstIteratorWithIndex<TImage>(ptr, region) 
{
  typename TImage::SpacingType spacing = ptr->GetSpacing();

  for (unsigned int dim = 0; dim < TImage::ImageDimension; dim++)
  {
    long skip_dim =  (long)floor(sampleSpacing / spacing[dim]);
    if (skip_dim < 1)
      skip_dim = 1;
    m_Skips[dim] = skip_dim;
  }
}

template<class TImage>
SubsampledImageRegionIterator<TImage> &
SubsampledImageRegionIterator<TImage>
::operator++()
{   

// TODO:
// Snap to grid locations

// Then add skips[i]

  this->m_Remaining = false;
  for(unsigned int in = 0; in < TImage::ImageDimension; in++ )
  {

    // Go to next location in the subsampled grid
    long old_in = this->m_PositionIndex[in];
    this->m_PositionIndex[in] =
      ((old_in+m_Skips[in]) / m_Skips[in]) * m_Skips[in];

    if( this->m_PositionIndex[in] < this->m_EndIndex[in] )
    {
      long diff_in = this->m_PositionIndex[in] - old_in;
      this->m_Position += diff_in*this->m_OffsetTable[in];
      this->m_Remaining = true;
      break;
    }
    else
    {
      // Go back to first position for this dimension
      this->m_Position -= this->m_OffsetTable[in] * old_in;
      this->m_PositionIndex[in] = this->m_BeginIndex[in];
    }
  }

  if( !this->m_Remaining ) // It will not advance here otherwise
  {
    this->m_Position = this->m_End;
  }

  return *this;
}

template<class TImage>
SubsampledImageRegionIterator<TImage> &
SubsampledImageRegionIterator<TImage> 
::operator--()
{
 
  this->m_Remaining = false;
  for( unsigned int in=0; in<TImage::ImageDimension; in++ )
  {
    long old_in = this->m_PositionIndex[in];

    if (this->m_PositionIndex[in] > this->m_BeginIndex[in])
    { 
      this->m_PositionIndex[in] =
        ((old_in-m_Skips[in]) / m_Skips[in]) * m_Skips[in];
      if (this->m_PositionIndex[in] < 0)
        this->m_PositionIndex[in] = 0;
      long diff_in = old_in - this->m_PositionIndex[in];
      this->m_Position -= diff_in*this->m_OffsetTable[in];
      this->m_Remaining = true;
      break;
    }
    else
    {
      // At the beginning
      this->m_Position -= this->m_OffsetTable[in];
      this->m_PositionIndex[ in ] = this->m_EndIndex[ in ] - 1;
    }
  }

  if( !this->m_Remaining ) // It will not advance here otherwise
  {
    this->m_Position = this->m_End;
  }

  return *this;
}


#endif
