
////////////////////////////////////////////////////////////////////////////////
//
// Customized itkSubsampledImageRegionIteratorWithIndex
//
// Supports offsets for doing subsampling
//
////////////////////////////////////////////////////////////////////////////////

#ifndef __SubsampledImageRegionIterator_h
#define __SubsampledImageRegionIterator_h

#include "itkImageRegionConstIteratorWithIndex.h"

template<typename TImage>
class SubsampledImageRegionIterator:
  public itk::ImageRegionConstIteratorWithIndex<TImage>
{
public:
  /** Standard class typedefs. */
  typedef SubsampledImageRegionIterator Self;
  typedef itk::ImageRegionConstIteratorWithIndex<TImage>  Superclass;
  
  /** Types inherited from the Superclass */
  typedef typename Superclass::IndexType              IndexType;
  typedef typename Superclass::IndexValueType         IndexValueType;
  typedef typename Superclass::SizeType               SizeType;
  typedef typename Superclass::SizeValueType          SizeValueType;
  typedef typename Superclass::OffsetType             OffsetType;
  typedef typename Superclass::OffsetValueType        OffsetValueType;
  typedef typename Superclass::RegionType             RegionType;
  typedef typename Superclass::ImageType              ImageType;
  typedef typename Superclass::PixelContainer         PixelContainer;
  typedef typename Superclass::PixelContainerPointer  PixelContainerPointer;
  typedef typename Superclass::InternalPixelType      InternalPixelType;
  typedef typename Superclass::PixelType              PixelType;
  typedef typename Superclass::AccessorType           AccessorType;

  /** Constructor establishes an iterator to walk a particular image and a
   * particular region of that image. */
  SubsampledImageRegionIterator(TImage *ptr, const RegionType& region, double sampleSpacing);

  /** Set the pixel value */
  void Set( const PixelType & value) const  
    { this->m_PixelAccessorFunctor.Set(*(const_cast<InternalPixelType *>(this->m_Position)),value); }

  /** Return a reference to the pixel 
   * This method will provide the fastest access to pixel
   * data, but it will NOT support ImageAdaptors. */
  PixelType & Value(void) 
    { return *(const_cast<InternalPixelType *>(this->m_Position)); }

  // Override the walk operators
  Self& operator++();
  Self& operator--();

protected:

  OffsetType m_Skips;
  
};

#ifndef MU_MANUAL_INSTANTIATION
#include "SubsampledImageRegionIterator.txx"
#endif

#endif 
