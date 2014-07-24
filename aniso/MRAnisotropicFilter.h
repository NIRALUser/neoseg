
////////////////////////////////////////////////////////////////////////////////
//
// Non-linear anisotropic diffusion filtering for multi-channel 3D MR data
//
// Input images are directly modified
//
// Gerig, G., K¨ubler, O., Kikinis, R., Jolesz, F.: Nonlinear anisotropic
// filtering of MRI data. IEEE TMI 11 (1992) 221-232
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 3/2004

#ifndef _MRAnisotropicFilter_h
#define _MRAnisotropicFilter_h

#include "itkImage.h"
#include "itkObject.h"

#include "DynArray.h"

template <class TInputImage>
class MRAnisotropicFilter : public itk::Object
{

public:

  /** Standard class typedefs. */
  typedef MRAnisotropicFilter Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** The dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  // Image types
  typedef TInputImage InputImageType;
  typedef typename TInputImage::Pointer InputImagePointer;
  typedef typename TInputImage::IndexType InputImageIndexType;
  typedef typename TInputImage::OffsetType InputImageOffsetType;
  typedef typename TInputImage::PixelType InputImagePixelType;
  typedef typename TInputImage::SizeType InputImageSizeType;
  typedef typename TInputImage::SpacingType InputImageSpacingType;
  typedef typename TInputImage::RegionType InputImageRegionType;

  typedef itk::Image<float, 3> InternalImageType;
  typedef InternalImageType::Pointer InternalImagePointer;
  typedef InternalImageType::IndexType InternalImageIndexType;
  typedef InternalImageType::OffsetType InternalImageOffsetType;
  typedef InternalImageType::PixelType InternalImagePixelType;
  typedef InternalImageType::RegionType InternalImageRegionType;
  typedef InternalImageType::SizeType InternalImageSizeType;

  itkGetMacro(NumberOfIterations, unsigned);
  itkSetMacro(NumberOfIterations, unsigned);

  itkGetMacro(TimeStep, double);
  itkSetMacro(TimeStep, double);

  void Filter(DynArray<InputImagePointer> images);

protected:

  MRAnisotropicFilter();
  ~MRAnisotropicFilter();

  void CheckInput();

  void CopyInput();

  double ComputeAverageGradientNorm();

  void Step();

private:

  double m_K;

  unsigned m_NumberOfIterations;

  double m_TimeStep;

  InternalImagePointer m_DiffusionImage;

  DynArray<InputImagePointer> m_InputImages;

  DynArray<InternalImagePointer> m_InternalImages;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "MRAnisotropicFilter.txx"
#endif

#endif
