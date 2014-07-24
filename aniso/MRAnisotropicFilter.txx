
#ifndef _MRAnisotropicFilter_txx
#define _MRAnisotropicFilter_txx

#include "MRAnisotropicFilter.h"

#include <float.h>
#include <math.h>


template <class TInputImage>
MRAnisotropicFilter<TInputImage>
::MRAnisotropicFilter()
{
  m_DiffusionImage = 0;

  m_K = 1.0;

  m_NumberOfIterations = 10;

  m_TimeStep = 0.01;
}

template <class TInputImage>
MRAnisotropicFilter<TInputImage>
::~MRAnisotropicFilter()
{

}

template <class TInputImage>
void
MRAnisotropicFilter<TInputImage>
::CheckInput()
{

  itkDebugMacro(<< "CheckInput");

  if (m_InputImages.GetSize() < 1)
    itkExceptionMacro(<< "Need at least one channel");

  if (m_InputImages[0].IsNull())
    itkExceptionMacro(<< "Image 0 is NULL");

  if (m_InputImages[0]->GetImageDimension() != 3)
    itkExceptionMacro(<< "Invalid image dimension: only supports 3D images");

  double maxdt = 0.5 / pow(2.0, 3.0);
  if (m_TimeStep > maxdt)
    itkExceptionMacro(<< "Time step is too high, should be below " << maxdt);

  InputImageSizeType size = m_InputImages[0]->GetLargestPossibleRegion().GetSize();

  InputImageSpacingType spacing = m_InputImages[0]->GetSpacing();

  for (unsigned int i = 1; i < m_InputImages.GetSize(); i++)
  {
    if (m_InputImages[i].IsNull())
      itkExceptionMacro(<< "Image " << i << " is NULL");

    InputImageSizeType size_i =
      m_InputImages[i]->GetLargestPossibleRegion().GetSize();
    if (size_i != size)
      itkExceptionMacro(<< "Image sizes do not match");

    InputImageSpacingType spacing_i = m_InputImages[i]->GetSpacing();
    if (spacing_i != spacing)
      itkExceptionMacro(<< "Voxel spacings do not match");
  }

}

template <class TInputImage>
void
MRAnisotropicFilter<TInputImage>
::CopyInput()
{

  itkDebugMacro(<< "Copying image to internal structure...");

  // Clear the internal data
  m_InternalImages.Clear();

  InputImageRegionType region = m_InputImages[0]->GetLargestPossibleRegion();

  InputImageIndexType ind;

  InputImageSizeType size = region.GetSize();

  for (unsigned int i = 0; i < m_InputImages.GetSize(); i++)
  {
    InternalImagePointer img = InternalImageType::New();
    img->SetRegions(region);
    img->Allocate();
    img->SetOrigin(m_InputImages[0]->GetOrigin());
    img->SetSpacing(m_InputImages[0]->GetSpacing());

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          img->SetPixel(ind, (InternalImagePixelType)
            m_InputImages[i]->GetPixel(ind));
        }

    m_InternalImages.Append(img);
  }

}

template <class TInputImage>
double
MRAnisotropicFilter<TInputImage>
::ComputeAverageGradientNorm()
{
  itkDebugMacro(<< "Computing average gradient norm...");

  InputImageRegionType region = m_InputImages[0]->GetLargestPossibleRegion();

  InputImageIndexType ind;

  InputImageSizeType size = region.GetSize();

  InputImageSpacingType spacing = m_InputImages[0]->GetSpacing();

  double dx = spacing[0] + FLT_EPSILON;
  double dy = spacing[1] + FLT_EPSILON;
  double dz = spacing[2] + FLT_EPSILON;

  double dx2 = 2.0 * dx;
  double dy2 = 2.0 * dy;
  double dz2 = 2.0 * dz;

  InternalImageOffsetType xofft = {{1, 0, 0}};
  InternalImageOffsetType yofft = {{0, 1, 0}};
  InternalImageOffsetType zofft = {{0, 0, 1}};

  double sumGradNorm = 0;
  unsigned int countNonZero = 0;
  for (ind[2] = 1; ind[2] < (size[2]-1); ind[2]++)
    for (ind[1] = 1; ind[1] < (size[1]-1); ind[1]++)
      for (ind[0] = 1; ind[0] < (size[0]-1); ind[0]++)
      {
        double sumSquaredNorm = 0;

        for (unsigned int i = 0; i < m_InternalImages.GetSize(); i++)
        {
          InternalImagePointer img = m_InternalImages[i];
          double gx = img->GetPixel(ind+xofft) - img->GetPixel(ind-xofft);
          gx /= dx2;
          double gy = img->GetPixel(ind+yofft) - img->GetPixel(ind-yofft);
          gy /= dy2;
          double gz = img->GetPixel(ind+zofft) - img->GetPixel(ind-zofft);
          gz /= dz2;
          sumSquaredNorm += gx*gx + gy*gy + gz*gz;
        }

        double gradNorm = sqrt(sumSquaredNorm);

        if (gradNorm > 0)
          countNonZero++;

        sumGradNorm += gradNorm;
      }

  return sumGradNorm / countNonZero;

}

template <class TInputImage>
void
MRAnisotropicFilter<TInputImage>
::Step()
{

  itkDebugMacro(<< "Step");

  InputImageRegionType region = m_InputImages[0]->GetLargestPossibleRegion();

  InputImageIndexType ind;

  InputImageSizeType size = region.GetSize();

  InputImageSpacingType spacing = m_InputImages[0]->GetSpacing();

  double dx = spacing[0] + FLT_EPSILON;
  double dy = spacing[1] + FLT_EPSILON;
  double dz = spacing[2] + FLT_EPSILON;

  double dx2 = 2.0 * dx;
  double dy2 = 2.0 * dy;
  double dz2 = 2.0 * dz;

  double dxdx = dx*dx;
  double dydy = dy*dy;
  double dzdz = dz*dz;

  InternalImageOffsetType xofft = {{1, 0, 0}};
  InternalImageOffsetType yofft = {{0, 1, 0}};
  InternalImageOffsetType zofft = {{0, 0, 1}};

  double kSquared = m_K*m_K;

  itkDebugMacro(<< "Computing diffusion terms");
  for (ind[2] = 1; ind[2] < (size[2]-1); ind[2]++)
    for (ind[1] = 1; ind[1] < (size[1]-1); ind[1]++)
      for (ind[0] = 1; ind[0] < (size[0]-1); ind[0]++)
      {
        double sumSquaredNorm = 0;

        for (unsigned int i = 0; i < m_InternalImages.GetSize(); i++)
        {
          InternalImagePointer img = m_InternalImages[i];

// Half-step derivative
// dI = I(x+dx/2) - I(x-dx/2)
//    = (I(x+dx)+I(x)) / 2 - (I(x-dx)+I(x)) / 2
//    = (I(x+dx) - I(x-dx)) / 2

          double gx = img->GetPixel(ind+xofft) - img->GetPixel(ind-xofft);
          gx /= dx2;
          double gy = img->GetPixel(ind+yofft) - img->GetPixel(ind-yofft);
          gy /= dy2;
          double gz = img->GetPixel(ind+zofft) - img->GetPixel(ind-zofft);
          gz /= dz2;

          sumSquaredNorm += gx*gx + gy*gy + gz*gz;
        }

#if 1
        double dTerm = exp(-sumSquaredNorm / kSquared);
#else
        double r = sqrt(sumSquaredNorm) / m_K;
        double dTerm = 1.0 / (1.0 + r*r);
#endif

        m_DiffusionImage->SetPixel(ind, (float)dTerm);
      }

  itkDebugMacro(<< "Flow loop");
  for (ind[2] = 1; ind[2] < (size[2]-1); ind[2]++)
    for (ind[1] = 1; ind[1] < (size[1]-1); ind[1]++)
      for (ind[0] = 1; ind[0] < (size[0]-1); ind[0]++)
      {
        // Get the diffusion terms
        double diff_center = m_DiffusionImage->GetPixel(ind);

        double diffx_fwd =
          (m_DiffusionImage->GetPixel(ind+xofft) + diff_center) / 2.0;
        double diffx_bwd =
          (m_DiffusionImage->GetPixel(ind-xofft) + diff_center) / 2.0;

        double diffy_fwd =
          (m_DiffusionImage->GetPixel(ind+yofft) + diff_center) / 2.0;
        double diffy_bwd =
          (m_DiffusionImage->GetPixel(ind-yofft) + diff_center) / 2.0;

        double diffz_fwd =
          (m_DiffusionImage->GetPixel(ind+zofft) + diff_center) / 2.0;
        double diffz_bwd =
          (m_DiffusionImage->GetPixel(ind-zofft) + diff_center) / 2.0;

        for (unsigned int i = 0; i < m_InternalImages.GetSize(); i++)
        {
          InternalImagePointer img = m_InternalImages[i];
          InputImagePointer img0 = m_InputImages[i];

          double phi = 0;

          // Flow along x
          phi += diffx_fwd
            * (img->GetPixel(ind+xofft)-img->GetPixel(ind)) / dxdx;
          phi -= diffx_bwd
            * (img->GetPixel(ind)-img->GetPixel(ind-xofft)) / dxdx;

          // Flow along y
          phi += diffy_fwd
            * (img->GetPixel(ind+yofft)-img->GetPixel(ind)) / dydy;
          phi -= diffy_bwd
            * (img->GetPixel(ind)-img->GetPixel(ind-yofft)) / dydy;

          // Flow along z
          phi += diffz_fwd
            * (img->GetPixel(ind+zofft)-img->GetPixel(ind)) / dzdz;
          phi -= diffz_bwd
            * (img->GetPixel(ind)-img->GetPixel(ind-zofft)) / dzdz;

          // Nordstr\"om's biased anisotropic diffusion scheme
          //phi += (img0->GetPixel(ind) - img->GetPixel(ind));

          // Update intensity
          img->SetPixel(ind, (float)(img->GetPixel(ind) + m_TimeStep*phi));
        }
      }

}

template <class TInputImage>
void
MRAnisotropicFilter<TInputImage>
::Filter(DynArray<InputImagePointer> images)
{

  m_InputImages = images;

  this->CheckInput();

  this->CopyInput();

  InputImageRegionType region = m_InputImages[0]->GetLargestPossibleRegion();

  itkDebugMacro(<< "Allocating space for diffusion terms");
  m_DiffusionImage = InternalImageType::New();
  m_DiffusionImage->SetRegions(region);
  m_DiffusionImage->Allocate();
  m_DiffusionImage->SetOrigin(m_InputImages[0]->GetOrigin());
  m_DiffusionImage->SetSpacing(m_InputImages[0]->GetSpacing());

  // Initialize to zero
  InternalImageIndexType ind;
  InternalImageSizeType size = region.GetSize();
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        m_DiffusionImage->SetPixel(ind, 0);

  m_K = this->ComputeAverageGradientNorm() + 1e-20;

  itkDebugMacro(<< "Average gradient norm = " << m_K);

  for (unsigned int iter = 1; iter <= m_NumberOfIterations; iter++)
  {
    itkDebugMacro(<< "Iteration " << iter+1);

    // Update image for one time step
    this->Step();
  }

  // Delete the allocated diffusion image
  m_DiffusionImage = 0;

  // Set input image values to filtered values
  itkDebugMacro(<< "Modifying input intensities...");
  for (unsigned int i = 0; i < m_InputImages.GetSize(); i++)
  {
    InputImageIndexType ind;

    InputImagePointer img = m_InputImages[i];
    InternalImagePointer filt_img = m_InternalImages[i];

    for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
      for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
        for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        {
          img->SetPixel(ind, (InputImagePixelType)
            filt_img->GetPixel(ind));
        }

  }

  // Clear the internal data
  for (unsigned int i = 0; i < m_InternalImages.GetSize(); i++)
    m_InternalImages[i] = 0;
  m_InternalImages.Clear();

  // Remove reference count to original images
  for (unsigned int i = 0; i < m_InputImages.GetSize(); i++)
    m_InputImages[i] = 0;
  m_InputImages.Clear();

}


#endif
