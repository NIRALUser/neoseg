
#ifndef _AtlasRegistrationMethod_txx
#define _AtlasRegistrationMethod_txx

#include "itkAffineTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

// MI registration module
#include "AtlasRegistrationMethod.h"
#include "PairRegistrationMethod.h"
#include "RegistrationParameters.h"

#include "Log.h"
#include "muFile.h"

#include "OrientImageModifier.h"

#include <fstream>

#define DO_BLUR 1

template <class TOutputPixel, class TProbabilityPixel>
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::AtlasRegistrationMethod()
{

  m_Suffix = "";

  m_OutputDirectory = "";

  m_TemplateFileName = "";

  m_ProbabilityFileNames.Clear();
  m_ImageFileNames.Clear();

  m_AtlasOrientation = "";
  m_ImageOrientations.Clear();

  m_Images.Clear();
  m_Probabilities.Clear();

  m_TemplateAffineTransform = AffineTransformType::New();
  m_AffineTransforms.Clear();

  m_AffineTransformReadFlags = FlagArrayType(1);
  m_AffineTransformReadFlags[0] = 0;

  m_UseNonLinearInterpolation = false;

  m_OutsideFOVCode = -65536;

  m_FOVMask = 0;

  m_DoneRegistration = false;
  m_DoneResample = false;

  m_WarpAtlas = false;

  m_WarpGridX = 5;
  m_WarpGridY = 5;
  m_WarpGridZ = 5;

  m_TemplateBSplineReadFlag = false;

  m_TemplateBSplineTransform = 0;

  m_Modified = false;

}

template <class TOutputPixel, class TProbabilityPixel>
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::~AtlasRegistrationMethod()
{

  m_ImageFileNames.Clear();
  m_ProbabilityFileNames.Clear();

  m_Probabilities.Clear();
  m_Images.Clear();
  m_AffineTransforms.Clear();

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::VerifyInitialization()
{

  if (m_ImageFileNames.GetSize() < 1)
    itkExceptionMacro(<< "No data images specified");

  if (m_ImageOrientations.GetSize() != m_ImageFileNames.GetSize())
    itkExceptionMacro(<< "Image - orientation info mismatch");

  /*
  // No atlas checks:
  // It's OK if we have no associated atlas files, then only do image-image
  // registrations
  if (m_TemplateFileName.length() == 0)
    itkExceptionMacro(<< "Template file name not specified");
  if (m_ProbabilityFileNames.GetSize() < 1)
    itkExceptionMacro(<< "No probability images specified");
  */

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::WriteParameters()
{

  itkDebugMacro(<< "Write parameters");

  if (!m_DoneRegistration)
    this->RegisterImages();

  std::string firststr =
    m_OutputDirectory + mu::get_name(m_ImageFileNames[0].c_str()) + std::string("_to_");

  std::string suffixstr;
  if (m_Suffix.length() != 0)
    suffixstr = std::string("_") + m_Suffix;
  else
    suffixstr = std::string("");

  if (m_TemplateFileName.length() != 0)
  {
    std::string name = mu::get_name(m_TemplateFileName.c_str());

    // Write recently computed affine transform
    if (m_AffineTransformReadFlags[0] == 0)
    {
      std::string affinefn =
        firststr + name + suffixstr + std::string(".affine");

      muLogMacro(<< "Writing " << affinefn << "...\n");

      PairRegistrationMethod<InternalImagePixelType>::
        WriteAffineTransform(affinefn.c_str(), m_TemplateAffineTransform);
    }

    // Write recently computed B-spline warp
    if (m_WarpAtlas && !m_TemplateBSplineReadFlag)
    {
      std::string warpfn =
        firststr + name + suffixstr + std::string(".bspline");

      muLogMacro(<< "Writing " << warpfn << "...\n");

      PairRegistrationMethod<InternalImagePixelType>::
        WriteBSplineTransform(warpfn.c_str(), m_TemplateBSplineTransform);
    }

  } // if template defined


  for (unsigned int i = 1; i < m_ImageFileNames.GetSize(); i++)
  {
    if (m_AffineTransformReadFlags[i] != 0)
      continue;

    std::string name = mu::get_name(m_ImageFileNames[i].c_str());

    std::string fn =
      firststr + name + suffixstr + std::string(".affine");

    muLogMacro(<< "Writing " << fn << "...\n");

    PairRegistrationMethod<InternalImagePixelType>::
      WriteAffineTransform(fn.c_str(), m_AffineTransforms[i]);
  }

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::ReadParameters()
{

  itkDebugMacro(<< "Read parameters");

  m_DoneRegistration = false;

  std::string firststr =
    m_OutputDirectory + mu::get_name(m_ImageFileNames[0].c_str()) +
    std::string("_to_");

  std::string suffixstr;
  if (m_Suffix.length() != 0)
    suffixstr = std::string("_") + m_Suffix;
  else
    suffixstr = std::string("");

  // Read template to image transforms
  if (m_TemplateFileName.length() != 0)
  {

    std::string name = mu::get_name(m_TemplateFileName.c_str());

    std::string fn =
      firststr + name + suffixstr + std::string(".affine");

    muLogMacro(<< "Reading " << fn << "...\n");

    try
    {
      m_TemplateAffineTransform =
        PairRegistrationMethod<InternalImagePixelType>::
          ReadAffineTransform(fn.c_str());
      m_AffineTransformReadFlags[0] = 1;
    }
    catch (...)
    {
      m_AffineTransformReadFlags[0] = 0;
    }

    if (m_WarpAtlas)
    {
      std::string warpfn =
        firststr + name + suffixstr + std::string(".bspline");

      muLogMacro(<< "Reading " << warpfn << "...\n");

      try   
      {     
        m_TemplateBSplineTransform =
          PairRegistrationMethod<InternalImagePixelType>::
            ReadBSplineTransform(warpfn.c_str());
        BSplineTransformType::SizeType splineGridSize =
          m_TemplateBSplineTransform->GetGridRegion().GetSize();
        if (splineGridSize[0] != (m_WarpGridX+3))
          itkExceptionMacro(<< "Grid size in file doesn't match input");
        if (splineGridSize[1] != (m_WarpGridY+3))
          itkExceptionMacro(<< "Grid size in file doesn't match input");
        if (splineGridSize[2] != (m_WarpGridZ+3))
          itkExceptionMacro(<< "Grid size in file doesn't match input");
        muLogMacro(<< "  Deformation obtained from file: " << warpfn << "\n");
        m_TemplateBSplineReadFlag = true;
      }
      catch(...)
      {
        m_TemplateBSplineReadFlag = false;
      }

    } // if do warping

  }

  // Read image to image transforms
  for (unsigned int i = 1; i < m_ImageFileNames.GetSize(); i++)
  {
    std::string name = mu::get_name(m_ImageFileNames[i].c_str());

    std::string fn =
      firststr + name + suffixstr + std::string(".affine");

    muLogMacro(<< "Reading " << fn << "...\n");

    try
    {
      m_AffineTransforms[i]  =
        PairRegistrationMethod<InternalImagePixelType>::
          ReadAffineTransform(fn.c_str());
      m_AffineTransformReadFlags[i] = 1;
    }
    catch (...)
    {
      m_AffineTransformReadFlags[i] = 0;
    }
  }

  bool allReadOK = true;

  for (unsigned i = 0; i < m_ImageFileNames.GetSize(); i++)
  {
    if (m_AffineTransformReadFlags[i] == 0)
    {
       allReadOK = false;
    }
  }
  // Can assume that registration has been done?
  if (allReadOK)
  {
    m_DoneRegistration = true;
  }

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetSuffix(std::string suffix)
{

  m_Suffix = suffix;

  m_Modified = true;

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetTemplateFileName(std::string filename)
{

  m_TemplateFileName = filename;

  m_TemplateAffineTransform = AffineTransformType::New();

  m_AffineTransformReadFlags[0] = 0;

  m_DoneRegistration = false;
  m_DoneResample = false;

  m_Modified = true;

}


template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetImageFileNames(StringList names)
{

  itkDebugMacro(<< "SetImageFileNames");

  unsigned int numImages = names.GetSize();

  if (numImages == 0)
    itkExceptionMacro(<< "No images specified");
  
  m_ImageFileNames = names;

  m_TemplateAffineTransform = AffineTransformType::New();

  // Clear previous transforms
  m_AffineTransforms.Clear();
  m_AffineTransforms.Allocate(numImages);
  for (unsigned int i = 0; i < numImages; i++)
  {
    // Also append identity matrix for each image
    AffineTransformPointer transform = AffineTransformType::New();
    m_AffineTransforms.Append(transform);
  }

  m_AffineTransformReadFlags = FlagArrayType(numImages);
  for (unsigned int i = 0; i < numImages; i++)
    m_AffineTransformReadFlags[i] = 0;

  m_DoneRegistration = false;
  m_DoneResample = false;

  m_Modified = true;

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetProbabilityFileNames(StringList names)
{

  unsigned int numProbabilities = names.GetSize();
  if (numProbabilities == 0)
  {
    itkExceptionMacro(<< "No probability images");
  }

  m_ProbabilityFileNames = names;

  m_DoneResample = false;

  m_Modified = true;

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetAtlasOrientation(std::string orient)
{

  m_AtlasOrientation = orient;

  m_DoneRegistration = false;
  m_DoneResample = false;

  m_Modified = true;

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetImageOrientations(StringList orientations)
{

  m_ImageOrientations = orientations;

  m_DoneRegistration = false;
  m_DoneResample = false;

  m_Modified = true;
}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::SetWarpGridSize(unsigned int nx, unsigned int ny, unsigned int nz)
{
  m_WarpGridX = nx;
  m_WarpGridY = ny;
  m_WarpGridZ = nz;

  m_Modified = true;
}


template <class TOutputPixel, class TProbabilityPixel>
typename AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::ProbabilityImageList
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::GetProbabilities()
{

  itkDebugMacro(<< "GetProbabilities");

  this->Update();

  return m_Probabilities;

}

template <class TOutputPixel, class TProbabilityPixel>
typename AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::OutputImageList
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::GetImages()
{

  itkDebugMacro(<< "GetImages");

  this->Update();

  return m_Images;

}

template <class TOutputPixel, class TProbabilityPixel>
typename AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::OutputImagePointer
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::GetTemplate()
{

  itkDebugMacro(<< "GetTemplate");

  this->Update();

  return m_Template;

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::Update()
{
  itkDebugMacro(<< "Update");

  if (m_Modified || !m_DoneRegistration)
    this->RegisterImages();

  if (m_Modified || !m_DoneResample)
    this->ResampleImages();

  if (m_Modified && m_WarpAtlas)
    this->DoBSplineWarp();

  m_Modified = false;

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::RegisterImages()
{

  itkDebugMacro(<< "RegisterImages");

  this->VerifyInitialization();

  typedef itk::ImageFileReader<InternalImageType> ReaderType;
  typedef typename ReaderType::Pointer ReaderPointer;

  typedef OrientImageModifier<InternalImageType> OrienterType;
  typedef typename OrienterType::Pointer OrienterPointer;

  // Define the blurring filter type
  typedef itk::DiscreteGaussianImageFilter<InternalImageType, InternalImageType>
    GaussianFilterType;
  typedef typename GaussianFilterType::Pointer GaussianFilterPointer;

  // Get the first image (for reference)
  InternalImagePointer first;
  {
    ReaderPointer reader = ReaderType::New();
    reader->SetFileName(m_ImageFileNames[0].c_str());

    reader->Update();

#if DO_BLUR
    itkDebugMacro(<< "Blurring first image with stddev = 0.1...");
    GaussianFilterPointer blurfilt = GaussianFilterType::New();
    blurfilt->SetMaximumError(1e-2);
    blurfilt->SetMaximumKernelWidth(32);
    blurfilt->SetVariance(0.01);
    blurfilt->SetInput(reader->GetOutput());
    blurfilt->Update();

    first = blurfilt->GetOutput();
#else
    first = reader->GetOutput();
#endif
  }

  // Register template to first image
  if ((m_TemplateFileName.length() != 0) && (m_AffineTransformReadFlags[0] == 0))
  {
    itkDebugMacro(<< "Registering template " << m_TemplateFileName << "...");
    ReaderPointer reader = ReaderType::New();
    reader->SetFileName(m_TemplateFileName.c_str());

    reader->Update();

    InternalImagePointer templateImg = reader->GetOutput();

    if (m_AtlasOrientation.length() != 0)
    {
      OrienterPointer orienter = OrienterType::New();
      orienter->SetSourceOrientation(m_AtlasOrientation);
      orienter->SetTargetOrientation(m_ImageOrientations[0]);
      orienter->Modify(templateImg);
    }

    muLogMacro(<< "Registering template to first image...\n");

    m_TemplateAffineTransform =
      PairRegistrationMethod<InternalImagePixelType>::
        RegisterAffine(first, templateImg);

  }

  // Register each image to first image
  for (unsigned int i = 1; i < m_ImageFileNames.GetSize(); i++)
  {
    muLogMacro(<< "Registering image " << i+1 << " to first image...\n");

    if (m_AffineTransformReadFlags[i] != 0)
      continue;

    ReaderPointer imgreader = ReaderType::New();
    imgreader->SetFileName(m_ImageFileNames[i].c_str());

    imgreader->Update();

    InternalImagePointer img_i = imgreader->GetOutput();

    OrienterPointer orienter = OrienterType::New();
    orienter->SetSourceOrientation(m_ImageOrientations[i]);
    orienter->SetTargetOrientation(m_ImageOrientations[0]);
    orienter->Modify(img_i);

#if DO_BLUR 
    itkDebugMacro(<< "Blurring image " << i+1 << " with stddev = 0.1...");
    GaussianFilterPointer blurfilt = GaussianFilterType::New();
    blurfilt->SetMaximumError(1e-2);
    blurfilt->SetMaximumKernelWidth(32);
    blurfilt->SetVariance(0.01);
    blurfilt->SetInput(img_i);
    blurfilt->Update();

    m_AffineTransforms[i] = 
      PairRegistrationMethod<InternalImagePixelType>::
        RegisterAffine(first, blurfilt->GetOutput());
#else
    m_AffineTransforms[i] = 
      PairRegistrationMethod<InternalImagePixelType>::
        RegisterAffine(first, img_i);
#endif

  }

  m_DoneRegistration = true;

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::ResampleImages()
{

  itkDebugMacro(<< "ResampleImages");

  if (!m_DoneRegistration)
    return;

  // Define the internal reader type
  typedef itk::ImageFileReader<InternalImageType> ReaderType;
  typedef typename ReaderType::Pointer ReaderPointer;

  // Orientation modifier
  typedef OrientImageModifier<InternalImageType> OrienterType;
  typedef typename OrienterType::Pointer OrienterPointer;

  // Get the first image (for reference)
  InternalImagePointer first;
  {
    ReaderPointer reader = ReaderType::New();
    reader->SetFileName(m_ImageFileNames[0].c_str());

    reader->Update();

    first = reader->GetOutput();
  }

  typedef itk::ResampleImageFilter<InternalImageType, InternalImageType>
    ResampleType;
  typedef typename ResampleType::Pointer ResamplePointer;

  typedef itk::LinearInterpolateImageFunction<InternalImageType, double>
    LinearInterpolatorType;
  typedef itk::BSplineInterpolateImageFunction<InternalImageType, double, double>
    SplineInterpolatorType;

  typename LinearInterpolatorType::Pointer linearInt =
    LinearInterpolatorType::New();

  // Spline interpolation, only available for input images, not atlas
  typename SplineInterpolatorType::Pointer splineInt =
    SplineInterpolatorType::New();
  splineInt->SetSplineOrder(3);

  // Resample the template
  if (m_TemplateFileName.length() != 0)
  {
    itkDebugMacro(<< "Resample template");

    ReaderPointer reader = ReaderType::New();
    reader->SetFileName(m_TemplateFileName.c_str());
    reader->Update();

    if (m_AtlasOrientation.length() != 0)
    {
      OrienterPointer orienter = OrienterType::New();
      orienter->SetSourceOrientation(m_AtlasOrientation);
      orienter->SetTargetOrientation(m_ImageOrientations[0]);
      orienter->Modify(reader->GetOutput());
    }

    ResamplePointer resampler = ResampleType::New();

    resampler->SetInput(reader->GetOutput());
    resampler->SetTransform(m_TemplateAffineTransform);

    resampler->SetInterpolator(linearInt);
    resampler->SetSize(first->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputOrigin(first->GetOrigin());
    resampler->SetOutputSpacing(first->GetSpacing());
    resampler->SetDefaultPixelValue(0);

    resampler->Update();

    m_Template = CopyOutputImage(resampler->GetOutput());
  }

  // Resample the probabilities
  for (unsigned int i = 0; i < m_Probabilities.GetSize(); i++)
    m_Probabilities[i] = 0;
  m_Probabilities.Clear();
  for (unsigned int i = 0; i < m_ProbabilityFileNames.GetSize(); i++)
  {
    itkDebugMacro(<< "Resample probability " << i);

    ReaderPointer reader = ReaderType::New();
    reader->SetFileName(m_ProbabilityFileNames[i].c_str());

    reader->Update();

    if (m_AtlasOrientation.length() != 0)
    {
      OrienterPointer orienter = OrienterType::New();
      orienter->SetSourceOrientation(m_AtlasOrientation);
      orienter->SetTargetOrientation(m_ImageOrientations[0]);
      orienter->Modify(reader->GetOutput());
    }

    ResamplePointer resampler = ResampleType::New();

    resampler->SetInput(reader->GetOutput());
    resampler->SetTransform(m_TemplateAffineTransform);

    resampler->SetInterpolator(linearInt);
    resampler->SetSize(first->GetLargestPossibleRegion().GetSize());
    resampler->SetOutputOrigin(first->GetOrigin());
    resampler->SetOutputSpacing(first->GetSpacing());
    resampler->SetDefaultPixelValue(0);

    resampler->Update();

    m_Probabilities.Append(CopyProbabilityImage(resampler->GetOutput()));
  }

  // Clear image list
  for (unsigned int i = 0; i < m_Images.GetSize(); i++)
    m_Images[i] = 0;
  m_Images.Clear();

  // Do nothing for first image
  m_Images.Append(CopyOutputImage(first));

  // The FOV mask, regions where intensities in all channels do not
  // match FOV code
  m_FOVMask = ByteImageType::New();
  m_FOVMask->SetRegions(m_Images[0]->GetLargestPossibleRegion());
  m_FOVMask->Allocate();
  m_FOVMask->SetOrigin(m_Images[0]->GetOrigin());
  m_FOVMask->SetSpacing(m_Images[0]->GetSpacing());

  typedef itk::ImageRegionIterator<ByteImageType> MaskIteratorType;
  MaskIteratorType maskIt(m_FOVMask, m_FOVMask->GetLargestPossibleRegion());

  maskIt.GoToBegin();
  while (!maskIt.IsAtEnd())
  {
    maskIt.Set(1);
    ++maskIt;
  }
   
  // Resample the other images
  for (unsigned int i = 1; i < m_ImageFileNames.GetSize(); i++)
  {
    ReaderPointer reader = ReaderType::New();
    reader->SetFileName(m_ImageFileNames[i].c_str());

    reader->Update();

    OrienterPointer orienter = OrienterType::New();
    orienter->SetSourceOrientation(m_ImageOrientations[i]);
    orienter->SetTargetOrientation(m_ImageOrientations[0]);
    orienter->Modify(reader->GetOutput());

    ResamplePointer resampler = ResampleType::New();

    resampler->SetInput(reader->GetOutput());
    resampler->SetTransform(m_AffineTransforms[i]);

    if (m_UseNonLinearInterpolation)
      resampler->SetInterpolator(splineInt);
    else
      resampler->SetInterpolator(linearInt);

    resampler->SetDefaultPixelValue(m_OutsideFOVCode);
    resampler->SetOutputOrigin(first->GetOrigin());
    resampler->SetOutputSpacing(first->GetSpacing());
    resampler->SetSize(first->GetLargestPossibleRegion().GetSize());

    resampler->Update();

    InternalImagePointer tmp = resampler->GetOutput();

    // Zero the mask region outside FOV and also the intensities with outside
    // FOV code
    typedef itk::ImageRegionIterator<InternalImageType> InternalIteratorType;

    InternalIteratorType tmpIt(tmp, first->GetLargestPossibleRegion());

    maskIt.GoToBegin();
    tmpIt.GoToBegin();
    while (!maskIt.IsAtEnd())
    {
      if (tmpIt.Get() == m_OutsideFOVCode)
      {
        maskIt.Set(0);
        tmpIt.Set(0);
      }
      ++maskIt;
      ++tmpIt;
    }

    // Add the image
    m_Images.Append(CopyOutputImage(tmp));
  }

  m_DoneResample = true;

}

template <class TOutputPixel, class TProbabilityPixel>
void
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::DoBSplineWarp()
{
  itkDebugMacro(<< "DoBSplineWarp");

  if ((m_TemplateFileName.length() == 0) || !m_WarpAtlas)
    return;

  typedef itk::ResampleImageFilter<InternalImageType, InternalImageType>
    ResampleType;
  typedef typename ResampleType::Pointer ResamplePointer;

  typedef itk::LinearInterpolateImageFunction<InternalImageType, double>
    LinearInterpolatorType;

  typename LinearInterpolatorType::Pointer linearInt =
    LinearInterpolatorType::New();

  // Compute B-spline transform if not read from file
  if (!m_TemplateBSplineReadFlag)
  {
    itkDebugMacro(<< "Warping template " << m_TemplateFileName << "...");

    muLogMacro(<< "Warping template to first image...\n");
    muLogMacro(<< "  Computing B-spline coefficients...\n");

    // TODO: crop images

    m_TemplateBSplineTransform =
      PairRegistrationMethod<InternalImagePixelType>::
        RegisterBSpline(m_Images[0], m_Template,
        m_WarpGridX, m_WarpGridY, m_WarpGridZ);

    // TODO: undo crop

  }

  // Warp template to first image
  {
    ResamplePointer warper = ResampleType::New();

    warper->SetInput(m_Template);
    warper->SetTransform(m_TemplateBSplineTransform);
    
    warper->SetInterpolator(linearInt);
    warper->SetSize(m_Images[0]->GetLargestPossibleRegion().GetSize());
    warper->SetOutputOrigin(m_Images[0]->GetOrigin());
    warper->SetOutputSpacing(m_Images[0]->GetSpacing());
    warper->SetDefaultPixelValue(0);

    warper->Update();

    m_Template = CopyOutputImage(warper->GetOutput());
  }

  // Warp probabilities
  DynArray<ProbabilityImagePointer> oldProbs = m_Probabilities;

  m_Probabilities.Clear();

  for (unsigned int i = 0; i < oldProbs.GetSize(); i++)
  {
    ResamplePointer warper = ResampleType::New();
  
    warper->SetInput(oldProbs[i]);
    warper->SetTransform(m_TemplateBSplineTransform);
    
    warper->SetInterpolator(linearInt);
    warper->SetSize(m_Images[0]->GetLargestPossibleRegion().GetSize());
    warper->SetOutputOrigin(m_Images[0]->GetOrigin());
    warper->SetOutputSpacing(m_Images[0]->GetSpacing());
    warper->SetDefaultPixelValue(0);

    warper->Update();

    m_Probabilities.Append(CopyProbabilityImage(warper->GetOutput()));
  }

}

template <class TOutputPixel, class TProbabilityPixel>
typename AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::OutputImagePointer
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::CopyOutputImage(
  typename AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
  ::InternalImagePointer img)
{

  itkDebugMacro(<< "CopyOutputImage");

  OutputImagePointer outimg = OutputImageType::New();

  outimg->SetRegions(img->GetLargestPossibleRegion());
  outimg->Allocate();

  outimg->SetOrigin(img->GetOrigin());
  outimg->SetSpacing(img->GetSpacing());

  typedef itk::ImageRegionIterator<InternalImageType> InternalIteratorType;
  InternalIteratorType inputIter(img, img->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<OutputImageType> OutputIteratorType;
  OutputIteratorType outputIter(outimg, outimg->GetLargestPossibleRegion());

  inputIter.GoToBegin();
  outputIter.GoToBegin();

  while (!inputIter.IsAtEnd())
  {
    outputIter.Set(static_cast<OutputImagePixelType>(inputIter.Get()));
    ++inputIter;
    ++outputIter;
  }

  return outimg;

}

template <class TOutputPixel, class TProbabilityPixel>
typename AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::ProbabilityImagePointer
AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
::CopyProbabilityImage(
  typename AtlasRegistrationMethod<TOutputPixel, TProbabilityPixel>
  ::InternalImagePointer img)
{

  itkDebugMacro(<< "CopyProbabilityImage");

  ProbabilityImagePointer outimg = ProbabilityImageType::New();

  outimg->SetRegions(img->GetLargestPossibleRegion());
  outimg->Allocate();

  outimg->SetOrigin(img->GetOrigin());
  outimg->SetSpacing(img->GetSpacing());

  typedef itk::ImageRegionIterator<InternalImageType> InternalIteratorType;
  InternalIteratorType inputIter(img, img->GetLargestPossibleRegion());

  typedef itk::ImageRegionIterator<ProbabilityImageType>
    ProbabilityIteratorType;
  ProbabilityIteratorType outputIter(outimg,
    outimg->GetLargestPossibleRegion());

  inputIter.GoToBegin();
  outputIter.GoToBegin();

  while (!inputIter.IsAtEnd())
  {
    double p = inputIter.Get();
    if (p < 0.0)
      p = 0.0;
    outputIter.Set(static_cast<ProbabilityImagePixelType>(p));
    ++inputIter;
    ++outputIter;
  }

  return outimg;

}


#endif
