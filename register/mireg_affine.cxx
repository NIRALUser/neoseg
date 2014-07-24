
#include "itkAffineTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMinimumMaximumImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include <iostream>
#include <fstream>

#include "Timer.h"
#include "muFile.h"

#include "PairRegistrationMethod.h"

int
main(int argc, char **argv)
{

  if (argc < 3)
  {
    std::cerr << "Affine MI registation for a pair of images" << std::endl;
    std::cerr << "Registers fixed image and first moving image" << std::endl;
    std::cerr << "Assumes all images have the same orientation" << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixed  moving [moving_2 ... moving_n]" << std::endl;
    return 1;
  }

  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  std::cout << "Fixed image: " << argv[1] << std::endl;
  std::cout << "Moving image: " << argv[2] << std::endl;

  const unsigned int Dimension = 3;
  typedef float PixelType;

  typedef itk::Image<PixelType, Dimension>  FixedImageType;
  typedef itk::Image<PixelType, Dimension>  MovingImageType;

  typedef PairRegistrationMethod<PixelType>::AffineTransformType TransformType;

  typedef itk::ImageFileReader<FixedImageType> FixedImageReaderType;
  typedef itk::ImageFileReader<MovingImageType> MovingImageReaderType;

  // Read the images
  FixedImageReaderType::Pointer  fixedImageReader =
    FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader =
    MovingImageReaderType::New();

  fixedImageReader->SetFileName(argv[1]);
  movingImageReader->SetFileName(argv[2]);

  std::cout << "Reading images..." << std::endl;
  try
  {
    fixedImageReader->Update();
    movingImageReader->Update();
  }
  catch (itk::ExceptionObject& exc)
  {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << exc << std::endl;
  }

/*
  std::cout << "Blurring first moving image..." << std::endl;

  typedef itk::DiscreteGaussianImageFilter<MovingImageType, MovingImageType>
    BlurFilterType;

  BlurFilterType::Pointer blur = BlurFilterType::New();
  blur->SetInput(movingImageReader->GetOutput());
  blur->SetVariance(0.25);
  blur->Update();
*/

  Timer* timer = new Timer();

  TransformType::Pointer finalTransform =
    PairRegistrationMethod<PixelType>::
      RegisterAffine(fixedImageReader->GetOutput(),
        movingImageReader->GetOutput());
//        blur->GetOutput());

  timer->Stop();

  std::cout << "Registration took " << timer->GetElapsedHours() << " hours, ";
  std::cout << timer->GetElapsedMinutes() << " minutes, ";
  std::cout << timer->GetElapsedSeconds() << " seconds";
  std::cout << std::endl;

  //
  // Write parameters to file
  //
  std::cout << "Writing parameters..." << std::endl;
  std::string paramfn;
  {
    std::ostringstream oss;
    oss << mu::get_name(argv[1]) << "_to_" << mu::get_name(argv[2])
      << ".affine" << std::ends;
    paramfn = oss.str();
  }
  PairRegistrationMethod<PixelType>::
    WriteAffineTransform(paramfn.c_str(), finalTransform);

  // Resample moving images
  for (unsigned int k = 2; k < argc; k++)
  {
    std::cout << "Resampling " << argv[k] << "..." << std::endl;

    // Read image
    MovingImageReaderType::Pointer reader = MovingImageReaderType::New();
    reader->SetFileName(argv[k]);
    reader->Update();

    MovingImageType::Pointer tmp = reader->GetOutput();

    typedef itk::ResampleImageFilter<MovingImageType, FixedImageType>
      ResampleFilterType;
    ResampleFilterType::Pointer resample = ResampleFilterType::New();

    resample->SetTransform( finalTransform );
    resample->SetInput(tmp);

    //typedef itk::BSplineInterpolateImageFunction<MovingImageType, double, double>
    //  InterpolatorType;
    typedef itk::LinearInterpolateImageFunction<MovingImageType, double>
      InterpolatorType;
    //typedef itk::NearestNeighborInterpolateImageFunction<MovingImageType, double>
    //  InterpolatorType;

    InterpolatorType::Pointer interp = InterpolatorType::New();
    //interp->SetSplineOrder(5);

    resample->SetInterpolator(interp);

    FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

    resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
    resample->SetOutputOrigin(fixedImage->GetOrigin());
    resample->SetOutputSpacing(fixedImage->GetSpacing());
    resample->SetDefaultPixelValue(0);

    typedef unsigned short OutputPixelType;
    typedef itk::Image<OutputPixelType, Dimension> OutputImageType;

    typedef itk::RescaleIntensityImageFilter<FixedImageType, OutputImageType >
      RescalerType;
    RescalerType::Pointer rescaler = RescalerType::New();

    // Rescale to original range
    typedef itk::MinimumMaximumImageFilter<MovingImageType> 
      MinMaxFilterType;
    MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
    minMaxFilter->SetInput(tmp);
    minMaxFilter->Update();

    rescaler->SetOutputMinimum((unsigned short)minMaxFilter->GetMinimum());
    rescaler->SetOutputMaximum((unsigned short)minMaxFilter->GetMaximum());
    rescaler->SetInput(resample->GetOutput());

    typedef itk::ImageFileWriter< OutputImageType >  WriterType;
    WriterType::Pointer writer = WriterType::New();

    std::ostringstream oss;
    oss << "mireg_affine_" << mu::get_name(argv[k]) << ".mha" << std::ends;

    std::cout << "Writing " << oss.str() << "..." << std::endl;
    writer->SetFileName(oss.str().c_str());
    writer->SetInput(rescaler->GetOutput() );
    writer->Update();

  } // for k

  return 0;

}

