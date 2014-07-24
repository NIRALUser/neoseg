
#include "itkAffineTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include <iostream>
#include <fstream>

#include "PairRegistrationMethod.h"

int
main(int argc, char **argv)
{

  if (argc != 5)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << "<affine> <target> <label> <out>" << std::endl;
    return 1;
  }

  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  std::cout << "Transform: " << argv[1] << std::endl;
  std::cout << "Target image: " << argv[2] << std::endl;
  std::cout << "Label image: " << argv[3] << std::endl;

  const unsigned int Dimension = 3;
  typedef unsigned char PixelType;

  typedef itk::Image<PixelType, Dimension>  LabelImageType;

  typedef PairRegistrationMethod<PixelType>::AffineTransformType TransformType;

  typedef itk::ImageFileReader<LabelImageType> ImageReaderType;

  // Read the affine transform
  TransformType::Pointer transform =
    PairRegistrationMethod<PixelType>::ReadAffineTransform(argv[1]);

  // Read the label image
  ImageReaderType::Pointer  reader1 = ImageReaderType::New();
  ImageReaderType::Pointer  reader2 = ImageReaderType::New();

  LabelImageType::Pointer targetImg;
  LabelImageType::Pointer labelImg;

  std::cout << "Reading images..." << std::endl;
  try
  {
    reader1->SetFileName(argv[2]);
    reader1->Update();

    targetImg = reader1->GetOutput();

    reader2->SetFileName(argv[3]);
    reader2->Update();

    labelImg = reader2->GetOutput();
  }
  catch (itk::ExceptionObject& exc)
  {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << exc << std::endl;
  }

  // Resample label image
  std::cout << "Resampling..." << std::endl;

  typedef itk::ResampleImageFilter<LabelImageType, LabelImageType>
    ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  resampler->SetTransform(transform);
  resampler->SetInput(labelImg);

  typedef itk::NearestNeighborInterpolateImageFunction<LabelImageType, double>
    InterpolatorType;

  InterpolatorType::Pointer interp = InterpolatorType::New();

  resampler->SetInterpolator(interp);
  resampler->SetSize(targetImg->GetLargestPossibleRegion().GetSize());
  resampler->SetOutputOrigin(targetImg->GetOrigin());
  resampler->SetOutputSpacing(targetImg->GetSpacing());
  resampler->SetDefaultPixelValue(0);

  typedef itk::ImageFileWriter<LabelImageType>  WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName(argv[4]);
  writer->SetInput(resampler->GetOutput() );
  writer->Update();

  return 0;

}

