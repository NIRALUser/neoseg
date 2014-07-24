
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include "MRAnisotropicFilter.h"

#include <iostream>
#include <sstream>
#include <string>

#include <math.h>

int
main(int argc, char** argv)
{

  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <image_1> ... <image_n>"
      << std::endl;
    return -1;
  }

  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  typedef itk::Image<unsigned short, 3> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  DynArray<ImageType::Pointer> images;

  try
  {

    for (unsigned int i = 1; i < argc; i++)
    {
      std::cout << "Reading " << argv[i] << std::endl;
      ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName(argv[i]);
      reader->Update();

      images.Append(reader->GetOutput());
    }

  }
  catch (itk::ExceptionObject& exp)
  {
    std::cerr << exp << std::endl;
    return -1;
  }

  typedef MRAnisotropicFilter<ImageType> FilterType;
  FilterType::Pointer filt = FilterType::New();

  filt->DebugOn();

  filt->SetNumberOfIterations(10);
  filt->SetTimeStep(0.01);

  filt->Filter(images);

  // Write output images
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  for (unsigned int i = 0; i < images.GetSize(); i++)
  {

    std::stringstream oss;
    oss << "anisofilt_" << i << ".gipl" << std::ends;
    writer->SetFileName(oss.str().c_str());
    writer->SetInput(images[i]);
    writer->Update();
  }

}
