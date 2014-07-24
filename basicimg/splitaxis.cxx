
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include "SplitAxisImageFilter.h"

#include <iostream>

#include <stdlib.h>

int
main(int argc, char** argv)
{

  if (argc != 5)
  {
    std::cerr << "Usage: " << argv[0]
      << " <axis> <imgfn> <out_evenfn> <out_oddfn>" << std::endl;
    return -1;
  }

  int splitAxis = atoi(argv[1]);
  char* imgfn = argv[2];

  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  typedef short PixelType;
  typedef itk::Image<PixelType, 3> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;

  ImageType::Pointer img = 0;
  {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(imgfn);
    reader->Update();

    img = reader->GetOutput();
  }

  typedef SplitAxisImageFilter<ImageType> SplitterType;
  SplitterType::Pointer splitter = SplitterType::New();

  splitter->SetSplitAxis(splitAxis);
  splitter->SetInput(img);

  std::cout << "Splitting image along axis " << splitAxis << "..." << std::endl;
  splitter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();


  std::cout << "Writing even image..." << std::endl;
  writer->SetFileName(argv[3]);
  writer->SetInput(splitter->GetEvenImage());
  writer->Update();

  std::cout << "Writing odd image..." << std::endl;
  writer->SetFileName(argv[4]);
  writer->SetInput(splitter->GetOddImage());
  writer->Update();

  return 0;

}
