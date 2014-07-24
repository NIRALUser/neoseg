
#include "mu.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include "ConnectedComponentsFilter.h"

#include <iostream>

int
main(int argc, char** argv)
{

  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " <imgfn>" << std::endl;
    return -1;
  }

  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  typedef itk::Image<unsigned char, 3> InputImageType;
  typedef itk::Image<unsigned short, 3> OutputImageType;

  typedef itk::ImageFileReader<InputImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(argv[1]);

  reader->Update();
  reader->GetImageIO()->SetByteOrderToBigEndian();
  reader->Update();

  typedef itk::BinaryBallStructuringElement<unsigned char, 3> StructElementType;
  typedef
    itk::BinaryDilateImageFilter<InputImageType, InputImageType,
      StructElementType> DilateType;
  typedef
    itk::BinaryErodeImageFilter<InputImageType, InputImageType,
      StructElementType> ErodeType;

  StructElementType structel;

  StructElementType::RadiusType radius;
  for (unsigned int i = 0; i < 3; i++)
    radius[i] = 1;
  structel.SetRadius(radius);

  ErodeType::Pointer erode = ErodeType::New();
  erode->SetKernel(structel);
  erode->SetInput(reader->GetOutput());

  erode->Update();

  DilateType::Pointer dil = DilateType::New();
  dil->SetKernel(structel);
  dil->SetInput(erode->GetOutput());

  dil->Update();

  typedef ConnectedComponentsFilter<InputImageType, OutputImageType> CCType;
  CCType::Pointer cc = CCType::New();

  cc->DebugOn();

  cc->SetInput(dil->GetOutput());
  cc->Update();

  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName("ccout.gipl");
  writer->SetInput(cc->GetOutput());
  writer->Update();

  return 0;

}

