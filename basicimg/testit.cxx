
#include "itkImage.h"
#include "itkImageFileWriter.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include "SubsampledImageRegionIterator.h"

#include <iostream>

int
main(int argc, char** argv)
{

  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  typedef itk::Image<unsigned char, 3> ByteImageType;

  ByteImageType::SizeType size;
  size[0] = 64;
  size[1] = 64;
  size[2] = 32;

  ByteImageType::SpacingType spacing;
  spacing[0] = 1;
  spacing[1] = 1;
  spacing[2] = 2;

  ByteImageType::RegionType region;
  region.SetSize(size);

  ByteImageType::Pointer img = ByteImageType::New();
  img->SetRegions(region);
  img->Allocate();
  img->SetSpacing(spacing);

  ByteImageType::IndexType ind;
  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        img->SetPixel(ind, 0);

  typedef SubsampledImageRegionIterator<ByteImageType> IteratorType;
  IteratorType it(img, region, 4.0);

  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    ByteImageType::IndexType ind = it.GetIndex();

    if (ind[0] < 32)
      it.Set(127);
    else
      it.Set(255);

std::cout << "Iterator ind = " << ind << std::endl;
  }

  typedef itk::ImageFileWriter<ByteImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFileName("testit.mha");
  writer->SetInput(img);
  writer->Update();

  return 0;

}
