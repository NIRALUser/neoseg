
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include "AtlasCropImageSource.h"

#include <iostream>

int
main(int argc, char** argv)
{

  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <maskfn> <imgfn>" << std::endl;
    return -1;
  }

  char* labelfn = argv[1];
  char* imgfn = argv[2];

  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  typedef itk::Image<unsigned char, 3> ByteImageType;

  typedef itk::ImageFileReader<ByteImageType> ReaderType;

  ByteImageType::Pointer labelImg = 0;
  {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(labelfn);
    reader->Update();

    labelImg = reader->GetOutput();
  }

  typedef AtlasCropImageSource<ByteImageType, ByteImageType> CropType;
  CropType::Pointer cropper = CropType::New();

  cropper->SetPadding(5.0);

  DynArray<ByteImageType::Pointer> plist;
  plist.Append(labelImg);

  cropper->UseProbabilities(plist);

  ByteImageType::Pointer img = 0;
  {
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(imgfn);
    reader->Update();

    img = reader->GetOutput();
  }

  typedef itk::ImageFileWriter<ByteImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  std::string fn;

/*
  fn = "crop_";
  fn += labelfn;

  writer->SetFileName(fn.c_str());
  writer->SetInput(cropper->Crop(labelImg));
  writer->Update();
*/

  fn = "crop_";
  fn += imgfn;

  writer->SetFileName(fn.c_str());
  writer->SetInput(cropper->Crop(img));
  writer->Update();

  fn = "restore_";
  fn += imgfn;

  writer->SetFileName(fn.c_str());
  writer->SetInput(cropper->Restore(cropper->Crop(img)));
  writer->Update();

/*
  ByteImageType::Pointer res = cropper->Restore(cropper->Crop(img));

  ByteImageType::IndexType ind;
  ByteImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();

  for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
    for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
      for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
        if (img->GetPixel(ind) != res->GetPixel(ind))
        {
          std::cerr << "Mismatch" << std::endl;
        }
*/

  return 0;

}
