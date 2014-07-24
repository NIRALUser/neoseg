
#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNumericTraits.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkVersion.h"

#include "AtlasCropImageSource.h"

#include "EMSParameters.h"
#include "EMSParametersXMLFile.h"

#include "DynArray.h"
#include "Log.h"
#include "Timer.h"

#include "muFile.h"

// Use manually instantiated classes for the big program chunks
#define MU_MANUAL_INSTANTIATION
#include "EMSegmentationFilter.h"
#include "AtlasRegistrationMethod.h"
#include "PairRegistrationMethod.h"
#undef MU_MANUAL_INSTANTIATION

#include "filterFloatImages.h"
#include "runEMS.h"

#include <iostream>
#include <string>
#include <sstream>

typedef itk::Image<float, 3> FloatImageType;
typedef itk::Image<unsigned char, 3> ByteImageType;
typedef itk::Image<short, 3> ShortImageType;

typedef FloatImageType::Pointer FloatImagePointer;
typedef ByteImageType::Pointer ByteImagePointer;
typedef ShortImageType::Pointer ShortImagePointer;

void
runEMS(EMSParameters* emsp, bool debugflag, bool writemoreflag)
{

  if (!emsp->CheckValues())
    throw std::string("Invalid segmentation parameter values");
  
  // Create and start a new timer (for the whole process)
  Timer* timer = new Timer();

  // Directory separator string
  std::string separator = std::string("/");
  separator[0] = MU_DIR_SEPARATOR;

  // Get output directory
  std::string outdir = emsp->GetOutputDirectory();
  // Make sure last character in output directory string is a separator
  if (outdir[outdir.size()-1] != MU_DIR_SEPARATOR)
    outdir += separator;

  // Create the output directory, stop if it does not exist
  if(!mu::create_dir(outdir.c_str()))
    return;

  // Set up the logger
  {
    std::string logfn = outdir + emsp->GetSuffix() + ".log";
    (Log::GetInstance())->EchoOn();
    (Log::GetInstance())->SetOutputFileName(logfn.c_str());
  }

  // Write out the parameters in XML
  {
    std::string xmlfn = outdir + emsp->GetSuffix() + ".xml";
    writeEMSParametersXML(xmlfn.c_str(), emsp);
  }

  // Set up suffix string for images
  std::string fmt = emsp->GetOutputFormat();
  std::string outext = ".mha";
  if (fmt.compare("Analyze") == 0)
    outext = ".hdr";
  if (fmt.compare("GIPL") == 0)
    outext = ".gipl";
  std::string suffstr =
    std::string("_") + std::string(emsp->GetSuffix()) + outext;
  std::string metasuffstr =
    std::string("_") + std::string(emsp->GetSuffix()) + std::string(".mha");

  muLogMacro(<< "mu::brainseg\n");
  muLogMacro(<< "========================================\n");
  muLogMacro(<< "Program compiled on: " << __DATE__ << "\n");
  muLogMacro(<< "\n");

  muLogMacro(<< "Using ITK version "
    << itk::Version::GetITKMajorVersion() << "."
    << itk::Version::GetITKMinorVersion() << "."
    << itk::Version::GetITKBuildVersion() <<  "\n");
  muLogMacro(<< "\n");

  // Write input parameters
  muLogMacro(<< "=== Parameters ===\n");
  muLogMacro(<< "Suffix: " << emsp->GetSuffix() << "\n");
  muLogMacro(<< "Atlas Directory: " << emsp->GetAtlasDirectory() << "\n");
  muLogMacro(<< "Atlas Orientation: " << emsp->GetAtlasOrientation() << "\n");
  muLogMacro(<< "Output Directory: " << emsp->GetOutputDirectory() << "\n");
  muLogMacro(<< "Output Format: " << emsp->GetOutputFormat() << "\n");
  muLogMacro(<< "Input images: \n");
  for (unsigned int i = 0; i < emsp->GetImages().GetSize(); i++)
    muLogMacro(<< "  " << "[" << (emsp->GetImageOrientations())[i] << "] " <<
      (emsp->GetImages())[i] << "\n");
  muLogMacro(
    << "Non-linear filtering, method: " << emsp->GetFilterMethod() << ", "
    << emsp->GetFilterIterations()
    << " iterations, dt = " << emsp->GetFilterTimeStep() << "\n");
  muLogMacro(
    << "Prior weight scales: " << emsp->GetPrior1() << ", "
    << emsp->GetPrior2() << ", " << emsp->GetPrior3()
    << ", " << emsp->GetPrior4() << "\n");
  muLogMacro(
    << "Max bias polynomial degree: " << emsp->GetMaxBiasDegree() << "\n");
  muLogMacro(<< "Atlas warping: " << emsp->GetDoAtlasWarp() << "\n");
  muLogMacro(
    << "Atlas warp spline grid size: " << emsp->GetAtlasWarpGridX() << " X "
    << emsp->GetAtlasWarpGridY() << " X "
    << emsp->GetAtlasWarpGridZ() << "\n");

  muLogMacro(<< "\n");

  muLogMacro(<< "=== Start ===\n");

  //Get atlas directory
  std::string atlasdir = emsp->GetAtlasDirectory();
  // Make sure last character is a separator
  if (atlasdir[atlasdir.size()-1] != MU_DIR_SEPARATOR)
    atlasdir += separator;

  muLogMacro(<< "Registering images using affine transform...\n");

  ByteImagePointer fovmask;

  FloatImagePointer croppedTemplate;

  DynArray<FloatImagePointer> images;
  DynArray<FloatImagePointer> priors;
  {
    typedef AtlasRegistrationMethod<float, float> AtlasRegType;
    AtlasRegType::Pointer atlasreg = AtlasRegType::New();

    if (debugflag)
      atlasreg->DebugOn();

    atlasreg->SetSuffix(emsp->GetSuffix());

    std::string templatefn = atlasdir + std::string("template.gipl");
    atlasreg->SetTemplateFileName(templatefn);

    atlasreg->SetAtlasOrientation(emsp->GetAtlasOrientation());

    atlasreg->SetImageFileNames(emsp->GetImages());
    atlasreg->SetImageOrientations(emsp->GetImageOrientations());
    atlasreg->SetOutputDirectory(outdir);

    // Compute list of file names for the priors
    DynArray<std::string> priorfnlist;
    {
      priorfnlist.Append(atlasdir + std::string("white.gipl"));
      priorfnlist.Append(atlasdir + std::string("gray.gipl"));
      priorfnlist.Append(atlasdir + std::string("csf.gipl"));
      priorfnlist.Append(atlasdir + std::string("rest.gipl"));
    }

    atlasreg->SetProbabilityFileNames(priorfnlist);

    atlasreg->SetWarpAtlas(emsp->GetDoAtlasWarp());
    atlasreg->SetWarpGridSize(
      emsp->GetAtlasWarpGridX(), 
      emsp->GetAtlasWarpGridY(), 
      emsp->GetAtlasWarpGridZ());

    muLogMacro(<< "Attempting to read previous registration results..."
      << std::endl);
    atlasreg->ReadParameters();

    muLogMacro(<< "Registering and resampling images..." << std::endl);
    Timer* regtimer = new Timer();
    atlasreg->Update();

    atlasreg->WriteParameters();
    regtimer->Stop();

    muLogMacro(<< "Registration took " << regtimer->GetElapsedHours() << " hours, ");
    muLogMacro(<< regtimer->GetElapsedMinutes() << " minutes, ");
    muLogMacro(<< regtimer->GetElapsedSeconds() << " seconds\n");
    delete regtimer;

    fovmask = atlasreg->GetFOVMask();

    images = atlasreg->GetImages();
    priors = atlasreg->GetProbabilities();

    croppedTemplate = atlasreg->GetTemplate();

    // Write the registered template and images
    if (writemoreflag)
    {
      muLogMacro(<< "Writing registered template...\n");

      typedef itk::RescaleIntensityImageFilter<FloatImageType, ByteImageType>
        ByteRescaleType;

      ByteRescaleType::Pointer rescaler = ByteRescaleType::New();
      rescaler->SetOutputMinimum(0);
      rescaler->SetOutputMaximum(255);
      rescaler->SetInput(atlasreg->GetTemplate());
      rescaler->Update();

      typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
      ByteWriterType::Pointer writer = ByteWriterType::New();

      std::string fn =
        outdir + mu::get_name((emsp->GetImages()[0]).c_str()) +
        std::string("_template") + suffstr;

      writer->SetInput(rescaler->GetOutput());
      writer->SetFileName(fn.c_str());
      writer->SetUseCompression(true);
      writer->Update();

      for (unsigned int i = 0; i < images.GetSize(); i++)
      {
        typedef itk::RescaleIntensityImageFilter<FloatImageType, ShortImageType>
          ShortRescaleType;

        ShortRescaleType::Pointer rescaler = ShortRescaleType::New();
        rescaler->SetOutputMinimum(0);
        rescaler->SetOutputMaximum(itk::NumericTraits<short>::max());
        rescaler->SetInput(images[i]);
        rescaler->Update();

        std::string fn =
          outdir + mu::get_name((emsp->GetImages()[i]).c_str()) +
          std::string("_registered") + suffstr;

        typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
        ShortWriterType::Pointer writer = ShortWriterType::New();

        writer->SetInput(rescaler->GetOutput());
        writer->SetFileName(fn.c_str());
        writer->SetUseCompression(true);
        writer->Update();
      }
    }
  }

  // Generate cropped images
  muLogMacro(<< "Cropping images based on atlas priors...\n");

  typedef AtlasCropImageSource<FloatImageType, FloatImageType> CropSourceType;
  CropSourceType::Pointer cropper = CropSourceType::New();

  typedef AtlasCropImageSource<ByteImageType, FloatImageType>
    ByteCropSourceType;
  ByteCropSourceType::Pointer byteCropper = ByteCropSourceType::New();

  typedef AtlasCropImageSource<ShortImageType, FloatImageType>
    ShortCropSourceType;
  ShortCropSourceType::Pointer shortCropper = ShortCropSourceType::New();

  DynArray<FloatImagePointer> croppedPriors;
  DynArray<FloatImagePointer> croppedImages;
  {
    DynArray<FloatImagePointer> brainprobs;
    for (unsigned int i = 0; i < (priors.GetSize()-1); i++)
      brainprobs.Append(priors[i]);

    cropper->UseProbabilities(brainprobs);
    byteCropper->UseProbabilities(brainprobs);
    shortCropper->UseProbabilities(brainprobs);

    for (unsigned int i = 0; i < images.GetSize(); i++)
      croppedImages.Append(cropper->Crop(images[i]));

    for (unsigned int i = 0; i < priors.GetSize(); i++)
      croppedPriors.Append(cropper->Crop(priors[i]));

    croppedTemplate = cropper->Crop(croppedTemplate);
  }

  // Crop FOV mask
  fovmask = byteCropper->Crop(fovmask);

  // Delete original images and atlas priors from memory, no longer needed
  images.Clear();
  priors.Clear();

  if (emsp->GetFilterIterations() > 0)
  {
    muLogMacro(<< "Non-linear filtering of registered images...\n");

    std::string method = emsp->GetFilterMethod();

    Timer* filtertimer = new Timer();

    filterFloatImages(croppedImages, method,
      emsp->GetFilterIterations(), emsp->GetFilterTimeStep());

    filtertimer->Stop();

    muLogMacro(<< "Non-linear filtering took " << filtertimer->GetElapsedHours() << " hours, ");
    muLogMacro(<< filtertimer->GetElapsedMinutes() << " minutes, ");
    muLogMacro(<< filtertimer->GetElapsedSeconds() << " seconds\n");

    delete filtertimer;
  }

  muLogMacro(<< "Rescale intensity of filtered cropped images...\n");
  {
    typedef itk::RescaleIntensityImageFilter<FloatImageType, FloatImageType>
      CropRescaleType;
    CropRescaleType::Pointer cropRescaler = CropRescaleType::New();
    cropRescaler->SetOutputMinimum(0);
    cropRescaler->SetOutputMaximum(4095);

    FloatImageType::SizeType size =
      croppedImages[0]->GetLargestPossibleRegion().GetSize();
    FloatImageType::IndexType ind;

    for (unsigned int i = 0; i < croppedImages.GetSize(); i++)
    {
      FloatImagePointer cImg = croppedImages[i];

      cropRescaler->SetInput(cImg);
      cropRescaler->Update();

      FloatImagePointer rImg = cropRescaler->GetOutput();
      for (ind[2] = 0; ind[2] < size[2]; ind[2]++)
        for (ind[1] = 0; ind[1] < size[1]; ind[1]++)
          for (ind[0] = 0; ind[0] < size[0]; ind[0]++)
          {
            cImg->SetPixel(ind, rImg->GetPixel(ind));
          }
    }
  }

  croppedTemplate = 0;

  muLogMacro(<< "Start segmentation...\n");
  typedef EMSegmentationFilter<FloatImageType, FloatImageType> SegFilterType;
  SegFilterType::Pointer segfilter = SegFilterType::New();

  if (debugflag)
    segfilter->DebugOn();

  segfilter->SetInputImages(croppedImages);
  segfilter->SetPriors(croppedPriors);

  segfilter->SetFOVMask(fovmask);

  SegFilterType::VectorType priorweights(4);
  priorweights[0] = emsp->GetPrior1();
  priorweights[1] = emsp->GetPrior2();
  priorweights[2] = emsp->GetPrior3();
  priorweights[3] = emsp->GetPrior4();
  segfilter->SetPriorWeights(priorweights);

  segfilter->SetMaxBiasDegree(emsp->GetMaxBiasDegree());

  segfilter->Update();

  DynArray<std::string> names = emsp->GetImages();

  // Write the labels
  muLogMacro(<< "Writing labels...\n");
  {
    typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
    ByteWriterType::Pointer writer = ByteWriterType::New();

    writer->SetInput(byteCropper->Restore(segfilter->GetOutput()));

    std::string fn = outdir + mu::get_name(names[0].c_str()) + std::string("_labels") + suffstr;
    writer->SetFileName(fn.c_str());
    writer->SetUseCompression(true);
    writer->Update();
  }

  // Write the secondary outputs
  if (writemoreflag)
  {
    muLogMacro(<< "Writing filtered and bias corrected images...\n");
    DynArray<FloatImagePointer> imgset = segfilter->GetCorrected();
    for (unsigned i = 0; i < imgset.GetSize(); i++)
    {
      typedef itk::CastImageFilter<FloatImageType, ShortImageType> CasterType;
      CasterType::Pointer caster = CasterType::New();

      caster->SetInput(cropper->Restore(imgset[i]));
      caster->Update();
      std::string fn =
        outdir + mu::get_name(names[i].c_str()) + std::string("_corrected")
        + suffstr;

      typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
      ShortWriterType::Pointer writer = ShortWriterType::New();

      writer->SetInput(caster->GetOutput());
      writer->SetFileName(fn.c_str());
      writer->SetUseCompression(true);
      writer->Update();
    }

    // Short posteriors
    muLogMacro(<< "Writing posterior images...\n");
    DynArray<ShortImagePointer> probset = segfilter->GetShortPosteriors();
    for (unsigned int i = 0; i < (probset.GetSize()-3); i++)
    {
      typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
      ShortWriterType::Pointer writer = ShortWriterType::New();

      std::string first = outdir + mu::get_name(names[0].c_str());

      std::ostringstream oss;
      oss << first << "_posterior" << i << suffstr << std::ends;

      writer->SetInput(shortCropper->Restore(probset[i]));
      writer->SetFileName(oss.str().c_str());
      writer->SetUseCompression(true);
      writer->Update();
    }

/*
    // Byte posteriors
    muLogMacro(<< "Writing posterior images...\n");
    DynArray<ByteImagePointer> probset = segfilter->GetBytePosteriors();
    for (unsigned i = 0; i < (probset.GetSize()-3); i++)
    {
      typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
      ByteWriterType::Pointer writer = ByteWriterType::New();

      std::string first = outdir + mu::get_name(names[0].c_str());

      std::ostringstream oss;
      oss << first << "_posterior" << i << suffstr << std::ends;

      writer->SetInput(byteCropper->Restore(probset[i]));
      writer->SetFileName(oss.str().c_str());
      writer->SetUseCompression(true);
      writer->Update();
    }

    // Floating point posteriors
    muLogMacro(<< "Writing posterior images...\n");
    DynArray<FloatImagePointer> probset = segfilter->GetPosteriors();
    for (unsigned i = 0; i < (probset.GetSize()-3); i++)
    {
      typedef itk::ImageFileWriter<FloatImageType> FloatWriterType;
      FloatWriterPointer writer = FloatWriterType::New();

      std::string first = outdir + mu::get_name(names[0].c_str());

      std::ostringstream oss;
      oss << first << "_posterior" << i << metasuffstr << std::ends;

      writer->SetInput(cropper->Restore(probset[i]));
      writer->SetFileName(oss.str().c_str());
      writer->SetUseCompression(true);
      writer->Update();
    }
*/

  }

  timer->Stop();

  muLogMacro(<< "All segmentation processes took " << timer->GetElapsedHours() << " hours, ");
  muLogMacro(<< timer->GetElapsedMinutes() << " minutes, ");
  muLogMacro(<< timer->GetElapsedSeconds() << " seconds\n");

  // Close the log file
  (Log::GetInstance())->CloseFile();

  delete timer;

}
