
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCastImageFilter.h"
#include "itkCommand.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMinMaxCurvatureFlowImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkNumericTraits.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkWarpImageFilter.h"

#include "AtlasCropImageSource.h"
#include "MRAnisotropicFilter.h"
#include "NeonateSegmentationFilter.h"
#include "TrimmedRescaleImageModifier.h"

// Use manually instantiated classes for the big program chunks
#define MU_MANUAL_INSTANTIATION
#include "AtlasRegistrationMethod.h"
#include "PairRegistrationMethod.h"
#undef MU_MANUAL_INSTANTIATION

#include "filterFloatImages.h"

#include "NeoSegParameters.h"
#include "NeoSegParametersXMLFile.h"
#include "runNeonateSeg.h"

#include "DynArray.h"
#include "Log.h"
#include "Timer.h"

#include "muException.h"
#include "muFile.h"

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

std::string
CheckImageFormat( std::string format , bool input = false )//if input format, we do not set default to ".nrrd"
{
  std::string ext ;
  format.erase( format.find_last_not_of( " \n\r\t") + 1 ) ;//if XML empty, format is just white spaces. We remove them.
  std::transform( format.begin(), format.end(), format.begin(), ::tolower ) ;//set string to lower case to simplify the comparisons
  if( !format.compare( "analyze" ) )
  {
    ext = ".hdr" ;
  }
  else if( !format.compare( "meta" ) )
  {
    ext = ".mha" ;
  }
  else if( !format.compare( "hdr" )
        || !format.compare( "gipl" )
        || !format.compare( "gipl.gz" )
        || !format.compare( "nrrd" )
        || !format.compare( "nhdr" )
        || !format.compare( "nii.gz" )
        || !format.compare( "nii" )
        || !format.compare( "mha" )
        || !format.compare( "mhd" )
         )
  {
    ext = "." + format ;
  }
  else if( format.compare(" ") && input == false )
  {
    ext = ".nrrd" ;
  }
  else
  {
    std::string str_error = "\'" + format + "\' is not a recognized format. Verify that your type is one of the following:\n\
    analyze\n\
    hdr\n\
    gipl\n\
    gipl.gz\n\
    nrrd\n\
    nhdr\n\
    nii\n\
    nii.gz" ;
    throw str_error ;
  }
  return ext ;
}

void
runNeonateSeg(NeoSegParameters* neop, bool debugflag, bool writemoreflag)
{

  if (!neop->CheckValues())
    muExceptionMacro(<< "Invalid parameters");

  // Create and start a new timer (for the whole process)
  Timer* timer = new Timer();

  // Directory separator string
  std::string separator = std::string("/");
  separator[0] = MU_DIR_SEPARATOR;

  //Get atlas directory
  std::string atlasdir = neop->GetAtlasDirectory();

  // Make sure last character is a separator
  if (atlasdir[atlasdir.size()-1] != MU_DIR_SEPARATOR)
  {
    atlasdir += separator;
  }

  // Get output directory
  std::string outdir = neop->GetOutputDirectory();

  // Make sure last character in output directory string is a separator
  if (outdir[outdir.size()-1] != MU_DIR_SEPARATOR)
  {
    outdir += separator;
  }

  // Create the output directory if not already there
  if (!mu::create_dir(outdir.c_str()))
    return;

  // Set up the logger
  {
    std::string logfn = outdir + neop->GetSuffix() + ".log";
    (Log::GetInstance())->EchoOn();
    (Log::GetInstance())->SetOutputFileName(logfn.c_str());
  }

  // Write out the parameters in XML
  {
    std::string xmlfn = outdir + neop->GetSuffix() + ".xml";
    writeNeoSegParametersXML(xmlfn.c_str(), neop);
  }

  // Set up suffix string for images
  std::string outext = CheckImageFormat( neop->GetOutputFormat() ) ;
  std::string suffstr =
    std::string("_") + std::string(neop->GetSuffix()) + outext;
  std::string metasuffstr =
    std::string("_") + std::string(neop->GetSuffix()) + std::string(".mha");

  muLogMacro(<< "mu::neoseg\n");
  muLogMacro(<< "========================================\n");
  muLogMacro(<< "Program compiled on: " << __DATE__ << "\n");
  muLogMacro(<< "\n");

  // Write input parameters
  muLogMacro(<< "=== Parameters ===\n");
  muLogMacro(<< "Suffix: " << neop->GetSuffix() << "\n");
  muLogMacro(<< "Atlas Directory: " << neop->GetAtlasDirectory() << "\n");
  muLogMacro(<< "Atlas Format: " << neop->GetAtlasFormat() << "\n");
  muLogMacro(<< "Atlas Orientation: " << neop->GetAtlasOrientation() << "\n");
  muLogMacro(<< "Output Directory: " << neop->GetOutputDirectory() << "\n");
  muLogMacro(<< "Output Format: " << neop->GetOutputFormat() << "\n");
  muLogMacro(<< "Input images: \n");
  for (unsigned int i = 0; i < neop->GetImages().GetSize(); i++)
    muLogMacro(<< "  " << "[" << (neop->GetImageOrientations())[i] << "] " <<
      (neop->GetImages())[i] << "\n");
  muLogMacro(
    << "Image filtering parameters: " << neop->GetFilterIterations()
    << " iterations, dt = " << neop->GetFilterTimeStep() << "\n");
  muLogMacro(
    << "Prior weight scales: " << neop->GetPrior1() << ", "
    << neop->GetPrior2() << ", " << neop->GetPrior3()
    << ", " << neop->GetPrior4() << ", " << neop->GetPrior5() << "\n");
  muLogMacro(
    << "Max bias polynomial degree: " << neop->GetMaxBiasDegree() << "\n");
  muLogMacro(<< "Atlas warping: " << neop->GetDoAtlasWarp() << "\n");
  muLogMacro(
    << "Atlas warp spline grid size: " << neop->GetAtlasWarpGridX() << " X "
    << neop->GetAtlasWarpGridY() << " X "
    << neop->GetAtlasWarpGridZ() << "\n");
  muLogMacro(<< "Reference image index = " << neop->GetReferenceImageIndex() << "\n");
  if (neop->UseT1())
    muLogMacro(<< "Reference image is a T1 image\n")
  if (neop->UseT2())
    muLogMacro(<< "Reference image is a T2 image\n")
  muLogMacro(
    << "Prior threshold for sampling = " << neop->GetPriorThreshold() << "\n");
  muLogMacro(
    << "Kernel width (in fraction) = " << neop->GetKernelWidth() << "\n");

// Deprecated parameters
#if 0
  muLogMacro(
    << "Mahalanobis threshold for MCD = " << neop->GetMahalanobisThreshold()
    << "\n");
#endif

  muLogMacro(<< "\n");

  muLogMacro(<< "=== Start ===\n");

  typedef itk::Image<float, 3> FloatImageType;
  typedef itk::Image<unsigned char, 3> ByteImageType;
  typedef itk::Image<short, 3> ShortImageType;

  typedef itk::NearestNeighborInterpolateImageFunction<ByteImageType, double>
    LabelInterpolatorType;

  DynArray<FloatImageType::Pointer> images;
  DynArray<FloatImageType::Pointer> priors;
  ByteImageType::Pointer fovmask;
  FloatImageType::Pointer templateImg;


  {
    typedef AtlasRegistrationMethod<float, float> AtlasRegType;
    AtlasRegType::Pointer atlasreg = AtlasRegType::New();

    atlasreg->SetSuffix(neop->GetSuffix().c_str());

    std::string atlasExt = CheckImageFormat( neop->GetAtlasFormat() , true ) ;

    std::string templatefn = atlasdir + std::string("template");
    if (neop->UseT1())
      templatefn += "T1" + atlasExt ;
    if (neop->UseT2())
      templatefn += "T2" + atlasExt ;
    atlasreg->SetTemplateFileName(templatefn.c_str());

    atlasreg->SetAtlasOrientation(neop->GetAtlasOrientation());

    atlasreg->SetImageFileNames(neop->GetImages());
    atlasreg->SetImageOrientations(neop->GetImageOrientations());
    atlasreg->SetOutputDirectory(outdir);
    // Compute list of file names for the priors
    DynArray<std::string> priorfnlist;
    {
      priorfnlist.Append(atlasdir + std::string("white") + atlasExt);
      priorfnlist.Append(atlasdir + std::string("gray") + atlasExt);
      priorfnlist.Append(atlasdir + std::string("csf") + atlasExt);
      priorfnlist.Append(atlasdir + std::string("rest") + atlasExt);
    }

    atlasreg->SetProbabilityFileNames(priorfnlist);

    atlasreg->SetWarpAtlas(neop->GetDoAtlasWarp());
    atlasreg->SetWarpGridSize(
      neop->GetAtlasWarpGridX(),
      neop->GetAtlasWarpGridY(), 
      neop->GetAtlasWarpGridZ());

    if (debugflag)
      atlasreg->DebugOn();

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

    images = atlasreg->GetImages();
    priors = atlasreg->GetProbabilities();
    fovmask = atlasreg->GetFOVMask();
    templateImg = atlasreg->GetTemplate();

    // Write the registered template
    if (writemoreflag)
    {
      muLogMacro(<< "Writing affine-registered template...\n");

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
        outdir + mu::get_name((neop->GetImages()[0]).c_str()) +
        std::string("_template") + suffstr;

      writer->SetInput(rescaler->GetOutput());
      writer->SetFileName(fn.c_str());
      writer->SetUseCompression(true);
      writer->Update();

// Write white matter prior (for debugging)
#if 0
      rescaler->SetInput(atlasreg->GetProbabilities()[0]);
      rescaler->Update();

      std::string fn2 =
        outdir + mu::get_name((neop->GetImages()[0]).c_str()) +
        std::string("_prior0") + suffstr;

      writer->SetFileName(fn2.c_str());
      writer->SetUseCompression(true);
      writer->Update();
#endif

    }

    // Attempt to read parcellation file and resample
    /*std::string parcfn = atlasdir + std::string("parcellation.gipl");

    muLogMacro(<< "Processing parcellation image...\n");
    try
    {
      typedef itk::ImageFileReader<ByteImageType> ByteReaderType;
      ByteReaderType::Pointer reader = ByteReaderType::New();
      reader->SetFileName(parcfn.c_str());
      reader->Update();

      typedef itk::ResampleImageFilter<ByteImageType, ByteImageType>
        LabelResamplerType;

      LabelResamplerType::Pointer resampler = LabelResamplerType::New();

      FloatImageType::Pointer first = images[0];

      resampler->SetInput(reader->GetOutput());
      resampler->SetTransform(atlasreg->GetTemplateAffineTransform());
      resampler->SetInterpolator(LabelInterpolatorType::New());
      resampler->SetSize(first->GetLargestPossibleRegion().GetSize());
      resampler->SetOutputOrigin(first->GetOrigin());
      resampler->SetOutputSpacing(first->GetSpacing());
      resampler->SetDefaultPixelValue(0);

      resampler->Update();
    
      LabelResamplerType::Pointer warper = LabelResamplerType::New();

      warper->SetInput(resampler->GetOutput());
      warper->SetTransform(atlasreg->GetTemplateBSplineTransform());
      warper->SetInterpolator(LabelInterpolatorType::New());
      warper->SetSize(first->GetLargestPossibleRegion().GetSize());
      warper->SetOutputOrigin(first->GetOrigin());
      warper->SetOutputSpacing(first->GetSpacing());
      warper->SetDefaultPixelValue(0);

      warper->Update();
       
      ByteImageType::Pointer parcellationImg = warper->GetOutput();

      muLogMacro(<< "Writing parcellation image...\n");

      typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
      ByteWriterType::Pointer writer = ByteWriterType::New();

      std::string fn =
        outdir + mu::get_name((neop->GetImages())[0].c_str()) +
        std::string("_parcellation") + suffstr;

      writer->SetInput(parcellationImg);
      writer->SetFileName(fn.c_str());
      writer->SetUseCompression(true);
      writer->Update();
    }
    catch (...)
    {
      muLogMacro(
        << "  Cannot read or register parcellation file: " << parcfn
        << ", skipping\n");
    }*/

  } // end of atlas reg block

  DynArray<std::string> names = neop->GetImages();

  // Write the affine registered images
  /*if (writemoreflag)
  {
    muLogMacro(<< "Writing registered images...\n");
    for (unsigned i = 0; i < images.GetSize(); i++)
    {

      typedef itk::ImageFileWriter<FloatImageType> FloatWriterType;
      FloatWriterType::Pointer floatWriter = FloatWriterType::New();

      std::string fn =
        outdir + mu::get_name(names[i].c_str()) + std::string("_registered_withoutRescaling") +
        suffstr;

      floatWriter->SetInput(images[i]);
      floatWriter->SetFileName(fn.c_str());
      floatWriter->SetUseCompression(true);
      floatWriter->Update();

      typedef itk::RescaleIntensityImageFilter<FloatImageType, ShortImageType>
        ShortRescaleType;
      ShortRescaleType::Pointer rescaler = ShortRescaleType::New();

      rescaler->SetOutputMinimum(0);
      rescaler->SetOutputMaximum(itk::NumericTraits<short>::max());
      rescaler->SetInput(images[i]);
      rescaler->Update();

      typedef itk::ImageFileWriter<ShortImageType> ShortWriterType;
      ShortWriterType::Pointer shortWriter = ShortWriterType::New();

      fn =
        outdir + mu::get_name(names[i].c_str()) + std::string("_registered") +
        suffstr;

      shortWriter->SetInput(rescaler->GetOutput());
      shortWriter->SetFileName(fn.c_str());
      shortWriter->SetUseCompression(true);
      shortWriter->Update();
    }

  }*/

  // Generate cropped images
  /*muLogMacro(<< "Cropping images based on atlas priors...\n");

  typedef AtlasCropImageSource<FloatImageType, FloatImageType> CropSourceType;
  CropSourceType::Pointer cropper = CropSourceType::New();

  cropper->SetPadding(5.0);

  typedef AtlasCropImageSource<ByteImageType, FloatImageType> ByteCropSourceType;
  ByteCropSourceType::Pointer byteCropper = ByteCropSourceType::New();

  byteCropper->SetPadding(5.0);

  DynArray<FloatImageType::Pointer> croppedPriors;
  DynArray<FloatImageType::Pointer> croppedImages;

  {
    DynArray<FloatImageType::Pointer> brainprobs;
    for (unsigned int i = 0; i < (priors.GetSize()-1); i++)
      brainprobs.Append(priors[i]);

    cropper->UseProbabilities(brainprobs);
    byteCropper->UseProbabilities(brainprobs);

    for (unsigned int i = 0; i < images.GetSize(); i++)
    {
      croppedImages.Append(cropper->Crop(images[i]));

      typedef itk::ImageFileWriter<FloatImageType> FloatImageType;
      FloatImageType::Pointer writer = FloatImageType::New();

      std::string fn = outdir + mu::get_name(names[i].c_str()) + std::string("_cropped") + suffstr;

      writer->SetInput(cropper->Crop(images[i]));
      writer->SetFileName(fn.c_str());
      writer->SetUseCompression(true);
      writer->Update();
   }

    for (unsigned int i = 0; i < priors.GetSize(); i++)
      croppedPriors.Append(cropper->Crop(priors[i]));
  }

  // Crop FOV mask
  fovmask = byteCropper->Crop(fovmask);

  typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
  ByteWriterType::Pointer writer = ByteWriterType::New();

  std::string fn = outdir + "/FOVMaskCropped.nrrd";

  writer->SetInput(fovmask);
  writer->SetFileName(fn.c_str());
  writer->SetUseCompression(true);
  writer->Update();*/

  //priors.Clear();
  //images.Clear();

  if (neop->GetFilterIterations() > 0)
  {
    muLogMacro(<< "Non-linear filtering of registered images...\n");

    std::string method = neop->GetFilterMethod();

    Timer* filtertimer = new Timer();

    filterFloatImages(images, method,
      neop->GetFilterIterations(), neop->GetFilterTimeStep());

    filtertimer->Stop();

    muLogMacro(
      << "Non-linear filtering took " << filtertimer->GetElapsedHours()
      << " hours, " << filtertimer->GetElapsedMinutes() << " minutes, "
      << filtertimer->GetElapsedSeconds() << " seconds\n");

    delete filtertimer;
  }

  
  // Rescale the intensities of the cropped images
  muLogMacro(<< "Rescaling intensity of cropped images... ");
  {
    typedef TrimmedRescaleImageModifier<FloatImageType> RescaleModifierType;

    for (unsigned int i = 0; i < images.GetSize(); i++)
    {
      RescaleModifierType::Pointer resm = RescaleModifierType::New();
      resm->SetOutputMinimum(0);
      resm->SetOutputMaximum(itk::NumericTraits<short>::max());
      resm->SetTrimFraction(0.01);

      // Turn both above/below rescale off for no trimming/clamping
      resm->TrimBelowOff();
      resm->TrimAboveOff();

      //resm->TrimAboveOn();

      muLogMacro(<< "1... ");
      resm->Rescale(images[i]);
      muLogMacro(<< "2... ");
    }
  }

  templateImg = 0;


  muLogMacro(<< "Start segmentation...\n");
  NeonateSegmentationFilter::Pointer segfilter =
    NeonateSegmentationFilter::New();

  if (debugflag)
    segfilter->DebugOn();

  //segfilter->SetInputImages(croppedImages);      // with cropping
  //segfilter->SetPriors(croppedPriors);

  segfilter->SetInputImages(images);               // without cropping
  segfilter->SetPriors(priors);

  segfilter->SetFOVMask(fovmask);

  segfilter->SetMahalanobisThreshold(neop->GetMahalanobisThreshold());
  segfilter->SetKernelWidthFraction(neop->GetKernelWidth());
  segfilter->SetPriorThresholdFraction(neop->GetPriorThreshold());

  if (neop->UseT1())
    segfilter->SetReferenceModality(NeonateSegmentationFilter::T1);
  if (neop->UseT2())
    segfilter->SetReferenceModality(NeonateSegmentationFilter::T2);
  segfilter->SetReferenceImageIndex(neop->GetReferenceImageIndex()-1);

  segfilter->SetMyelinPriorWeight(neop->GetPrior1());

  NeonateSegmentationFilter::VectorType priorweights(4);
  priorweights[0] = neop->GetPrior2();
  priorweights[1] = neop->GetPrior3();
  priorweights[2] = neop->GetPrior4();
  priorweights[3] = neop->GetPrior5();
  segfilter->SetPriorWeights(priorweights);

  segfilter->SetMaxBiasDegree(neop->GetMaxBiasDegree());

  segfilter->Update();

  muLogMacro(<< "Writing EM classification labels...\n");
  {
    //ByteImageType::Pointer labelImg =
      //byteCropper->Restore(segfilter->GetOutput());             // with cropping

    ByteImageType::Pointer labelImg = segfilter->GetOutput();     // without cropping

    typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
    ByteWriterType::Pointer writer = ByteWriterType::New();

    std::string fn =
      outdir + mu::get_name(names[0].c_str()) + std::string("_EMonly_labels") +
      suffstr;

    writer->SetInput(labelImg);
    writer->SetFileName(fn.c_str());
    writer->SetUseCompression(true);
    writer->Update();
  }

  // Write the secondary outputs
  if (writemoreflag)
  {

    muLogMacro(<< "Writing filtered and bias corrected images...\n");

    DynArray<FloatImageType::Pointer> imgset = segfilter->GetCorrected();
    for (unsigned i = 0; i < imgset.GetSize(); i++)
    {
      typedef itk::CastImageFilter<FloatImageType, ShortImageType>
        CasterType;
      CasterType::Pointer caster = CasterType::New();

      //caster->SetInput(cropper->Restore(imgset[i]));      // with cropping
      caster->SetInput(imgset[i]);                          // without cropping
      caster->Update();

      typedef itk::ImageFileWriter<ShortImageType> ShortWriterType; 
      ShortWriterType::Pointer writer = ShortWriterType::New();

      std::string fn =
        outdir + mu::get_name(names[i].c_str()) + std::string("_corrected") +
        suffstr;

      writer->SetInput(caster->GetOutput());
      writer->SetFileName(fn.c_str());
      writer->SetUseCompression(true);
      writer->Update();
    }

    muLogMacro(<< "Writing posterior images...\n");

    DynArray<ByteImageType::Pointer> probset = segfilter->GetBytePosteriors();
    for (unsigned i = 0; i < (probset.GetSize()-3); i++)
    {
      ByteImageType::Pointer img = probset[i];

      typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
      ByteWriterType::Pointer writer = ByteWriterType::New();

      std::string first = outdir + mu::get_name(names[0].c_str());
      std::ostringstream oss;
      oss << first << "_EMonly_posterior" << i << suffstr << std::ends;

      //writer->SetInput(byteCropper->Restore(img));    // with cropping
      writer->SetInput(img);                            // without cropping
      writer->SetFileName(oss.str().c_str());
      writer->SetUseCompression(true);
      writer->Update();
    }
  }

/*
  // K-nearest refinement
  muLogMacro(<< "KNN refinement...\n");
  segfilter->DoKNearestNeighbors();
  muLogMacro(<< "Writing KNN classification labels...\n");
  {
    ByteImageType::Pointer labelImg =
      byteCropper->Restore(segfilter->GetOutput());

    typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
    ByteWriterType::Pointer writer = ByteWriterType::New();

    std::string fn =
      outdir + mu::get_name(names[0].c_str()) + std::string("_labels_knn") +
      suffstr;

    writer->SetInput(labelImg);
    writer->SetFileName(fn.c_str());
    writer->SetUseCompression(true);
    writer->Update();
  }
*/

/*
  // Kernel density refinement
  muLogMacro(<< "Non-parametric kernel density refinement...\n");
  segfilter->DoKernelDensityEstimation();
  muLogMacro(<< "Writing kernel density classification labels...\n");
  {
    ByteImageType::Pointer labelImg =
      byteCropper->Restore(segfilter->GetOutput());

    typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
    ByteWriterType::Pointer writer = ByteWriterType::New();

    std::string fn =
      outdir + mu::get_name(names[0].c_str()) + std::string("_parzen_labels") +
      suffstr;

    writer->SetInput(labelImg);
    writer->SetFileName(fn.c_str());
    writer->SetUseCompression(true);
    writer->Update();
  }
  // Write Parzen posteriors
  if (writemoreflag)
  {
    muLogMacro(<< "Writing Parzen posterior images...\n");
    DynArray<ByteImageType::Pointer> probset = segfilter->GetBytePosteriors();
    for (unsigned i = 0; i < (probset.GetSize()-3); i++)
    {
      ByteImageType::Pointer img = probset[i];

      typedef itk::ImageFileWriter<ByteImageType> ByteWriterType;
      ByteWriterType::Pointer writer = ByteWriterType::New();

      std::string first = outdir + mu::get_name(names[0].c_str());
      std::ostringstream oss;
      oss << first << "_parzen_posterior" << i << suffstr << std::ends;

      writer->SetInput(byteCropper->Restore(img));
      writer->SetFileName(oss.str().c_str());
      writer->SetUseCompression(true);
      writer->Update();
    }
  }
*/
  timer->Stop();

  muLogMacro(
    << "All segmentation processes took "
    << timer->GetElapsedHours() << " hours, "
    << timer->GetElapsedMinutes() << " minutes, "
    << timer->GetElapsedSeconds() << " seconds\n");

  // Close the log file
  (Log::GetInstance())->CloseFile();

  delete timer;

}
