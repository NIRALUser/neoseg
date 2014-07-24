
#include "itkCurvatureFlowImageFilter.h"

#include "MRAnisotropicFilter.h"

#include "filterFloatImages.h"

void
filterFloatImages(
  DynArray<itk::Image<float, 3>::Pointer>& images,
  std::string& method,
  unsigned int iters,
  double dt)
{
  typedef itk::Image<float, 3> FloatImageType;

  if (method.compare("Grad aniso diffusion") == 0)
  {
std::cout << "Grad aniso" << std::endl;
    typedef MRAnisotropicFilter<FloatImageType> AnisoFilterType;
    AnisoFilterType::Pointer anisofilt = AnisoFilterType::New();

    //if (debugflag)
    //  anisofilt->DebugOn();

    anisofilt->SetNumberOfIterations(iters);
    anisofilt->SetTimeStep(dt);

    anisofilt->Filter(images);
  }
  // else if? Default is curvature flow
  else
  {
std::cout << "K flow" << std::endl;
    typedef itk::CurvatureFlowImageFilter<
      FloatImageType, FloatImageType> CurvatureFilterType;

    for (unsigned int k = 0; k < images.GetSize(); k++)
    { 
      CurvatureFilterType::Pointer cfilt = CurvatureFilterType::New();
      
      cfilt->SetNumberOfIterations(iters);
      cfilt->SetTimeStep(dt);
      cfilt->SetInput(images[k]);
      cfilt->Update();
      
      images[k] = cfilt->GetOutput();
    }
  }
}
