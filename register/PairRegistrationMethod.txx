
#ifndef _PairRegistrationMethod_txx
#define _PairRegistrationMethod_txx

#include "PairRegistrationMethod.h"

#include "itkCommand.h"
#include "itkDemonsRegistrationFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMultiResolutionImageRegistrationMethod.h"

#include "itkLBFGSBOptimizer.h"

#include "vnl/vnl_math.h"

#include "PairRegistrationMethod.h"
#include "NegativeMIImageMatchMetric.h"
#include "GradientDescentOptimizer.h"
#include "PowellOptimizer.h"
#include "RegistrationParameters.h"

#include "DynArray.h"
#include "Log.h"
#include "muException.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

#include <math.h>
#include <stdlib.h>

#ifndef LINE_MAX
#define LINE_MAX 1024
#endif

/*
// Observer for the amoeba optimizer
class AmoebaIterationUpdate : public itk::Command
{
public:
  typedef  AmoebaIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
  typedef itk::AmoebaOptimizer OptimizerType;
  typedef const OptimizerType* OptimizerPointer;
protected:
  AmoebaIterationUpdate() { m_CurrentIteration = 0; };

  unsigned int m_CurrentIteration;
public:
  void ResetIterationCounter() { m_CurrentIteration = 0; }
  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }
  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >(object);
    if (!itk::IterationEvent().CheckEvent(&event))
    {
      return;
    }
    m_CurrentIteration++;
    muLogMacro(
      << "  Iter: " << m_CurrentIteration << " ||  "
      << "-MI: " << optimizer->GetCachedValue() << "\n");
    muLogMacro(<< "  " << optimizer->GetCachedCurrentPosition() << "\n");
  }
};
*/

// Observer for the affine optimizer iterations
class AffineIterationUpdate : public itk::Command
{
public:
  typedef  AffineIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
  //typedef AffineOptimizer OptimizerType;
  typedef PowellOptimizer OptimizerType;
  typedef const OptimizerType* OptimizerPointer;
protected:
  AffineIterationUpdate() { }
public:
  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }
  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >(object);
    if (!itk::IterationEvent().CheckEvent(&event))
    {
      return;
    }
    muLogMacro(
      << "  Iter: " << optimizer->GetCurrentIteration() << " ||  "
      << "-MI: " << optimizer->GetValue() << "\n");
    muLogMacro(<< "  " << optimizer->GetCurrentPosition() << "\n");
  }
};

// Update parameters at change of resolution level
template <class TRegistration>
class AffineLevelUpdate : public itk::Command
{
public:
  typedef  AffineLevelUpdate   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );
protected:
  AffineLevelUpdate() {};
public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  //typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
  //typedef   itk::LBFGSBOptimizer   OptimizerType;
  typedef   PowellOptimizer   OptimizerType;
  //typedef   GradientDescentOptimizer   OptimizerType;
  typedef   OptimizerType *                            OptimizerPointer;
  void Execute(itk::Object * object, const itk::EventObject & event)
  {
    if( typeid( event ) != typeid( itk::IterationEvent ) )
    {
      return;
    }
    RegistrationPointer registration =
                        dynamic_cast<RegistrationPointer>( object );
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >(
                       registration->GetOptimizer() );

    muLogMacro(<< "  Registration at level "
      << registration->GetCurrentLevel() + 1 << "\n");

    if ( registration->GetCurrentLevel() == 0 )
    {
      //optimizer->SetMaximumStepLength(10.0);
      //optimizer->SetMinimumStepLength(1.0);
    }
    else
    {
      //optimizer->SetLearningRate(optimizer->GetLearningRate() / 2.0);
      //optimizer->SetNumberOfIterations(
      //  optimizer->GetNumberOfIterations() +
      //  registration->GetCurrentLevel()*500);
      //optimizer->SetMaximumStepLength(optimizer->GetCurrentStepLength());
      //optimizer->SetMinimumStepLength(optimizer->GetMinimumStepLength() / 10.0);
      //optimizer->SetStepLength(optimizer->GetStepLength() / 10.0);
    }
  }
  void Execute(const itk::Object * , const itk::EventObject & )
  { return; }
};

// Observer for the BSpline deformable registration
class BSplineIterationUpdate : public itk::Command
{
public:
  typedef  BSplineIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  BSplineIterationUpdate() {};
public:
  typedef itk::LBFGSBOptimizer     OptimizerType;
  typedef const OptimizerType*    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer =
      dynamic_cast< OptimizerPointer >( object );
    if( typeid( event ) != typeid( itk::IterationEvent ) )
    {
      return;
    }
    muLogMacro(
      << "  B-spline iter: " << optimizer->GetCurrentIteration()+1 << " ||  "
      << "-MI: " << optimizer->GetValue() << "\n");
    //muLogMacro(<< optimizer->GetInfinityNormOfProjectedGradient() << "\n");
  }
};


// Observer for the Demons deformable registration
class DemonsIterationUpdate : public itk::Command
{
public:
  typedef  DemonsIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro(Self);
protected:
  DemonsIterationUpdate() { m_Iterations = 0; };

  typedef itk::Image<float, 3> InternalImageType;
  typedef itk::Vector<float, 3> VectorPixelType;
  typedef itk::Image<VectorPixelType, 3> DeformationFieldType;

  typedef itk::DemonsRegistrationFilter<
    InternalImageType,
    InternalImageType,
    DeformationFieldType> RegistrationFilterType;

private:
  unsigned int m_Iterations;

public:
  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }
  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    const RegistrationFilterType * filter =
      dynamic_cast< const RegistrationFilterType * >( object );
    if( typeid( event ) != typeid( itk::IterationEvent ) )
    {
      return;
    }
    m_Iterations++;
    muLogMacro(
      << "  Demons iteration " << m_Iterations << ", metric = "
      << filter->GetMetric() << "\n");
  }
};

template <class TPixel>
void
PairRegistrationMethod<TPixel>
::ReadNextLine(char* s, std::ifstream& infile)
{
  while (true)
  {
    infile.getline(s, LINE_MAX);

    if (infile.fail())
      break;

    if (s[0] != '#')
      break;
  }
}

template <class TPixel>
PairRegistrationMethod<TPixel>::AffineTransformType::Pointer
PairRegistrationMethod<TPixel>
::RegisterAffine(ImageType* fixedImg, ImageType* movingImg)
{
  if (fixedImg == NULL || movingImg == NULL)
    muExceptionMacro(<< "One of input images is NULL");

  // Get image info
  typename ImageType::SizeType fixedSize =
    fixedImg->GetLargestPossibleRegion().GetSize();
  typename ImageType::SizeType movingSize =
    movingImg->GetLargestPossibleRegion().GetSize();

  typename ImageType::SpacingType fixedSpacing = fixedImg->GetSpacing();
  typename ImageType::SpacingType movingSpacing = movingImg->GetSpacing();

  AffineTransformType::CenterType fixedCenter;
  fixedCenter[0] = (fixedSize[0] - 1.0) * fixedSpacing[0] / 2.0;
  fixedCenter[1] = (fixedSize[1] - 1.0) * fixedSpacing[1] / 2.0;
  fixedCenter[2] = (fixedSize[2] - 1.0) * fixedSpacing[2] / 2.0;

  AffineTransformType::CenterType movingCenter;
  movingCenter[0] = (movingSize[0] - 1.0) * movingSpacing[0] / 2.0;
  movingCenter[1] = (movingSize[1] - 1.0) * movingSpacing[1] / 2.0;
  movingCenter[2] = (movingSize[2] - 1.0) * movingSpacing[2] / 2.0;

  // Define framework
  //typedef PowellOptimizer OptimizerType;
  //typedef GradientDescentOptimizer GradientOptimizerType;
  typedef itk::LinearInterpolateImageFunction<
    ImageType, double> InterpolatorType;
  typedef NegativeMIImageMatchMetric<ImageType, ImageType> MetricType;    
  //typedef itk::MattesMutualInformationImageToImageMetric<
  //  ImageType, ImageType> MetricType;
  
  typedef itk::MultiResolutionImageRegistrationMethod<
    ImageType, ImageType> RegistrationType;

  // Create objects
  PowellOptimizer::Pointer optimizer = PowellOptimizer::New();
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  typename RegistrationType::Pointer registration = RegistrationType::New();
  typename MetricType::Pointer metric = MetricType::New();

  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(interpolator);
  registration->SetMetric(metric);

  interpolator->SetInputImage(movingImg);

  AffineTransformType::Pointer affine = AffineTransformType::New();
  affine->SetSourceCenter(fixedCenter[0], fixedCenter[1], fixedCenter[2]);
  affine->SetTargetCenter(movingCenter[0], movingCenter[1], movingCenter[2]);

  // Use initial identity transform (centered at mid-image)
  AffineTransformType::ParametersType initp = affine->GetParameters();

  registration->SetTransform(affine);
  registration->SetInitialTransformParameters(initp);

  registration->SetFixedImage(fixedImg);
  registration->SetMovingImage(movingImg);

  registration->SetFixedImageRegion(fixedImg->GetLargestPossibleRegion());

  // Set steps for optimization
  {
    PowellOptimizer::ParametersType steps(12);
    steps[0] = MU_AFFINE_STEP_TRANSLATE;
    steps[1] = MU_AFFINE_STEP_TRANSLATE;
    steps[2] = MU_AFFINE_STEP_TRANSLATE;
    steps[3] = MU_AFFINE_STEP_ROTATE;
    steps[4] = MU_AFFINE_STEP_ROTATE;
    steps[5] = MU_AFFINE_STEP_ROTATE;
    steps[6] = MU_AFFINE_STEP_SCALE;
    steps[7] = MU_AFFINE_STEP_SCALE;
    steps[8] = MU_AFFINE_STEP_SCALE;
    steps[9] = MU_AFFINE_STEP_SKEW;
    steps[10] = MU_AFFINE_STEP_SKEW;
    steps[11] = MU_AFFINE_STEP_SKEW;
    optimizer->SetInitialSteps(steps);

    PowellOptimizer::OrderType order(12);
    order[0] = MU_AFFINE_ORDER0;
    order[1] = MU_AFFINE_ORDER1;
    order[2] = MU_AFFINE_ORDER2;
    order[3] = MU_AFFINE_ORDER3;
    order[4] = MU_AFFINE_ORDER4;
    order[5] = MU_AFFINE_ORDER5;
    order[6] = MU_AFFINE_ORDER6;
    order[7] = MU_AFFINE_ORDER7;
    order[8] = MU_AFFINE_ORDER8;
    order[9] = MU_AFFINE_ORDER9;
    order[10] = MU_AFFINE_ORDER10;
    order[11] = MU_AFFINE_ORDER11;
    optimizer->SetOrder(order);
  }

  optimizer->SetMaximumIterations(15);

/*
  typename MetricType::ParametersType derivSteps(12);
  derivSteps.Fill(1e-8);
  for (int i = 0; i < 3; i++)
    derivSteps[i] = 1e-4;
  metric->SetDerivativeStepLengths(derivSteps);
*/

  metric->SetNumberOfBins(255);
  metric->SetNormalized(false);
  metric->SetUseKMeans(false);

/*
  metric->SetNumberOfHistogramBins(128);
  unsigned int numSamples =
    fixedImg->GetLargestPossibleRegion().GetNumberOfPixels() / 5;
std::cout << "20% of #pixels = " << numSamples << std::endl;
  if (numSamples > 100000)
    numSamples = 100000;
//TODO: need to be adjusted based on pyramid level
  metric->SetNumberOfSpatialSamples(numSamples);
  metric->ReinitializeSeed( 76926294 );
*/

  // Create the Command observer and register it with the optimizer
  AffineIterationUpdate::Pointer observer = AffineIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  // Create the Command interface observer and register it with the optimizer.
  typedef AffineLevelUpdate<RegistrationType> LevelUpdaterType;
  typename LevelUpdaterType::Pointer levelUpd = LevelUpdaterType::New();
  registration->AddObserver(itk::IterationEvent(), levelUpd);

  muLogMacro(<< "Beginning affine registration...\n");
#if 0
  // Use ITK framework to handle multi resolution registration
  metric->SetSampleSpacing(1.0);

  registration->SetNumberOfLevels(3);

#if ITK_VERSION_MAJOR <= 4 and ITK_VERSION_MINOR <= 8
  registration->StartRegistration();
#else
  registration->Update();
#endif

  affine->SetParameters(registration->GetLastTransformParameters());
#else
  // Manage the multi resolution registration here
  metric->SetFixedImage(fixedImg);
  metric->SetMovingImage(movingImg);
  metric->SetTransform(affine);

  optimizer->SetCostFunction(metric);
  optimizer->SetInitialPosition(affine->GetParameters());

  muLogMacro(<< "Registering at [2x2x2]...\n");
  metric->SetSampleSpacing(2.0);
  optimizer->StartOptimization();

  muLogMacro(<< "Registering at [1x1x1]...\n");
  metric->SetSampleSpacing(1.0);
  optimizer->SetInitialPosition(optimizer->GetCurrentPosition());
  optimizer->StartOptimization();

  affine->SetParameters(optimizer->GetCurrentPosition());
#endif

  muLogMacro(<< "Done with affine registration\n");

  return affine;
}

template <class TPixel>
PairRegistrationMethod<TPixel>::BSplineTransformType::Pointer
PairRegistrationMethod<TPixel>
::RegisterBSpline(ImageType* fixedImg, ImageType* movingImg,
  unsigned int nx, unsigned int ny, unsigned int nz)
{
  typedef itk::MattesMutualInformationImageToImageMetric<
    ImageType, ImageType> MIMetricType;

  typedef itk::LinearInterpolateImageFunction<ImageType, double>
    InterpolatorType;

  typedef itk::LBFGSBOptimizer OptimizerType;

  typedef itk::ImageRegistrationMethod<ImageType, ImageType>
    RegistrationType;

  typename MIMetricType::Pointer metric = MIMetricType::New();
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  OptimizerType::Pointer optimizer = OptimizerType::New();

  typename RegistrationType::Pointer registration = RegistrationType::New();
  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetInterpolator(interpolator);

  interpolator->SetInputImage(movingImg);

  BSplineTransformType::Pointer btrafo = BSplineTransformType::New();
  registration->SetTransform(btrafo);

  registration->SetFixedImage(fixedImg);
  registration->SetMovingImage(movingImg);

  registration->SetFixedImageRegion(fixedImg->GetLargestPossibleRegion());

  double gridSize[3] = {(double)nx,(double)ny,(double)nz};

  typedef BSplineTransformType::RegionType RegionType;
  RegionType bsplineRegion;
  // Pad grid for border regions
  BSplineTransformType::SizeType totalGridSize;
  for (unsigned int i = 0; i < 3; i++)
    totalGridSize[i] = gridSize[i] + 3;
  bsplineRegion.SetSize(totalGridSize);

  typedef BSplineTransformType::SpacingType SpacingType;
  SpacingType spacing = fixedImg->GetSpacing();

  typedef BSplineTransformType::OriginType OriginType;
  OriginType origin = fixedImg->GetOrigin();

  BSplineTransformType::SizeType fixedImgSize =
    fixedImg->GetLargestPossibleRegion().GetSize();

  for (unsigned int i = 0; i < 3; i++)
  {
    spacing[i] *= floor( static_cast<double>(fixedImgSize[i] - 1)  /
                  static_cast<double>(gridSize[i] - 1) );
    origin[i]  -=  spacing[i];
  }

  btrafo->SetGridSpacing(spacing);
  btrafo->SetGridOrigin(origin);
  btrafo->SetGridRegion(bsplineRegion);

  typedef BSplineTransformType::ParametersType ParametersType;

  unsigned int numParams = btrafo->GetNumberOfParameters();

  ParametersType initp(numParams);
  initp.Fill(0.0);

  // Force assignment by value, otherwise will need to maintain p
  // outside of this function
  btrafo->SetParametersByValue(initp);

  registration->SetInitialTransformParameters(initp);

  OptimizerType::BoundSelectionType boundSelect(numParams);
  OptimizerType::BoundValueType upperBound(numParams);
  OptimizerType::BoundValueType lowerBound(numParams);

  boundSelect.Fill( 0 );
  upperBound.Fill( 0.0 );
  lowerBound.Fill( 0.0 );

  optimizer->SetBoundSelection( boundSelect );
  optimizer->SetUpperBound( upperBound );
  optimizer->SetLowerBound( lowerBound );

  optimizer->SetCostFunctionConvergenceFactor(1e+7);
  optimizer->SetProjectedGradientTolerance(1e-4);
  optimizer->SetMaximumNumberOfIterations(400);
  optimizer->SetMaximumNumberOfEvaluations(800);
  optimizer->SetMaximumNumberOfCorrections(5);

  BSplineIterationUpdate::Pointer obs = BSplineIterationUpdate::New();
  optimizer->AddObserver(itk::IterationEvent(), obs);

  metric->SetNumberOfHistogramBins(128);
  metric->SetNumberOfSpatialSamples(
    fixedImg->GetLargestPossibleRegion().GetNumberOfPixels() / 8);

  metric->ReinitializeSeed( 76926294 );

  #if ITK_VERSION_MAJOR <= 4 and ITK_VERSION_MINOR <= 8
    registration->StartRegistration();
  #else
    registration->Update();
  #endif

  btrafo->SetParametersByValue(registration->GetLastTransformParameters());

  return btrafo;
}

template <class TPixel>
PairRegistrationMethod<TPixel>::DeformationFieldType::Pointer
PairRegistrationMethod<TPixel>
::RegisterDemons(
  ImageType* fixedImg, ImageType* movingImg, unsigned int numIters)
{
  typedef itk::HistogramMatchingImageFilter<ImageType, ImageType>
    MatchFilterType;
  typename MatchFilterType::Pointer matcher = MatchFilterType::New();

  matcher->SetInput(movingImg);
  matcher->SetReferenceImage(fixedImg);

  matcher->SetNumberOfHistogramLevels(1024);
  matcher->SetNumberOfMatchPoints(8);
  matcher->ThresholdAtMeanIntensityOn();

  matcher->Update();

  typedef itk::DemonsRegistrationFilter<
    ImageType, ImageType, DeformationFieldType> DemonsType;

  typename DemonsType::Pointer demons = DemonsType::New();

  demons->AddObserver(itk::IterationEvent(), DemonsIterationUpdate::New());

  demons->SetFixedImage(fixedImg);
  demons->SetMovingImage(matcher->GetOutput());

  demons->SetNumberOfIterations(numIters);
  demons->SetStandardDeviations(1.0);

  demons->Update();

  return demons->GetOutput();
}

template <class TPixel>
PairRegistrationMethod<TPixel>::AffineTransformType::Pointer
PairRegistrationMethod<TPixel>
::ReadAffineTransform(const char* fn)
{
  std::ifstream infile;
  infile.open(fn);

  if (infile.fail())
    muExceptionMacro(<< "Failed opening " << fn);

  AffineTransformType::Pointer trafo = AffineTransformType::New();

  AffineTransformType::ParametersType p = trafo->GetParameters();
  for (unsigned int i = 0; i < 12; i++)
  {
    char s[LINE_MAX];
    ReadNextLine(s, infile);
    p[i] = atof(s);
  }

  AffineTransformType::CenterType sourceC;
  for (unsigned int i = 0; i < 3; i++)
  {
    char s[LINE_MAX];
    ReadNextLine(s, infile);
    sourceC[i] = atof(s);
  }

  AffineTransformType::CenterType targetC;
  for (unsigned int i = 0; i < 3; i++)
  {
    char s[LINE_MAX];
    ReadNextLine(s, infile);
    targetC[i] = atof(s);
  }

  bool isForwardOrder;
  {
    char s[LINE_MAX];
    ReadNextLine(s, infile);
    isForwardOrder= atoi(s);
  }

  infile.close();

  trafo->SetAllParameters(
    p,
    sourceC[0], sourceC[1], sourceC[2],
    targetC[0], targetC[1], targetC[2],
    isForwardOrder);

  return trafo;
}

template <class TPixel>
void
PairRegistrationMethod<TPixel>
::WriteAffineTransform(const char* fn, const AffineTransformType* affine)
{
  std::ofstream outfile;
  outfile.open(fn);
  outfile.setf(std::ios::fixed, std::ios::floatfield);
  outfile.precision(50);

  if (outfile.fail())
    muExceptionMacro(<< "Failed opening " << fn);

  outfile << "# Chained affine transform" << std::endl;
  outfile << "#" << std::endl;

  DynArray<std::string> commentsList;
  commentsList.Append("Translations:");
  commentsList.Append("Rotations:");
  commentsList.Append("Scalings:");
  commentsList.Append("Skews:");

  AffineTransformType::ParametersType p = affine->GetParameters();
  for (unsigned int i = 0; i < 12; i++)
  {
    if ((i % 3) == 0)
      outfile << "# " << commentsList[i/3] << std::endl;
    outfile << p[i] << std::endl;
  }

  outfile << "# Center of rotation (source): " << std::endl;
  AffineTransformType::CenterType sourceC = affine->GetSourceCenter();
  for (unsigned int i = 0; i < 3; i++)
    outfile << sourceC[i] << std::endl;

  outfile << "# Center of rotation (target): " << std::endl;
  AffineTransformType::CenterType targetC = affine->GetTargetCenter();
  for (unsigned int i = 0; i < 3; i++)
    outfile << targetC[i] << std::endl;

  outfile << "# Forward composition order?" << std::endl;
  outfile << (int)affine->IsForwardEvaluation() << std::endl;

  outfile.close();
}

template <class TPixel>
PairRegistrationMethod<TPixel>::BSplineTransformType::Pointer
PairRegistrationMethod<TPixel>
::ReadBSplineTransform(const char* fn)
{
  std::ifstream infile;
  infile.open(fn);

  if (infile.fail())
    muExceptionMacro(<< "Failed reading " << fn);

  BSplineTransformType::Pointer btrafo = BSplineTransformType::New();

  typedef BSplineTransformType::RegionType RegionType;
  RegionType bsplineRegion;

  // Read size
  BSplineTransformType::SizeType size;
  for (unsigned int i = 0; i < 3; i++)
  {
    char s[LINE_MAX];
    ReadNextLine(s, infile);
    size[i] = atoi(s);
  }

  // Pad
  for (unsigned int i = 0; i < 3; i++)
  {
    size[i] += 2;
  }

  bsplineRegion.SetSize(size);

  // Read spacing
  BSplineTransformType::SpacingType spacing;
  for (unsigned int i = 0; i < 3; i++)
  {
    char s[LINE_MAX];
    ReadNextLine(s, infile);
    spacing[i] = atof(s);
  }

  // Assume default origin
  typedef BSplineTransformType::OriginType OriginType;
  OriginType origin;
  for (unsigned int i = 0; i < 3; i++)
    origin[i] = -spacing[i];

  btrafo->SetGridSpacing(spacing);
  btrafo->SetGridOrigin(origin);
  btrafo->SetGridRegion(bsplineRegion);

  // Read parameters
  unsigned int numParams = btrafo->GetNumberOfParameters();

  BSplineTransformType::ParametersType p(numParams);
  for (unsigned int i = 0; i < numParams; i++)
  {
    char s[LINE_MAX];
    ReadNextLine(s, infile);
    p[i] = atof(s);
  }
  btrafo->SetParametersByValue(p);

  infile.close();

  return btrafo;
}

template <class TPixel>
void
PairRegistrationMethod<TPixel>
::WriteBSplineTransform(const char* fn, const BSplineTransformType* btrafo)
{
  if (btrafo == 0)
    muExceptionMacro(<< "NULL B-spline");

  if (fn == 0)
    muExceptionMacro(<< "NULL file name");

  std::ofstream outfile;
  outfile.open(fn);
  outfile.setf(std::ios::fixed, std::ios::floatfield);
  outfile.precision(50);

  if (outfile.fail())
    muExceptionMacro(<< "Error writing file: " << fn);

  outfile << "# B-spline warp transform" << std::endl;
  outfile << "#" << std::endl;

  // Write size
  outfile << "# Grid size:" << std::endl;
  BSplineTransformType::SizeType size = btrafo->GetGridRegion().GetSize();
  for (unsigned int i = 0; i < 3; i++)
  {
    outfile << size[i]-2 << std::endl;
  }

  // Write spacing
  outfile << "# Grid spacings:" << std::endl;
  BSplineTransformType::SpacingType spacing = btrafo->GetGridSpacing();
  for (unsigned int i = 0; i < 3; i++)
  {
    outfile << spacing[i] << std::endl;
  }

  outfile << "# B-spline coefficients:" << std::endl;
  BSplineTransformType::ParametersType p = btrafo->GetParameters();
  for (unsigned int i = 0; i < p.GetSize(); i++)
  {
    outfile << p[i] << std::endl;
  }

  outfile.close();
}

template <class TPixel>
PairRegistrationMethod<TPixel>::DeformationFieldType::Pointer
PairRegistrationMethod<TPixel>
::ReadDeformationField(const char* fn)
{
  typedef itk::ImageFileReader<DeformationFieldType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName(fn);
  reader->Update();

  DeformationFieldType::Pointer ret = reader->GetOutput();

  return ret;
}

template <class TPixel>
void
PairRegistrationMethod<TPixel>
::WriteDeformationField(const char* fn, const DeformationFieldType* def)
{
  typedef itk::ImageFileWriter<DeformationFieldType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetInput(def);
  writer->SetFileName(fn);
  writer->Update();
}


#endif
