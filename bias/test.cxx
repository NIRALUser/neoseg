
#include "LLSBiasCorrector.h"

#include "itkDiscreteGaussianImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkOutputWindow.h"
#include "itkTextOutput.h"

using namespace itk;

#include <iostream>

#include <math.h>

int main(int argc, char** argv)
{

  // Use text output
  itk::TextOutput::Pointer textout = itk::TextOutput::New();
  itk::OutputWindow::SetInstance(textout);

  typedef Image<float, 3> FloatImageType;
  typedef Image<unsigned char, 3> ByteImageType;

  typedef LLSBiasCorrector<FloatImageType, FloatImageType> BiasCorrectorType;

  typedef itk::RescaleIntensityImageFilter<FloatImageType, ByteImageType>
    ConverterType;

  typedef ImageFileWriter<ByteImageType> WriterType;

  // Generate images
  FloatImageType::IndexType index;

  FloatImageType::RegionType region;
  FloatImageType::SizeType size;

  size[0] = 128;
  size[1] = 128;
  size[2] = 128;
  region.SetSize(size);

  FloatImageType::SpacingType spacing;
  spacing[0] = 1.0;
  spacing[1] = 1.0;
  spacing[2] = 1.0;

  // Image data (cube)
  std::cout << "Generate test image..." << std::endl;
  FloatImageType::Pointer img = FloatImageType::New();

  img->SetRegions(region);
  img->SetSpacing(spacing);
  img->Allocate();

  FloatImageType::Pointer corrImg = FloatImageType::New();

  corrImg->SetRegions(region);
  corrImg->SetSpacing(spacing);
  corrImg->Allocate();

  for (index[0] = 0; index[0] < size[0]; index[0]++)
    for (index[1] = 0; index[1] < size[1]; index[1]++)
      for (index[2] = 0; index[2] < size[2]; index[2]++)
      {
        if ((index[0] >= 32 && index[0] <= 96)
          &&
          (index[1] >= 32 && index[1] <= 96)
          &&
          (index[2] >= 32 && index[2] <= 96))
          img->SetPixel(index, 200.0);
        else
          img->SetPixel(index, 150.0);
      }

  // Foreground probability
  std::cout << "Generating probabilities..." << std::endl;
  FloatImageType::Pointer fgprob = FloatImageType::New();

  fgprob->SetRegions(region);
  fgprob->Allocate();

  fgprob->SetSpacing(spacing);

  fgprob->FillBuffer(0);

  for (index[0] = 0; index[0] < size[0]; index[0]++)
    for (index[1] = 0; index[1] < size[1]; index[1]++)
      for (index[2] = 0; index[2] < size[2]; index[2]++)
      {
        if (index[0] < 32 || index[0] > 96)
          continue;
        if (index[1] < 32 || index[1] > 96)
          continue;
        if (index[2] < 32 || index[2] > 96)
          continue;
        fgprob->SetPixel(index, 1.0);
      }

  // Background probabilities
  FloatImageType::Pointer bgprob = FloatImageType::New();

  bgprob->SetRegions(region);
  bgprob->Allocate();

  bgprob->SetSpacing(spacing);

  bgprob->FillBuffer(0);

  for (index[0] = 0; index[0] < size[0]; index[0]++)
    for (index[1] = 0; index[1] < size[1]; index[1]++)
      for (index[2] = 0; index[2] < size[2]; index[2]++)
      {
        if ((index[0] >= 32 && index[0] <= 96)
          &&
          (index[1] >= 32 && index[1] <= 96)
          &&
          (index[2] >= 32 && index[2] <= 96))
          continue;
        bgprob->SetPixel(index, 1.0);
      }

  BiasCorrectorType::MaskImageType::Pointer maskImg = 
    BiasCorrectorType::MaskImageType::New();
  maskImg->SetRegions(region);
  maskImg->Allocate();
  maskImg->SetSpacing(spacing);
  maskImg->FillBuffer(1);

  // Blur probabilities
  std::cout << "Blurring probabilities..." << std::endl;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType, FloatImageType>
    BlurType;

  BlurType::Pointer blur1 = BlurType::New();
  blur1->SetInput(fgprob);
  blur1->SetVariance(0.5);
  blur1->Update();

  BlurType::Pointer blur2 = BlurType::New();
  blur2->SetInput(bgprob);
  blur2->SetVariance(0.5);
  blur2->Update();

  // Generate and apply bias field
  std::cout << "Applying bias..." << std::endl;

  unsigned int maxdegree = 2;
  unsigned int numcoeffs = (maxdegree+1)*(maxdegree+2)*(maxdegree+3)/6;
  //unsigned int numcoeffs = 10;

  double coeffs[numcoeffs];
  for (unsigned int i = 0; i < numcoeffs; i++)
    coeffs[i] = 0.0;

  coeffs[0] = 0.01;
  coeffs[1] = 0.03;
  coeffs[2] = -0.14;
  coeffs[3] = 0.046;
  coeffs[4] = 0.09;
  coeffs[5] = -0.04;
  coeffs[6] = 0.2005;
  coeffs[7] = 0.003;
  coeffs[8] = 0.0115;
  coeffs[9] = -0.0128;

  double xc, yc, zc;
  double xmid = 64;
  double ymid = 64;
  double zmid = 64;

  for (index[0] = 0; index[0] < size[0]; index[0]++)
    for (index[1] = 0; index[1] < size[1]; index[1]++)
      for (index[2] = 0; index[2] < size[2]; index[2]++)
      {
        double logbias = 0;

        unsigned int c = 0;

        for (int order = 0; order <= maxdegree; order++) {
          for (int xorder = 0; xorder <= order; xorder++) {
            for (int yorder = 0; yorder <= (order-xorder); yorder++) {

              int zorder = order - xorder - yorder;

              xc = (index[0] - xmid) / xmid;
              yc = (index[1] - ymid) / ymid;
              zc = (index[2] - zmid) / zmid;

              double poly =
                (double)(pow(xc,xorder) * pow(yc,yorder) * pow(zc,zorder));

              logbias += coeffs[c] * poly;

              c++;

            }
          }
        }

        double bias = exp(logbias);

/*
#define TEST_BIAS_CLAMP_MIN 0.6667
#define TEST_BIAS_CLAMP_MAX 1.5

        if (bias < TEST_BIAS_CLAMP_MIN)
          bias = TEST_BIAS_CLAMP_MIN;
        if (bias > TEST_BIAS_CLAMP_MAX)
          bias = TEST_BIAS_CLAMP_MAX;
*/

        img->SetPixel(index, img->GetPixel(index)*bias);

      }


  // Write image
  std::cout << "Write bias-ed image..." << std::endl;
  ConverterType::Pointer converter = ConverterType::New();
  converter->SetInput(img);
  converter->SetOutputMinimum(0);
  converter->SetOutputMaximum(255);
  converter->Update();
 
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("bias.gipl");
  writer->SetInput(converter->GetOutput());

  writer->Update();

  converter->SetInput(blur1->GetOutput());
  writer->SetFileName("p0.gipl");
  writer->Update();
  converter->SetInput(blur2->GetOutput());
  writer->SetFileName("p1.gipl");
  writer->Update();

  // Setup bias correction
  DynArray<FloatImageType::Pointer> probs;
  probs.Append(blur1->GetOutput());
  probs.Append(blur2->GetOutput());

  std::cout << "Create bias corrector..." << std::endl;
  BiasCorrectorType::Pointer biascorr = BiasCorrectorType::New();

  std::cout << "Set bias parameters..." << std::endl;
  biascorr->DebugOn();
  biascorr->SetClampBias(true);
  biascorr->SetSampleSpacing(2.0);
  biascorr->SetWorkingSpacing(2.0);
  biascorr->SetMask(maskImg);
  biascorr->SetMaxDegree(2);
  biascorr->SetProbabilities(probs);

  std::cout << "Bias correction... (max degree = ";
  std::cout << biascorr->GetMaxDegree() << ")" << std::endl;

  biascorr->Correct(img, corrImg, true);

  // Write image
  std::cout << "Write corrected image..." << std::endl;
  ConverterType::Pointer converter_corr = ConverterType::New();
  converter_corr->SetInput(corrImg);
  converter_corr->SetOutputMinimum(0);
  converter_corr->SetOutputMaximum(255);
  converter_corr->Update();
 
  writer->SetFileName("corr.gipl");
  writer->SetInput(converter_corr->GetOutput());

  writer->Update();

  return 0;

}
