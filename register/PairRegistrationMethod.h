
////////////////////////////////////////////////////////////////////////////////
//
// Registration of a pair of 3D images using different metrics and 
// transformations
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 5/2005

#ifndef _PairRegistrationMethod_h
#define _PairRegistrationMethod_h

#include "itkAffineTransform.h"
#include "itkBSplineDeformableTransform.h"
#include "itkImage.h"
#include "itkVector.h"

#include "ChainedAffineTransform3D.h"

#include <fstream>

template <class TPixel>
class PairRegistrationMethod
{

public:

  // typedefs
  typedef ChainedAffineTransform3D AffineTransformType;
  // 3-D spline, order 3
  typedef itk::BSplineDeformableTransform<double, 3, 3> BSplineTransformType;

  typedef itk::Image<TPixel, 3> ImageType;

  typedef itk::Vector<float, 3> VectorType;
  typedef itk::Image<VectorType, 3> DeformationFieldType;

  // Registration functions
  static AffineTransformType::Pointer
    RegisterAffine(ImageType* fixedImg, ImageType* movingImg);

  static BSplineTransformType::Pointer
    RegisterBSpline(ImageType* fixedImg, ImageType* movingImg,
      unsigned int nx, unsigned int ny, unsigned int nz);

  static DeformationFieldType::Pointer
    RegisterDemons(ImageType* fixedImg, ImageType* movingImg,
      unsigned int numIters=100);

  // Read / write functions
  static AffineTransformType::Pointer
    ReadAffineTransform(const char* fn);
  static void
    WriteAffineTransform(const char* fn, const AffineTransformType* affine);

  static BSplineTransformType::Pointer
    ReadBSplineTransform(const char* fn);
  static void
    WriteBSplineTransform(const char* fn, const BSplineTransformType* bspline);

  static DeformationFieldType::Pointer
    ReadDeformationField(const char* fn);
  static void
    WriteDeformationField(const char* fn, const DeformationFieldType* def);

private:

  // Find next uncommented line (doesn't begin with #)
  static void ReadNextLine(char* s, std::ifstream& infile);


};

#ifndef MU_MANUAL_INSTANTIATION
#include "PairRegistrationMethod.txx"
#endif

#endif
