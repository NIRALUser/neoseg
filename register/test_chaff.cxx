
#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include <iostream>
#include <fstream>

#include "ChainedAffineTransform3D.h"
//#include "PairRegistrationMethod.h"

int
main(int argc, char **argv)
{
  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  ChainedAffineTransform3D::Pointer trafo = ChainedAffineTransform3D::New();

  ChainedAffineTransform3D::ParametersType p(12);
  p[0] = 4;
  p[1] = 5;
  p[2] = 6;
  p[3] = 3.14 / 4.0;
  p[4] = -3.14 / 8.0;
  p[5] = 3.14 / 16.0;
  p[6] = 1.25;
  p[7] = 1.25;
  p[8] = 1.25;
  p[9] = 0.01;
  p[10] = 0.2;
  p[11] = 0.05;
  trafo->SetAllParameters(p, 64, 60, 70, 128, 100, 120, true);

  ChainedAffineTransform3D::Pointer invT = trafo->GetInverse();

  ChainedAffineTransform3D::InputPointType point;
  point[0] = 67.5;
  point[1] = 52.8;
  point[2] = 80.0;

  std::cout << "point = " << point << std::endl;

  point = trafo->TransformPoint(point);

  std::cout << "T(point) = " << point << std::endl;

  point = invT->TransformPoint(point);

  std::cout << "invT(T(point)) = " << point << std::endl;

  //PairRegistrationMethod<float>::WriteAffineTransform("test_chaff1.out", trafo);
  //PairRegistrationMethod<float>::WriteAffineTransform("test_chaff2.out", invT);

  std::cout << "------------" << std::endl;

  point[0] = -1.5;
  point[1] = 0.01;
  point[2] = 32.00209;

  std::cout << "point = " << point << std::endl;

  point = trafo->TransformPoint(point);

  std::cout << "T(point) = " << point << std::endl;

  point = invT->TransformPoint(point);

  std::cout << "invT(T(point)) = " << point << std::endl;

  std::cout << "------------" << std::endl;

  point[0] = 0.0;
  point[1] = 0.0;
  point[2] = 0.0;

  std::cout << "point = " << point << std::endl;

  point = trafo->TransformPoint(point);

  std::cout << "T(point) = " << point << std::endl;

  point = invT->TransformPoint(point);

  std::cout << "invT(T(point)) = " << point << std::endl;

  return 0;

}

