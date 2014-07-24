
#include "FastMCDSampleFilter.h"
#include "MSTClusteringProcess.h"

#include <iostream>

int main()
{

  //
  // MCD test case
  //

  // 1-D test
  FastMCDSampleFilter::MatrixType input1d(12, 1);
  input1d(0, 0) = 3;
  input1d(1, 0) = 1;
  input1d(2, 0) = 2;
  input1d(3, 0) = 4.9;
  input1d(4, 0) = 5;
  input1d(5, 0) = 6;
  input1d(6, 0) = 4;
  input1d(7, 0) = 5000;
  input1d(8, 0) = 1.7;
  input1d(9, 0) = 5;
  input1d(10, 0) = 3.5;
  input1d(11, 0) = 1200;

  FastMCDSampleFilter::MatrixType input2d(12, 2);
  // First dimension
  input2d(0, 0) = 3;
  input2d(1, 0) = 1;
  input2d(2, 0) = 2;
  input2d(3, 0) = 4.9;
  input2d(4, 0) = 5;
  input2d(5, 0) = 6;
  input2d(6, 0) = 4;
  input2d(7, 0) = 5000;
  input2d(8, 0) = 1.7;
  input2d(9, 0) = 5;
  input2d(10, 0) = 3.5;
  input2d(11, 0) = 1200;
  // Second dimension
  input2d(0, 1) = 7;
  input2d(1, 1) = 4;
  input2d(2, 1) = 9000;
  input2d(3, 1) = 3;
  input2d(4, 1) = 5.1;
  input2d(5, 1) = 5;
  input2d(6, 1) = 7;
  input2d(7, 1) = 30;
  input2d(8, 1) = 6.1;
  input2d(9, 1) = 300;
  input2d(10, 1) = 8;
  input2d(11, 1) = 15;

  std::cout << "Before" << std::endl;
  std::cout << "------" << std::endl;
  std::cout << "Input for 1-D: \n" << input1d << std::endl;
  std::cout << "Input for 2-D: \n" << input2d << std::endl;

  FastMCDSampleFilter mcdfilt;

  FastMCDSampleFilter::MatrixType mu;
  FastMCDSampleFilter::MatrixType cov;

  mcdfilt.GetRobustEstimate(mu, cov, input1d);

  std::cout << "Robust estimate for 1-D: " << std::endl;
  std::cout << "Mean = " << mu;
  std::cout << "Covariance = \n" << cov << std::endl;

  mcdfilt.GetRobustEstimate(mu, cov, input2d);

  std::cout << "Robust estimate for 2-D: " << std::endl;
  std::cout << "Mean = " << mu;
  std::cout << "Covariance = \n" << cov << std::endl;

  mcdfilt.SetMahalanobisThreshold(4.0);
  FastMCDSampleFilter::MatrixType output1d = mcdfilt.GetInliers(input1d);
  FastMCDSampleFilter::MatrixType output2d = mcdfilt.GetInliers(input2d);

  std::cout << "Input points > 4.0 stddev removed" << std::endl;
  std::cout << "------" << std::endl;
  std::cout << "1-D case:\n" <<  output1d << std::endl;
  std::cout << "2-D case:\n" <<  output2d << std::endl;

  // Test case for MST clustering
  DynArray<MSTClusteringProcess::VertexType> mstsamples;

  MSTClusteringProcess::VertexType x(2);

  x[0] = 1;
  x[1] = 1;
  mstsamples.Append(x);

  x[0] = 2;
  x[1] = 5;
  mstsamples.Append(x);

  x[0] = 3;
  x[1] = 2;
  mstsamples.Append(x);

  x[0] = 30;
  x[1] = 21;
  mstsamples.Append(x);

  x[0] = 1;
  x[1] = 10000;
  mstsamples.Append(x);

  x[0] = 25;
  x[1] = 17;
  mstsamples.Append(x);

  x[0] = 27;
  x[1] = 18;
  mstsamples.Append(x);

  MSTClusteringProcess mstProc;
  mstProc.SetInputVertices(mstsamples);
  mstProc.SortOn();

  unsigned int* cmap = new unsigned int[mstsamples.GetSize()];

  for (double T = 2.0; T >= 1.0; T -= 0.2)
  {
    std::cout << "Cluster labels for T = " << T << ":" << std::endl;

    unsigned int numC = mstProc.GetClusters(cmap, T);
    std::cout << "Obtained " << numC << " separate clusters" << std::endl;

    if (numC == 0)
      continue;

    for (unsigned int i = 0; i < mstsamples.GetSize(); i++)
    {
      std::cout << mstsamples[i][0] << ", " << mstsamples[i][1] << " -> ";
      std::cout <<  cmap[i] << std::endl;
    }
    std::cout << std::endl;
  }

  return 0;
}
