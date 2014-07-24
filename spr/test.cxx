
#include "DynArray.h"
#include "ReducedSetDensityEstimator.h"
#include "KernelWidthEstimator.h"
#include "KMeansEstimator.h"
#include "KNNClassifier.h"

#include "MersenneTwisterRNG.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>

int
main()
{

  srand(time(NULL));

  unsigned int n = 100;

  KNNClassifier::VectorType u1(2);
  KNNClassifier::VectorType u2(2);

  u1[0] = -2;
  u1[1] = -1;

  u2[0] = 2;
  u2[1] = 2;

  unsigned int nhalf = n / 2;

  // Create input set (two gaussians)
  KNNClassifier::MatrixType set(n, 2);
  DynArray<unsigned char> labels;
  labels.Initialize(n, 0);
  for (unsigned int i = 0; i < nhalf; i++)
  {
    float r1 = (float)rand() / (float)RAND_MAX - 0.5;
    float r2 = (float)rand() / (float)RAND_MAX - 0.5;

    KNNClassifier::VectorType s = u1;
    s[0] += r1;
    s[1] += r2;

    set.set_row(i, s);
    labels[i] = 1;
  }
  for (unsigned int i = nhalf; i < n; i++)
  {
    float r1 = (float)rand() / (float)RAND_MAX - 0.5;
    float r2 = (float)rand() / (float)RAND_MAX - 0.5;

    KNNClassifier::VectorType s = u2;
    s[0] += r1;
    s[1] += r2;

    set.set_row(i, s);
    labels[i] = 2;
  }

  KMeansEstimator kmeans;
  kmeans.SetNumberOfClusters(2);
  kmeans.SetNumberOfStarts(1000);
  kmeans.SetInput(set);

  KMeansEstimator::MatrixType kmu = kmeans.GetMeans();

  std::cout << "K-means: " << std::endl;
  std::cout << kmu.get_row(0) << " and " << kmu.get_row(1) << std::endl;

  // Create input set (one gaussians)
  KNNClassifier::MatrixType set2(n, 2);
  for (unsigned int i = 0; i < n; i++)
  {
    float r1 = (float)rand() / (float)RAND_MAX;
    float r2 = (float)rand() / (float)RAND_MAX;

    set2(i, 0) = r1;
    set2(i, 1) = r2;
  }

  // Create test set
  KNNClassifier::MatrixType testset(7, 2);

  testset(0, 0) = 0;
  testset(0, 1) = 0;

  testset.set_row(1, u1);
  testset.set_row(2, u2);

  testset(3, 0) = -5;
  testset(3, 1) = -5;

  testset(4, 0) = 5;
  testset(4, 1) = 5;

  testset(5, 0) = -0.5;
  testset(5, 1) = 3.0;

  testset(6, 0) = 0.5;
  testset(6, 1) = -1.0;

  unsigned int q = 10;
  KNNClassifier::MatrixType testset2(q, 2);
  for (unsigned int i = 0; i < q; i++)
  {
    float r1 = (float)rand() / (float)RAND_MAX;
    float r2 = (float)rand() / (float)RAND_MAX;

    testset2(i, 0) = r1;
    testset2(i, 1) = r2;
  }

  // KNN test
  KNNClassifier knn;
  knn.SetDimension(2);
  knn.SetKNeighbors(3);
  knn.SetTrainingData(set, labels);
  for (unsigned int i = 0; i < testset.rows(); i++)
  {
    KNNClassifier::VectorType x = testset.get_row(i);
    unsigned int c1 = knn.Classify(x);
    unsigned int c2 = knn.ClassifyWithoutCondense(x);
    std::cout << x << ":\t" << c1 << " || " << c2 << std::endl;
  }

  // RSDE test
  ReducedSetDensityEstimator rsde;
  rsde.SetDimension(2);
  rsde.SetInputSet(set2);
  rsde.SetKernelWidth(0.2);
  //float s = 1.0 / pow(2*M_PI, d/2.0); 
  //float s = 1.0 / (2*M_PI);
  for (unsigned int i = 0; i < testset2.rows(); i++)
  {
    KNNClassifier::VectorType x = testset2.get_row(i);
    float f1 = rsde.Evaluate(x);
    float f2 = rsde.EvaluateWithoutReduce(x);

    //float fideal = s*exp(-0.5 * x.squared_magnitude());
    float fideal = 0;
    if (x[0] >= 0 && x[0] <= 1 && x[1] >= 0 && x[1] <= 1)
      fideal = 1.0;

/*
    std::cout << x << ":\t" <<  f1 << "\t|| " << f2 << "\t|| fideal = " << fideal << std::endl;
    std::cout << "DIFF =\t" << fabs(f1-f2) << "\t|| ERR1 = " << fabs(f1-fideal) << "\t|| ERR2 = " << fabs(f2-fideal) << std::endl;
*/
    std::cout << "DIFF =\t" << fabs(f1-f2) << std::endl;
  }

  // Bandwidth estimator test
  KernelWidthEstimator hEst;

  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  KernelWidthEstimator::MatrixType hInput(50, 2);
  for (int i = 0; i < hInput.rows(); i++)
  {
    hInput(i, 0) = rng->GenerateNormal(0, 4.0);
    hInput(i, 1) = rng->GenerateNormal(0, 4.0);
  }
  std::cout << "Input for width est = \n" << hInput << std::endl;

  KernelWidthEstimator::MatrixType cov(2, 2);
  cov.set_identity();
  cov *= 4.0;

  KernelWidthEstimator::VectorType mu(2, 0.0);

  hEst.SetGaussianParameters(mu, cov);

  // Need to start from different locations?
  std::cout << "h = " << hEst.GetKernelWidth(hInput, 0.1) << std::endl;
  std::cout << "h = " << hEst.GetKernelWidth(hInput, 1.0) << std::endl;
  std::cout << "h = " << hEst.GetKernelWidth(hInput, 4.0) << std::endl;

}
