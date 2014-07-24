#include "MersenneTwisterRNG.h"

#include <stdio.h>

int main(void)
{
  MersenneTwisterRNG rng;

  int i;
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
  rng.InitializeByArray(init, length);

  printf("1000 outputs of genrand_int32()\n");
  for (i=0; i<1000; i++)
  {
    printf("%10lu ", rng.GenerateUniformInteger());
    if (i%5==4) printf("\n");
  }
  printf("\n1000 outputs of genrand_real2()\n");
  for (i=0; i<1000; i++)
  {
    printf("%10.8f ", rng.GenerateUniformRealClosedOpenInterval());
    if (i%5==4) printf("\n");
  }

  MersenneTwisterRNG seqRNG;
  unsigned int* seq = seqRNG.GenerateIntegerSequence(100, 5000);
  printf("\nSequence of 100 with max = 5000\n");
  for (unsigned int i = 0; i < 100; i++)
  {
    printf("%d ", seq[i]);
    if (i%5==4) printf("\n");
  } 

  return 0;
}
