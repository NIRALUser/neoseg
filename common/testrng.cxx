
#include "RNG.h"

#include <iostream>

int main()
{

  RNG<unsigned int> rng;

  rng.SetMinimum(2);
  rng.SetMaximum(15);

  DynArray<unsigned int> rlist = rng.GetArray(8);

  for (unsigned int i = 0; i < rlist.GetSize(); i++)
    std::cout << rlist[i] << std::endl;

  for (unsigned int i = 0; i < 10; i++)
    std::cout << rng.GetRandomProbability() << std::endl;

  return 0;

}
