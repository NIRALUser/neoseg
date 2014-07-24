
#include "Pair.h"

#include <iostream>

int main()
{

  typedef Pair<int> PairType;

  PairType p(1, 2);

  std::cout << p.GetFirst() << ", " << p.GetSecond() << std::endl;

}
