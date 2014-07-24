
#include "MRUTable.h"

#include <iostream>

class WrapInt
{
public:
  int x;
  WrapInt() { x = 0; }
  WrapInt(int i) { x = i; }
};

int main()
{

  MRUTable<WrapInt> mru(4);

  mru.Insert(WrapInt(1));
  mru.Insert(WrapInt(2));
  mru.Insert(WrapInt(3));
  mru.Insert(WrapInt(4));

  mru.GoToBegin();
  for (unsigned int i = 0; i < 4; i++)
    std::cout << mru.GetElement(i).x << " ";
  std::cout << std::endl;

  std::cout << "MRU size = " << mru.GetSize() << std::endl;

  mru.Insert(WrapInt(5));
  mru.Insert(WrapInt(6));
  mru.Insert(WrapInt(7));

  mru.GoToBegin();
  for (unsigned int i = 0; i < 4; i++)
    std::cout << mru.ReadNext().x << " ";
  std::cout << std::endl;

  std::cout << "MRU size = " << mru.GetSize() << std::endl;

  mru.Insert(WrapInt(8));
  mru.Insert(WrapInt(9));
  mru.Insert(WrapInt(10));

  mru.GoToBegin();
  for (unsigned int i = 0; i < 4; i++)
    std::cout << mru.GetElement(i).x << " ";
  std::cout << std::endl;

  std::cout << "MRU size = " << mru.GetSize() << std::endl;

  return 0;
}
