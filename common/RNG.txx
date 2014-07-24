
#ifndef _RNG_txx
#define _RNG_txx

#include "RNG.h"

#include <iostream>

#include <stdlib.h>
#include <time.h>

template <class T>
RNG<T>
::RNG()
{

  m_RoundFlag = true;
  m_Minimum = 0;
  m_Maximum = 10;

  this->Initialize();

}

template <class T>
RNG<T>
::~RNG()
{

}

template <class T>
void
RNG<T>
::Initialize()
{

  // Let rand() take care of itself
  // Otherwise can have problems with multiple instances of RNG
  //srand(time(NULL));

}

template <class T>
void
RNG<T>
::SetMinimum(T i)
{

  if (i >= m_Maximum)
  {
    std::cerr << "[RNG::SetMinimum] Invalid range" << std::endl;
    return;
  }

  m_Minimum = i;

}

template <class T>
void
RNG<T>
::SetMaximum(T i)
{

  if (i <= m_Minimum)
  {
    std::cerr << "[RNG::SetMaximum] Invalid range" << std::endl;
    return;
  }

  m_Maximum = i;

}

template <class T>
double
RNG<T>
::GetRandomProbability()
{
  return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

template <class T>
T
RNG<T>
::GetValue()
{

  double r = this->GetRandomProbability();
  double x = r*(m_Maximum-m_Minimum) + m_Minimum;

  // Round x?
  if (m_RoundFlag)
    x = floor(x + 0.5);

  return static_cast<T>(x);

}

template <class T>
T*
RNG<T>
::GetRawArray(unsigned int size)
{

  if (size == 0)
    return NULL;

  unsigned int range = static_cast<unsigned int>(m_Maximum - m_Minimum + 1);

  // Should we generate array with unique elements (size less than range)
  bool unique = false;
  if (size <= range)
    unique = true;

  T* array = new T[size];

  if (unique)
  {
    // Mask for values already randomly selected
    unsigned char* selmask = new unsigned char[range];
    for (unsigned int i = 0; i < range; i++)
      selmask[i] = 0;

    for (unsigned int i = 0; i < size; i++)
    {
      T r;
      bool  insertOK;
      do
      {
        r = this->GetValue();
        unsigned int isel = static_cast<unsigned int>(r - m_Minimum);
        insertOK = (selmask[isel] == 0);
        selmask[isel] = 1;
      }
      while (!insertOK);
      array[i] = r;
    }

    delete [] selmask;
  }
  else
  {
    for (unsigned int i = 0; i < size; i++)
    {
      array[i] = this->GetValue();
    }
  }

  return array;

}

template <class T>
DynArray<T>
RNG<T>
::GetArray(unsigned int size)
{

  DynArray<T> array;

  T* rawarray = this->GetRawArray(size);

  // Return empty dynamic array if null
  if (rawarray == NULL)
    return array;

  array.Allocate(size);
  for (unsigned int i = 0; i < size; i++)
    array.Append(rawarray[i]);

  delete [] rawarray;

  return array;

}

#endif
