
////////////////////////////////////////////////////////////////////////////////
//
// Simple random number generator
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 8/2003

#ifndef _RNG_h
#define _RNG_h

#include "DynArray.h"

template <class T>
class RNG
{

public:
  RNG();
  ~RNG();

  void Initialize();

  // Set the range of values for the random numbers
  void SetMinimum(T i);
  void SetMaximum(T i);

  // Random real 0.0 - 1.0
  double GetRandomProbability();

  // Return a random number within the given range
  T GetValue();

  // Return an array of random numbers within the given range
  T* GetRawArray(unsigned int size);
  DynArray<T> GetArray(unsigned int size);

  inline bool RoundOn() { m_RoundFlag = true; }
  inline bool RoundOff() { m_RoundFlag = false; }

private:

  bool m_RoundFlag;
  T m_Minimum;
  T m_Maximum;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "RNG.txx"
#endif

#endif
