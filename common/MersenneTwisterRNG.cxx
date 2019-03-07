
#include "MersenneTwisterRNG.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>


// Constants
#define MT_MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define MT_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define MT_LOWER_MASK 0x7fffffffUL /* least significant r bits */

MersenneTwisterRNG
::MersenneTwisterRNG()
{

  m_StateIterator = MT_N+1; // Not initialized

}

MersenneTwisterRNG
::~MersenneTwisterRNG()
{

}

/* Initializes state vector with a seed */
void
MersenneTwisterRNG
::Initialize(unsigned long s)
{
  m_StateVector[0]= s & 0xffffffffUL;
  for (m_StateIterator=1; m_StateIterator<MT_N; m_StateIterator++)
  {
    m_StateVector[m_StateIterator] = 
     (1812433253UL *
       (m_StateVector[m_StateIterator-1] ^
       (m_StateVector[m_StateIterator-1] >> 30)) + m_StateIterator); 
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, MSBs of the seed affect   */
    /* only MSBs of the array m_StateVector[].                        */
    /* 2002/01/09 modified by Makoto Matsumoto             */
    m_StateVector[m_StateIterator] &= 0xffffffffUL;
        /* for >32 bit machines */
  }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void
MersenneTwisterRNG
::InitializeByArray(unsigned long init_key[], int key_length)
{
    this->Initialize(19650218UL);

    int i, j, k;

    i=1; j=0;
    k = ((MT_N > key_length) ? MT_N : key_length);
    for (; k; k--)
    {
      m_StateVector[i] =
        (m_StateVector[i] ^
          ((m_StateVector[i-1] ^ (m_StateVector[i-1] >> 30)) * 1664525UL))
        + init_key[j] + j; /* non linear */
      m_StateVector[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++; j++;
      if (i>=MT_N)
      {
        m_StateVector[0] = m_StateVector[MT_N-1];
        i=1;
      }
      if (j>=key_length)
        j=0;
    }
    for (k=MT_N-1; k; k--)
    {
      m_StateVector[i] =
        (m_StateVector[i] ^
          ((m_StateVector[i-1] ^ (m_StateVector[i-1] >> 30)) * 1566083941UL))
        - i; /* non linear */
      m_StateVector[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      if (i>=MT_N)
      {
        m_StateVector[0] = m_StateVector[MT_N-1];
        i=1;
      }
    }

    m_StateVector[0] = 0x80000000UL;
    /* MSB is 1; assuring non-zero initial array */ 
}

void
MersenneTwisterRNG
::InitializeUsingTime()
{
  unsigned long keys[2] = {(unsigned long)time(NULL), clock()};

  this->InitializeByArray(keys, 2);
}

MersenneTwisterRNG*
MersenneTwisterRNG
::GetGlobalInstance()
{
  static MersenneTwisterRNG globalRNG;

  return &globalRNG;
}

unsigned long 
MersenneTwisterRNG
::GenerateUniformInteger()
{
  unsigned long y;
  static unsigned long mag01[2]={0x0UL, MT_MATRIX_A};
  // mag01[x] = x * MT_MATRIX_A  for x=0,1

  // generate N words at one time
  if (m_StateIterator >= MT_N)
  {
    int kk;

    if (m_StateIterator == MT_N+1)   // if not initialized
      this->Initialize(5489UL); // a default initial seed is used

    for (kk=0;kk<MT_N-MT_M;kk++)
    {
      y = (m_StateVector[kk]&MT_UPPER_MASK)|(m_StateVector[kk+1]&MT_LOWER_MASK);
      m_StateVector[kk] = m_StateVector[kk+MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (;kk<MT_N-1;kk++)
    {
      y = (m_StateVector[kk]&MT_UPPER_MASK)|(m_StateVector[kk+1]&MT_LOWER_MASK);
      m_StateVector[kk] =
        m_StateVector[kk+(MT_M-MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (m_StateVector[MT_N-1]&MT_UPPER_MASK)|(m_StateVector[0]&MT_LOWER_MASK);
    m_StateVector[MT_N-1] = m_StateVector[MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    m_StateIterator = 0;
  }
  
  y = m_StateVector[m_StateIterator++];

  // Tempering
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

unsigned long 
MersenneTwisterRNG
::GenerateUniformIntegerUpToK(unsigned long k)
{
  // Find which bits are used in n
  // Optimized by Magnus Jonsson (magnus@smartelectronix.com)
  unsigned long used = k;
  used |= used >> 1;
  used |= used >> 2;
  used |= used >> 4;
  used |= used >> 8;
  used |= used >> 16;

  // Draw numbers until one is found in [0,k]
  unsigned long i;
  do
  {
    unsigned long u = this->GenerateUniformInteger();
    i = u & used;  // toss unused bits to shorten search
  }
  while (i > k);

  return i;
}

unsigned int*
MersenneTwisterRNG
::GenerateIntegerSequence(unsigned int size, unsigned int max)
{
  if (size == 0)
    return 0;

  // Build the array
  unsigned int* sequence = new unsigned int[size];

  unsigned int range = max+1; // 0 is included

  // Can we get unique elements? (size less than range)
  bool unique = false;
  if (size < range)
    unique = true;

  if (unique)
  {
    // Mask for values already selected
    unsigned char* selmask = new unsigned char[range];
    for (unsigned int i = 0; i < range; i++)
      selmask[i] = 0;

    for (unsigned int i = 0; i < size; i++)
    {
      unsigned int u = 0;
      bool insertOK;
      do
      {
        u = this->GenerateUniformIntegerUpToK(max);
        insertOK = (selmask[u] == 0);
        selmask[u] = 1;
      }
      while (!insertOK);
      sequence[i] = u;
    }

    delete [] selmask;
  }
  else
  {
    for (unsigned int i = 0; i < size; i++)
      sequence[i] = this->GenerateUniformIntegerUpToK(max);
  }

  return sequence;
}

// Note from original author:
// These real versions are due to Isaku Wada, 2002/01/09 added

double 
MersenneTwisterRNG
::GenerateUniformRealClosedInterval()
{
  return this->GenerateUniformInteger()*(1.0/4294967295.0); 
  /* divided by 2^32-1 */ 
}

double 
MersenneTwisterRNG
::GenerateUniformRealClosedOpenInterval()
{
#if 1
  // Uses 53-bit resolution
  unsigned long a = this->GenerateUniformInteger()>>5;
  unsigned long b = this->GenerateUniformInteger()>>6; 
  return(a*67108864.0+b)*(1.0/9007199254740992.0); 
#else
  return this->GenerateUniformInteger()*(1.0/4294967296.0); 
  /* divided by 2^32 */
#endif
}

double
MersenneTwisterRNG
::GenerateUniformRealOpenInterval()
{
  return (((double)this->GenerateUniformInteger()) + 0.5)*(1.0/4294967296.0); 
  /* divided by 2^32 */
}

double
MersenneTwisterRNG
::GenerateNormal(double mean, double variance)
{
  double u1 = this->GenerateUniformRealOpenInterval();
  double u2 = this->GenerateUniformRealClosedOpenInterval();
  double r = sqrt( -2.0 * log(1.0-u1)) * variance;
  double phi = 2.0 * 3.14159265358979323846264338328 * u2;
  return mean + r * cos(phi);
}
