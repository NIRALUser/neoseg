
//
// C++ class for the Mersenne Twister random number generator
// Converted from the original C implementation
//

// prastawa@cs.unc.edu 3/2005

// Original copyright text (from mt19937ar.c)
/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/em_StateVector.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#ifndef _MersenneTwisterRNG_h
#define _MersenneTwisterRNG_h

/* Period parameters */  
#define MT_N 624
#define MT_M 397

class MersenneTwisterRNG
{

public:
  MersenneTwisterRNG();
  ~MersenneTwisterRNG();

  void Initialize(unsigned long s);
  void InitializeByArray(unsigned long init_key[], int key_length);
  void InitializeUsingTime();

  // Get the pointer to the single global instance of the RNG class
  static MersenneTwisterRNG* GetGlobalInstance();

  unsigned long GenerateUniformInteger();

  // Integers in [0, k]
  unsigned long GenerateUniformIntegerUpToK(unsigned long k);

  unsigned int* GenerateIntegerSequence(unsigned int size, unsigned int max);

  // Generate values in [0, 1]
  double GenerateUniformRealClosedInterval();

  // Generate values in (0, 1)
  double GenerateUniformRealOpenInterval();

  // Generate values in [0, 1)
  double GenerateUniformRealClosedOpenInterval();

  double GenerateNormal(double mean, double variance);

protected:
  unsigned long m_StateVector[MT_N];
  int m_StateIterator;

};

#endif
