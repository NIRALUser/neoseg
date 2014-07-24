
////////////////////////////////////////////////////////////////////////////////
//
// Handles simple timings
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 10/2003

#ifndef _Timer_h
#define _Timer_h

#include <time.h>

class Timer
{

public:

  Timer();
  ~Timer() {}

  void Start();
  void Stop();

  inline unsigned int GetElapsedHours()
  { this->Stop(); return m_ElapsedHours; }
  inline unsigned int GetElapsedMinutes()
  { this->Stop(); return m_ElapsedMinutes; }
  inline unsigned int GetElapsedSeconds()
  { this->Stop(); return m_ElapsedSeconds; }

private:

  time_t m_StartTime;

  unsigned int m_ElapsedHours;
  unsigned int m_ElapsedMinutes;
  unsigned int m_ElapsedSeconds;

  bool m_Started;

};

#endif
