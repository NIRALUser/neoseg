
#include "Timer.h"

#include <math.h>

Timer
::Timer()
{
  this->Start();
}

void
Timer
::Start()
{
  m_StartTime = time(NULL);

  m_ElapsedHours = 0;
  m_ElapsedMinutes = 0;
  m_ElapsedSeconds = 0;

  m_Started = true;
}

void
Timer
::Stop()
{

  if (!m_Started)
    return;

  time_t stoptime = time(NULL);

  double secs = difftime(stoptime, m_StartTime);

  double hours = floor(secs / 3600.0);
  secs -= hours * 3600.0;

  double mins = floor(secs / 60.0);
  secs -= mins * 60.0;

  m_ElapsedHours = (unsigned int)hours;
  m_ElapsedMinutes = (unsigned int)mins;
  m_ElapsedSeconds = (unsigned int)secs;

  m_Started = false;

}
