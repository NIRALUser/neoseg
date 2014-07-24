
////////////////////////////////////////////////////////////////////////////////
//
// Handles log messages to terminal and disk, follows singleton pattern
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 2/2004

#ifndef _Log_h
#define _Log_h

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

class Log
{

public:

  static Log* GetInstance();

  void CloseFile();

  // Enable / disable writes to terminal
  inline void EchoOn() { m_EchoFlag = true; }
  inline void EchoOff() { m_EchoFlag = false; }

  void SetOutputFileName(const char* s);
  void SetOutputFileName(const std::string& s);

  void WriteString(const char* s);
  void WriteString(const std::string& s);

private:

  // Restrict access to constructors
  Log();
  ~Log();
  Log(const Log& l);

  bool m_EchoFlag;

  std::ofstream m_Output;

  std::string m_OutputFileName;

};

// Allows declarations such as: muLogMacro(<< "Message: " << 1.1234);
#define muLogMacro(x) \
  { \
    std::ostringstream outss; \
    outss << "" x << std::ends; \
    (Log::GetInstance())->WriteString(outss.str().c_str()); \
  }

#endif
