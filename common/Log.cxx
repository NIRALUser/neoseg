
#include "Log.h"

#include "muException.h"


Log*
Log
::GetInstance()
{

  // Allow only one instance
  static Log instance;

  return &instance;

}

Log
::Log()
{

  m_EchoFlag = true;
  m_OutputFileName = "";

}

Log
::~Log()
{

  if (m_Output.is_open())
    m_Output.close();

}

Log
::Log(const Log& l)
{

  m_EchoFlag = l.m_EchoFlag;
  m_OutputFileName = l.m_OutputFileName;

}

void
Log
::CloseFile()
{
  if (m_Output.is_open())
    m_Output.close();
}

void
Log
::SetOutputFileName(const char* s)
{
  if (m_Output.is_open())
    m_Output.close();

  m_Output.open(s);

  if (m_Output.fail())
  {
    muExceptionMacro(
      << "[Log::SetOutputFileName] Failed to open " << s);
  }
}

void
Log
::SetOutputFileName(const std::string& s)
{
  this->SetOutputFileName(s.c_str());
}

void
Log
::WriteString(const char* s)
{

  if (s == NULL)
  {
    std::cout << "[Log::WriteString] NULL argument" << std::endl;
    return;
  }

  if (m_Output.good())
  {
    m_Output << s;
    m_Output.flush();
  }

  if (m_EchoFlag)
  {
    std::cout << s;
    (std::cout).flush();
  }

}

void
Log
::WriteString(const std::string& s)
{

  this->WriteString(s.c_str());

}

// Test
#if 0

int main()
{

  Log* log = Log::GetInstance();

  log->SetOutputFileName("foo.log");

  log->WriteString("abc\n");
  log->WriteString(std::string("def\n"));

  muLogMacro(<< "foo " << 3 << "\n");

  return 0;

}

#endif
