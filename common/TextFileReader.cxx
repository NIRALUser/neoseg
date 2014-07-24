
#include "TextFileReader.h"

#include <fstream>

#define BUFFER_SIZE 16384

TextFileReader
::TextFileReader()
{

  m_FileName = "";
  m_CurrentLineNumber = 0;

}

TextFileReader
::~TextFileReader()
{

}

std::string
TextFileReader
::GetNextLine()
{
  if (m_CurrentLineNumber >= m_Lines.GetSize())
    return "";

  ++m_CurrentLineNumber;

  return m_Lines[m_CurrentLineNumber-1];
}

std::string
TextFileReader
::GetLine(unsigned int i)
{
  if (i >= m_Lines.GetSize())
    return "";

  return m_Lines[i];
}

bool
TextFileReader
::Remaining()
{
  if (m_CurrentLineNumber >= m_Lines.GetSize())
    return false;

  return true;
}

void
TextFileReader
::SetFileName(const char* s)
{
  // Reset
  m_CurrentLineNumber = 0;
  m_Lines.Clear();

  if (s == NULL)
    return;

  std::ifstream fin;

  fin.open(s);
  if (fin.fail())
    return;

  char* buffer = new char[BUFFER_SIZE];
  while (!fin.eof())
  {
    fin.getline(buffer, BUFFER_SIZE); 
    m_Lines.Append(std::string(buffer));
  }
  delete [] buffer;

  fin.close();
}
