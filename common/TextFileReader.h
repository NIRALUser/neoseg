
// Converts a text file to a list of strings, each string representing
// a line in the file

#ifndef _TextFileReader_h
#define _TextFileReader_h

#include <string>

#include "DynArray.h"

class TextFileReader
{

public:
  TextFileReader();
  ~TextFileReader();

  inline std::string GetFileName() { return m_FileName; }
  void SetFileName(const char* s);

  std::string GetNextLine();
  std::string GetLine(unsigned int i);

  bool Remaining();

  inline unsigned int GetCurrentLineNumber() { return m_CurrentLineNumber; }

private:

  std::string m_FileName;

  unsigned int m_CurrentLineNumber;

  DynArray<std::string> m_Lines;

};

#endif
