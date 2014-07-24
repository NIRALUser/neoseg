
#ifndef _NeoSegParametersXMLFile_h
#define _NeoSegParametersXMLFile_h

#include "itkXMLFile.h"

#include "NeoSegParameters.h"

class NeoSegParametersXMLFileReader: public itk::XMLReader<NeoSegParameters>
{

public:

  // Standard typedefs
  typedef NeoSegParametersXMLFileReader Self;
  typedef itk::XMLReader<NeoSegParameters> Superclass;
  typedef itk::SmartPointer<Self> Pointer;

  // RTTI
  itkTypeMacro(Self, Superclass)

  itkNewMacro(Self);

  virtual int CanReadFile(const char* name);

protected:
  NeoSegParametersXMLFileReader();
  ~NeoSegParametersXMLFileReader();

  virtual void StartElement(const char * name, const char **atts);
  virtual void EndElement(const char *name);
  virtual void CharacterDataHandler(const char *inData, int inLength);

  NeoSegParameters::Pointer m_PObject;
  std::string m_CurrentString;

  std::string m_LastFile;
  std::string m_LastOrient;

};

class NeoSegParametersXMLFileWriter: public itk::XMLWriterBase<NeoSegParameters>
{

public:
  // Standard typedefs
  typedef NeoSegParametersXMLFileWriter Self;
  typedef itk::XMLWriterBase<NeoSegParameters> Superclass;
  typedef itk::SmartPointer<Self> Pointer;

  // RTTI
  itkTypeMacro(Self, Superclass)

  itkNewMacro(Self);

  virtual int CanWriteFile(const char* name);

  virtual int WriteFile();

protected:
  NeoSegParametersXMLFileWriter();
  ~NeoSegParametersXMLFileWriter();

};

// Convenience functions
NeoSegParameters::Pointer readNeoSegParametersXML(const char* fn);
bool writeNeoSegParametersXML(const char* fn, NeoSegParameters* p);

#endif
