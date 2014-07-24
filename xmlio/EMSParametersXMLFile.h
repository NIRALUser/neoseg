
#ifndef _EMSParametersXMLFile_h
#define _EMSParametersXMLFile_h

#include "itkXMLFile.h"

#include "EMSParameters.h"

class EMSParametersXMLFileReader: public itk::XMLReader<EMSParameters>
{

public:

  // Standard typedefs
  typedef EMSParametersXMLFileReader Self;
  typedef itk::XMLReader<EMSParameters> Superclass;
  typedef itk::SmartPointer<Self> Pointer;

  // RTTI
  itkTypeMacro(Self, Superclass)

  itkNewMacro(Self);

  virtual int CanReadFile(const char* name);

protected:
  EMSParametersXMLFileReader();
  ~EMSParametersXMLFileReader();

  virtual void StartElement(const char * name, const char **atts);
  virtual void EndElement(const char *name);
  virtual void CharacterDataHandler(const char *inData, int inLength);

  EMSParameters::Pointer m_PObject;
  std::string m_CurrentString;

  std::string m_LastFile;
  std::string m_LastOrient;

};

class EMSParametersXMLFileWriter: public itk::XMLWriterBase<EMSParameters>
{

public:
  // Standard typedefs
  typedef EMSParametersXMLFileWriter Self;
  typedef itk::XMLWriterBase<EMSParameters> Superclass;
  typedef itk::SmartPointer<Self> Pointer;

  // RTTI
  itkTypeMacro(Self, Superclass)

  itkNewMacro(Self);

  virtual int CanWriteFile(const char* name);

  virtual int WriteFile();

protected:
  EMSParametersXMLFileWriter();
  ~EMSParametersXMLFileWriter();

};

// Convenience functions
EMSParameters::Pointer readEMSParametersXML(const char* fn);
bool writeEMSParametersXML(const char* fn, EMSParameters* p);

#endif
