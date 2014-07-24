
#include "muFile.h"

#include "NeoSegParametersXMLFile.h"

#include "itksys/SystemTools.hxx"

#include <fstream>
#include <sstream>

#include <stdlib.h>

NeoSegParametersXMLFileReader
::NeoSegParametersXMLFileReader()
{
  m_PObject = 0;
  m_LastFile = "";
  m_LastOrient = "";
}

NeoSegParametersXMLFileReader
::~NeoSegParametersXMLFileReader()
{

}

int
NeoSegParametersXMLFileReader
::CanReadFile(const char* name)
{
  if(!itksys::SystemTools::FileExists(name) ||
     itksys::SystemTools::FileIsDirectory(name) ||
     itksys::SystemTools::FileLength(name) == 0)
    return 0;
  return 1;
}

void
NeoSegParametersXMLFileReader
::StartElement(const char* name, const char** atts)
{
  if(itksys::SystemTools::Strucmp(name,"SEGMENTATION-PARAMETERS") == 0)
  {
    m_PObject = NeoSegParameters::New();
  }
  else if(itksys::SystemTools::Strucmp(name,"IMAGE") == 0)
  {
    m_LastFile = "";
    m_LastOrient = "";
  }
}

void
NeoSegParametersXMLFileReader
::EndElement(const char* name)
{
  if(itksys::SystemTools::Strucmp(name,"SEGMENTATION-PARAMETERS") == 0)
  {
    m_OutputObject = &(*m_PObject);
  }
  else if(itksys::SystemTools::Strucmp(name,"SUFFIX") == 0)
  {
    m_PObject->SetSuffix(m_CurrentString);
  }
  else if(itksys::SystemTools::Strucmp(name,"ATLAS-DIRECTORY") == 0)
  {
    m_PObject->SetAtlasDirectory(m_CurrentString);
  }
  else if(itksys::SystemTools::Strucmp(name,"ATLAS-ORIENTATION") == 0)
  {
    m_PObject->SetAtlasOrientation(m_CurrentString);
  }
  else if(itksys::SystemTools::Strucmp(name,"OUTPUT-DIRECTORY") == 0)
  {
    m_PObject->SetOutputDirectory(m_CurrentString);
  }
  else if(itksys::SystemTools::Strucmp(name,"ATLAS-FORMAT") == 0)
  {
    m_PObject->SetAtlasFormat(m_CurrentString);
  }
  else if(itksys::SystemTools::Strucmp(name,"OUTPUT-FORMAT") == 0)
  {
    m_PObject->SetOutputFormat(m_CurrentString);
  }
  else if(itksys::SystemTools::Strucmp(name,"FILE") == 0)
  {
    m_LastFile = m_CurrentString;
  }
  else if(itksys::SystemTools::Strucmp(name,"ORIENTATION") == 0)
  {
    m_LastOrient = m_CurrentString;
  }
  else if(itksys::SystemTools::Strucmp(name,"IMAGE") == 0)
  {
    m_PObject->AddImage(m_LastFile, m_LastOrient);
  }
  else if(itksys::SystemTools::Strucmp(name,"FILTER-ITERATIONS") == 0)
  {
    int iter = atoi(m_CurrentString.c_str());
    if (iter < 0)
      itkExceptionMacro(<< "Error: negative #iterations for filtering");
    m_PObject->SetFilterIterations(iter);
  }
  else if(itksys::SystemTools::Strucmp(name,"FILTER-TIME-STEP") == 0)
  {
    double dt = atof(m_CurrentString.c_str());
    m_PObject->SetFilterTimeStep(dt);
  }
  else if(itksys::SystemTools::Strucmp(name,"FILTER-METHOD") == 0)
  {
    m_PObject->SetFilterMethod(m_CurrentString);
  }
  else if(itksys::SystemTools::Strucmp(name,"MAX-BIAS-DEGREE") == 0)
  {
    int degree = atoi(m_CurrentString.c_str());
    if (degree < 0)
      itkExceptionMacro(<< "Error: negative bias degree");
    m_PObject->SetMaxBiasDegree(degree);
  }
  else if(itksys::SystemTools::Strucmp(name,"PRIOR-1") == 0)
  {
    double p = atof(m_CurrentString.c_str());
    m_PObject->SetPrior1(p);
  }
  else if(itksys::SystemTools::Strucmp(name,"PRIOR-2") == 0)
  {
    double p = atof(m_CurrentString.c_str());
    m_PObject->SetPrior2(p);
  }
  else if(itksys::SystemTools::Strucmp(name,"PRIOR-3") == 0)
  {
    double p = atof(m_CurrentString.c_str());
    m_PObject->SetPrior3(p);
  }
  else if(itksys::SystemTools::Strucmp(name,"PRIOR-4") == 0)
  {
    double p = atof(m_CurrentString.c_str());
    m_PObject->SetPrior4(p);
  }
  else if(itksys::SystemTools::Strucmp(name,"PRIOR-5") == 0)
  {
    double p = atof(m_CurrentString.c_str());
    m_PObject->SetPrior5(p);
  }
  else if(itksys::SystemTools::Strucmp(name,"DO-ATLAS-WARP") == 0)
  {
    unsigned int i = atoi(m_CurrentString.c_str());
    m_PObject->SetDoAtlasWarp(i != 0);
  }
  else if(itksys::SystemTools::Strucmp(name,"ATLAS-WARP-GRID-X") == 0)
  {
    unsigned int i = atoi(m_CurrentString.c_str());
    m_PObject->SetAtlasWarpGridX(i);
  }
  else if(itksys::SystemTools::Strucmp(name,"ATLAS-WARP-GRID-Y") == 0)
  {
    unsigned int i = atoi(m_CurrentString.c_str());
    m_PObject->SetAtlasWarpGridY(i);
  }
  else if(itksys::SystemTools::Strucmp(name,"ATLAS-WARP-GRID-Z") == 0)
  {
    unsigned int i = atoi(m_CurrentString.c_str());
    m_PObject->SetAtlasWarpGridZ(i);
  }
  else if(itksys::SystemTools::Strucmp(name,"MAHALANOBIS-THRESHOLD") == 0)
  {
    double t = atof(m_CurrentString.c_str());
    m_PObject->SetMahalanobisThreshold(t);
  }
  else if(itksys::SystemTools::Strucmp(name,"PARZEN-KERNEL-WIDTH") == 0)
  {
    double h = atof(m_CurrentString.c_str());
    m_PObject->SetKernelWidth(h);
  }
  else if(itksys::SystemTools::Strucmp(name,"PRIOR-THRESHOLD") == 0)
  {
    double p = atof(m_CurrentString.c_str());
    m_PObject->SetPriorThreshold(p);
  }
  else if(itksys::SystemTools::Strucmp(name,"REFERENCE-IMAGE-INDEX") == 0)
  {
    int i = atoi(m_CurrentString.c_str());
    if (i <= 0)
      itkExceptionMacro(<< "Error: index <= 0");
    m_PObject->SetReferenceImageIndex(i);
  }
  else if(itksys::SystemTools::Strucmp(name,"REFERENCE-MODALITY") == 0)
  {
    if (m_CurrentString.compare("T1") == 0)
      m_PObject->UseT1On();
    if (m_CurrentString.compare("T2") == 0)
      m_PObject->UseT2On();
  }
  // Support for older version of XML file
  else if(itksys::SystemTools::Strucmp(name,"T2-INDEX") == 0)
  {
    int i = atoi(m_CurrentString.c_str());
    if (i <= 0)
      itkExceptionMacro(<< "Error: index <= 0");
    m_PObject->SetReferenceImageIndex(i);
    m_PObject->UseT2On();
  }
}

void
NeoSegParametersXMLFileReader
::CharacterDataHandler(const char* inData, int inLength)
{
  m_CurrentString = "";
  for (int i = 0; i < inLength; i++)
    m_CurrentString += inData[i];
}

NeoSegParametersXMLFileWriter
::NeoSegParametersXMLFileWriter()
{

}

NeoSegParametersXMLFileWriter
::~NeoSegParametersXMLFileWriter()
{

}

int
NeoSegParametersXMLFileWriter
::CanWriteFile(const char* name)
{
  return true;
}

// Support function
template <typename T>
static void
WriteField(NeoSegParametersXMLFileWriter* writer, const char* attname,
  T value, std::ofstream& output)
{
  writer->WriteStartElement(attname, output);
  output << value;
  writer->WriteEndElement(attname, output);
  output << std::endl;
}

int
NeoSegParametersXMLFileWriter
::WriteFile()
{
  if (m_InputObject == 0)
    itkExceptionMacro(<< "No object to write");

  if (!m_InputObject->CheckValues())
    itkExceptionMacro(<< "Invalid values");

  if (m_Filename.length() == 0)
    itkExceptionMacro(<< "No file name specified");

  std::ofstream output(m_Filename.c_str());
  if (output.fail())
    itkExceptionMacro(<< "Can not open " << m_Filename);

  // Header
  WriteStartElement("?xml version=\"1.0\"?",output);
  output << std::endl;
  WriteStartElement("!DOCTYPE SEGMENTATION-PARAMETERS",output);
  output << std::endl;

  WriteStartElement("SEGMENTATION-PARAMETERS", output);
  output << std::endl;

  NeoSegParameters::Pointer p = m_InputObject;

  WriteField<std::string>(this, "SUFFIX", p->GetSuffix(), output);

  WriteField<std::string>(this, "ATLAS-DIRECTORY", p->GetAtlasDirectory(), output);

  WriteField<std::string>(this, "ATLAS-ORIENTATION", p->GetAtlasOrientation(), output);

  WriteField<std::string>(this, "OUTPUT-DIRECTORY", p->GetOutputDirectory(), output);

  WriteField<std::string>(this, "OUTPUT-FORMAT", p->GetOutputFormat(), output);

  // Write the list of images
  for (unsigned int k = 0; k < p->GetImages().GetSize(); k++)
  {
    this->WriteStartElement("IMAGE", output);

    output << std::endl;
    output << "  ";
    WriteField<std::string>(this, "FILE", p->GetImages()[k], output);
    output << "  ";
    WriteField<std::string>(this, "ORIENTATION", p->GetImageOrientations()[k],
      output);

    this->WriteEndElement("IMAGE", output);
    output << std::endl;
  }

  WriteField<unsigned int>(this, "FILTER-ITERATIONS", p->GetFilterIterations(), output);

  WriteField<float>(this, "FILTER-TIME-STEP", p->GetFilterTimeStep(), output);

  WriteField<std::string>(this, "FILTER-METHOD", p->GetFilterMethod(), output);

  WriteField<unsigned int>(this, "MAX-BIAS-DEGREE", p->GetMaxBiasDegree(), output);

  WriteField<float>(this, "PRIOR-1", p->GetPrior1(), output);
  WriteField<float>(this, "PRIOR-2", p->GetPrior2(), output);
  WriteField<float>(this, "PRIOR-3", p->GetPrior3(), output);
  WriteField<float>(this, "PRIOR-4", p->GetPrior4(), output);
  WriteField<float>(this, "PRIOR-5", p->GetPrior5(), output);

  WriteField<bool>(this, "DO-ATLAS-WARP", p->GetDoAtlasWarp(), output);

  WriteField<unsigned int>(this, "ATLAS-WARP-GRID-X", p->GetAtlasWarpGridX(), output);
  WriteField<unsigned int>(this, "ATLAS-WARP-GRID-Y", p->GetAtlasWarpGridY(), output);
  WriteField<unsigned int>(this, "ATLAS-WARP-GRID-Z", p->GetAtlasWarpGridZ(), output);

  WriteField<float>(this, "MAHALANOBIS-THRESHOLD", p->GetMahalanobisThreshold(), output);
  WriteField<float>(this, "PARZEN-KERNEL-WIDTH", p->GetKernelWidth(), output);
  WriteField<float>(this, "PRIOR-THRESHOLD", p->GetPriorThreshold(), output);

  WriteField<unsigned int>(this, "REFERENCE-IMAGE-INDEX", p->GetReferenceImageIndex(), output);
  if (p->UseT1())
    WriteField<std::string>(this, "REFERENCE-MODALITY", "T1", output);
  if (p->UseT2())
    WriteField<std::string>(this, "REFERENCE-MODALITY", "T2", output);

  // Finish
  WriteEndElement("SEGMENTATION-PARAMETERS", output);
  output << std::endl;
  output.close();

  return 0;
}

// Definition of some convenience functions
NeoSegParameters::Pointer
readNeoSegParametersXML(const char* fn)
{
  NeoSegParametersXMLFileReader::Pointer reader =
    NeoSegParametersXMLFileReader::New();
  try
  {
    reader->SetFilename(fn);
    reader->GenerateOutputInformation();
  }
  catch (...)
  {
    return 0;
  }
  return reader->GetOutputObject();
}

bool
writeNeoSegParametersXML(const char* fn, NeoSegParameters* p)
{
  if (p == 0)
    return false;
  if (!p->CheckValues())
    return false;

  // Enforce XML file extension
  std::string outfn = fn;
  std::string ext = mu::get_ext(fn);
  if (ext.compare("xml") != 0)
    outfn += std::string(".xml");

  NeoSegParametersXMLFileWriter::Pointer writer =
    NeoSegParametersXMLFileWriter::New();
  writer->SetFilename(outfn.c_str());
  writer->SetObject(p);
  writer->WriteFile();

  return true;
}

