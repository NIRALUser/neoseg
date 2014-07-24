
#ifndef _muMacros_h
#define _muMacros_h

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <sstream>
#include <string>

#define muEcho(varname) \
  std::cout << #varname << " = " << varname << std::endl;

#define muStringMacro(strname, s) \
  std::string strname; \
  { \
    std::ostringstream outss; \
    outss << "" s << std::ends; \
    strname = outss.str(); \
  }

#define muSelfFilterMacro(filter, obj) \
  { \
    filter->SetInput(obj); \
    iterator copy \
  }

#define muReadMacro(type, filename, image) \
  { \
    typedef itk::ImageFileReader<type> ReaderType; \
    typename ReaderType::Pointer reader = ReaderType::New(); \
    reader->SetFileName(filename); \
    reader->Update();\
    image = reader->GetOutput(); \
  }

#define muWriteMacro(type, filename, image) \
  { \
    typedef itk::ImageFileWriter<type> WriterType; \
    typename WriterType::Pointer writer = WriterType::New(); \
    writer->SetFileName(filename); \
    writer->SetInput(image); \
    writer->Update();\
  }

#endif
