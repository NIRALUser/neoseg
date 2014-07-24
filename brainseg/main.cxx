
#include "mu.h"

#include "itkOutputWindow.h"
#include "itkTextOutput.h"

#include "EMSParameters.h"
#include "EMSParametersXMLFile.h"
#include "runEMS.h"

#include <exception>
#include <iostream>


void
printUsage(char* progname)
{
  std::cerr << "Usage: " << progname << " <segfile> [options]" << std::endl;
  std::cerr << "Available options:" << std::endl;
  std::cerr << "--debug:\tdisplay debug messages" << std::endl;
  std::cerr << "--write-less:\twrite posteriors and filtered, bias corrected images";
  std::cerr << std::endl;
}

int
main(int argc, char** argv)
{

  std::cerr << "Run without any arguments to see command line options" << std::endl;

  if (argc < 2)
  {
    printUsage(argv[0]);
    return -1;
  }

  bool validargs = true;

  bool debugflag = false;
  bool writeflag = true;

  for (unsigned int i = 3; i < argc; i++)
  {
    if (strcmp(argv[i], "--debug") == 0)
      debugflag = true;
    else if (strcmp(argv[i], "--write-less") == 0)
      writeflag = false;
    else
      validargs = false;
  }

  if (!validargs)
  {
    printUsage(argv[0]);
    return -1;
  }

  itk::OutputWindow::SetInstance(itk::TextOutput::New());

  try
  {
    std::cout << "Reading " << argv[1] << "..." << std::endl;
    EMSParameters::Pointer emsp = readEMSParametersXML(argv[1]);
    runEMS(emsp, debugflag, writeflag);
  }
  catch (itk::ExceptionObject& e)
  {
    std::cerr << e << std::endl;
    return -1;
  }
  catch (std::exception& e)
  {
    std::cerr << "Exception: " << e.what() << std::endl;
    return -1;
  }
  catch (std::string& s)
  {
    std::cerr << "Exception: " << s << std::endl;
    return -1;
  }
  catch (...)
  {
    std::cerr << "Unknown exception" << std::endl;
    return -1;
  }

  return 0;

}
