
#include "EMSParameters.h"

EMSParameters
::EMSParameters()
{
  m_Suffix = "";

  m_AtlasDirectory = "";

  m_AtlasOrientation = "RAI";

  m_DoAtlasWarp = true;

  m_OutputDirectory = "";
  m_OutputFormat = "Meta";

  m_Images.Clear();
  m_ImageOrientations.Clear();

  m_FilterMethod = "Curvature flow";
  m_FilterIterations = 1;
  m_FilterTimeStep = 0.01;

  m_MaxBiasDegree = 4;

  m_AtlasWarpGridX = 5;
  m_AtlasWarpGridY = 5;
  m_AtlasWarpGridZ = 5;

  m_Prior1 = 1.0;
  m_Prior2 = 1.0;
  m_Prior3 = 1.0;
  m_Prior4 = 1.0;
}

EMSParameters
::~EMSParameters()
{

}

void
EMSParameters
::AddImage(std::string s, std::string orient)
{
  m_Images.Append(s);
  m_ImageOrientations.Append(orient);
}

void
EMSParameters
::ClearImages()
{
  m_Images.Clear();
  m_ImageOrientations.Clear();
}

bool
EMSParameters
::CheckValues()
{
  if (m_Suffix.length() == 0)
    return false;

  if (m_AtlasDirectory.length() == 0)
    return false;

  if (m_OutputDirectory.length() == 0)
    return false;

  bool validFormat = false;

  m_OutputFormat.erase( m_OutputFormat.find_last_not_of( " \n\r\t") + 1 ) ;//if XML empty, format is just white spaces. We remove them.
  std::transform( m_OutputFormat.begin(), m_OutputFormat.end(), m_OutputFormat.begin(), ::tolower ) ;//set string to lower case to simplify the comparisons
  if( !m_OutputFormat.compare( "analyze" ) 
   || !m_OutputFormat.compare( "hdr" )
   || !m_OutputFormat.compare( "gipl" )
   || !m_OutputFormat.compare( "gipl.gz" )
   || !m_OutputFormat.compare( "nrrd" )
   || !m_OutputFormat.compare( "nhdr" )
   || !m_OutputFormat.compare( "nii.gz" )
   || !m_OutputFormat.compare( "nii" )
   || !m_OutputFormat.compare( "Meta" )
   || !m_OutputFormat.compare( "mha" )
   || !m_OutputFormat.compare( "mhd" )
         )
  {
    validFormat = true ;
  }



  if (!validFormat)
    return false;

  if (m_Images.GetSize() == 0)
    return false;

  return true;
}

void
EMSParameters
::PrintSelf(std::ostream& os)
{
  os << "Suffix = " << m_Suffix << std::endl;
  os << "Atlas directory = " << m_AtlasDirectory << std::endl;
  os << "Atlas orientation = " << m_AtlasOrientation << std::endl;
  os << "Output directory = " << m_OutputDirectory << std::endl;
  os << "Output format = " << m_OutputFormat << std::endl;
  os << "Images:" << std::endl;
  for (unsigned int k = 0; k < m_Images.GetSize(); k++)
    os << "  " << m_Images[k] << " --- " << m_ImageOrientations[k] << std::endl;
  os << "Filter iterations = " << m_FilterIterations << std::endl;
  os << "Filter time step = " << m_FilterTimeStep << std::endl;
  os << "Max bias degree = " << m_MaxBiasDegree << std::endl;
  os << "Prior 1 = " << m_Prior1 << std::endl;
  os << "Prior 2 = " << m_Prior2 << std::endl;
  os << "Prior 3 = " << m_Prior3 << std::endl;
  os << "Prior 4 = " << m_Prior4 << std::endl;
  if (m_DoAtlasWarp)
  {
    os << "Atlas warping, grid = "
       << m_AtlasWarpGridX << "x"
       << m_AtlasWarpGridY << "x"
       << m_AtlasWarpGridZ << std::endl;
  }
  else
  {
    os << "No atlas warping..." << std::endl;
  }
}
