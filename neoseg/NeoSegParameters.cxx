
#include "NeoSegParameters.h"

NeoSegParameters
::NeoSegParameters()
{
  m_ReferenceImageIndex = 1;
  m_UseT1Flag = false;
  m_UseT2Flag = true;

  m_MahalanobisThreshold = 5.0;
  m_KernelWidth = 0.05;

  m_PriorThreshold = 0.9;

  m_Prior1 = 0.5;
  m_Prior5 = 1.0;

  m_AtlasFormat = "gipl" ; //In case the format is not specified, we use the old default format 
}

NeoSegParameters
::~NeoSegParameters()
{

}


bool
NeoSegParameters
::CheckValues()
{
  if (!Superclass::CheckValues())
    return false;

  if (m_ReferenceImageIndex < 1)
    return false;

  if (m_ReferenceImageIndex > m_Images.GetSize())
    return false;

  if (m_MahalanobisThreshold <= 0)
    return false;

  if (m_KernelWidth <= 0)
    return false;

  if (m_KernelWidth > 1)
    return false;

  if (m_PriorThreshold <= 0)
    return false;

  if (m_PriorThreshold >= 1.0)
    return false;

  return true;
}

void
NeoSegParameters
::PrintSelf(std::ostream& os)
{
  Superclass::PrintSelf(os);

  os << "Reference image index = " << m_ReferenceImageIndex << std::endl;
  os << "Use T1 = " << m_UseT1Flag << std::endl;
  os << "Use T2 = " << m_UseT2Flag << std::endl;
  os << "Mahalanobis threshold = " << m_MahalanobisThreshold << std::endl;
  os << "Kernel width = " << m_KernelWidth << std::endl;
  os << "Prior threshold = " << m_PriorThreshold << std::endl;
  os << "Atlas Format = " << m_AtlasFormat << std::endl ;
}
