
#ifndef _NeoSegParameters_h
#define _NeoSegParameters_h

#include "EMSParameters.h"

class NeoSegParameters: public EMSParameters
{

public:

  typedef NeoSegParameters Self;
  typedef EMSParameters Superclass;

  typedef itk::SmartPointer<Self> Pointer;

  itkNewMacro(Self);

  virtual bool CheckValues();
  virtual void PrintSelf(std::ostream& os);

  void UseT1On() { m_UseT1Flag = true; m_UseT2Flag = false; }
  void UseT2On() { m_UseT1Flag = false; m_UseT2Flag = true; }

  bool UseT1() { return m_UseT1Flag; }
  bool UseT2() { return m_UseT2Flag; }

  itkGetMacro(AtlasFormat, std::string);
  itkSetMacro(AtlasFormat, std::string);

  itkGetMacro(ReferenceImageIndex, unsigned int);
  itkSetMacro(ReferenceImageIndex, unsigned int);

  itkGetMacro(MahalanobisThreshold, float);
  itkSetMacro(MahalanobisThreshold, float);

  itkGetMacro(PriorThreshold, float);
  itkSetMacro(PriorThreshold, float);

  itkGetMacro(KernelWidth, float);
  itkSetMacro(KernelWidth, float);

  itkGetMacro(Prior5, float);
  itkSetMacro(Prior5, float);

protected:

  NeoSegParameters();
  ~NeoSegParameters();

private:

  bool m_UseT1Flag;
  bool m_UseT2Flag;

  std::string m_AtlasFormat;

  unsigned int m_ReferenceImageIndex;

  float m_MahalanobisThreshold;

  float m_PriorThreshold;

  float m_KernelWidth;

  float m_Prior5;

};

#endif
