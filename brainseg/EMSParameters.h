
#ifndef _EMSParameters_h
#define _EMSParameters_h

#include "itkObject.h"
#include "itkObjectFactory.h"

#include <iostream>
#include <string>
#include <algorithm>

#include "DynArray.h"

class EMSParameters: public itk::Object
{

public:

  typedef EMSParameters Self;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkNewMacro(Self);

  // Make sure all values are OK
  virtual bool CheckValues();

  virtual void PrintSelf(std::ostream& os);

  itkGetMacro(Suffix, std::string);
  itkSetMacro(Suffix, std::string);

  itkGetMacro(AtlasDirectory, std::string);
  itkSetMacro(AtlasDirectory, std::string);

  itkGetMacro(AtlasOrientation, std::string);
  itkSetMacro(AtlasOrientation, std::string);

  itkGetMacro(DoAtlasWarp, bool);
  itkSetMacro(DoAtlasWarp, bool);

  itkGetMacro(OutputDirectory, std::string);
  itkSetMacro(OutputDirectory, std::string);

  itkGetMacro(OutputFormat, std::string);
  itkSetMacro(OutputFormat, std::string);

  void AddImage(std::string s, std::string orientation);
  void ClearImages();

  DynArray<std::string> GetImages() { return m_Images; }
  DynArray<std::string> GetImageOrientations() { return m_ImageOrientations; }

  itkGetMacro(FilterMethod, std::string);
  itkSetMacro(FilterMethod, std::string);

  itkGetMacro(FilterIterations, unsigned int);
  itkSetMacro(FilterIterations, unsigned int);

  itkGetMacro(FilterTimeStep, float);
  itkSetMacro(FilterTimeStep, float);

  itkGetMacro(MaxBiasDegree, unsigned int);
  itkSetMacro(MaxBiasDegree, unsigned int);

  itkGetMacro(AtlasWarpGridX, unsigned int);
  itkSetMacro(AtlasWarpGridX, unsigned int);
  itkGetMacro(AtlasWarpGridY, unsigned int);
  itkSetMacro(AtlasWarpGridY, unsigned int);
  itkGetMacro(AtlasWarpGridZ, unsigned int);
  itkSetMacro(AtlasWarpGridZ, unsigned int);

  itkGetMacro(Prior1, float);
  itkSetMacro(Prior1, float);
  itkGetMacro(Prior2, float);
  itkSetMacro(Prior2, float);
  itkGetMacro(Prior3, float);
  itkSetMacro(Prior3, float);
  itkGetMacro(Prior4, float);
  itkSetMacro(Prior4, float);

protected:

  EMSParameters();
  ~EMSParameters();

  std::string m_Suffix;

  std::string m_AtlasDirectory;
  std::string m_AtlasOrientation;

  bool m_DoAtlasWarp;

  unsigned int m_AtlasWarpGridX;
  unsigned int m_AtlasWarpGridY;
  unsigned int m_AtlasWarpGridZ;

  std::string m_OutputDirectory;
  std::string m_OutputFormat;

  DynArray<std::string> m_Images;
  DynArray<std::string> m_ImageOrientations;

  std::string m_FilterMethod;
  unsigned int m_FilterIterations;
  float m_FilterTimeStep;

  unsigned int m_MaxBiasDegree;

  float m_Prior1;
  float m_Prior2;
  float m_Prior3;
  float m_Prior4;

};

#endif
