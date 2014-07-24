
#include "EMSParametersXMLFile.h"
#include "NeoSegParametersXMLFile.h"

#include <iostream>

int
main()
{

  EMSParameters::Pointer p0 = EMSParameters::New();

  p0->SetSuffix("test0");
  p0->SetAtlasDirectory("/atlas");
  p0->SetOutputDirectory("/output");
  p0->SetOutputFormat("GIPL");
  p0->SetFilterIterations(10);
  p0->SetFilterTimeStep(0.05);
  p0->SetMaxBiasDegree(6);
  p0->SetPrior1(0.7);
  p0->SetPrior2(0.8);
  p0->SetPrior3(0.9);
  p0->SetPrior4(1.1);

  p0->AddImage("T1.gipl", "RAI");
  p0->AddImage("T2.mha", "ASR");

  writeEMSParametersXML("test.xml", p0);

  EMSParameters::Pointer p1 = readEMSParametersXML("test.xml");

  p1->PrintSelf(std::cout);

  NeoSegParameters::Pointer p2 = NeoSegParameters::New();

  p2->SetSuffix("test0");
  p2->SetAtlasDirectory("/atlasneo");
  p2->SetOutputDirectory("/outputneo");
  p2->SetOutputFormat("GIPL");
  p2->SetFilterIterations(10);
  p2->SetFilterTimeStep(0.05);
  p2->SetMaxBiasDegree(6);
  p2->SetPrior1(0.7);
  p2->SetPrior2(0.8);
  p2->SetPrior3(0.9);
  p2->SetPrior4(1.1);

  p2->AddImage("neoT1.gipl", "RPS");
  p2->AddImage("neoT2.mha", "AIR");

  p2->SetParzenKernelWidth(0.13);
  p2->SetPriorThreshold(0.79);
  p2->SetT2Index(2);

  writeNeoSegParametersXML("testneo.xml", p2);

  NeoSegParameters::Pointer p3 = readNeoSegParametersXML("testneo.xml");

  p3->PrintSelf(std::cout);

}
