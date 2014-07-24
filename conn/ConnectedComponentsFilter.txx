
#ifndef _ConnectedComponentsFilter_txx
#define _ConnectedComponentsFilter_txx

#include "itkNumericTraits.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"

#include "ConnectedComponentsFilter.h"

#include "Heap.h"

class ClusterInfo
{
public:
  ClusterInfo() { id = 0; order = 0; size = 0; }
  bool operator<(const ClusterInfo& c) const { return this->size > c.size; }
  const ClusterInfo& operator=(const ClusterInfo& c)
  { this->id = c.id; this->order = c.order; this->size = c.size; return *this; }
public:
  unsigned int id;
  unsigned int order;
  unsigned long size;
};

template< class TInputImage, class TOutputImage >
void
ConnectedComponentsFilter< TInputImage, TOutputImage >
::GenerateData()
{

  OutputImagePixelType outmax = itk::NumericTraits<OutputImagePixelType>::max();

  unsigned int maxNumberOfClusters = (unsigned int)(outmax - 1);

  itkDebugMacro(<< "maxNumberOfClusters = " << maxNumberOfClusters);

  unsigned int actualmax = (unsigned int)
    (itk::NumericTraits<unsigned int>::max() -
      itk::NumericTraits<unsigned int>::min());
  if (maxNumberOfClusters > actualmax)
    itkExceptionMacro(<< "Output pixel type has higher resolution than filter output");

  // Allocate space for the output image
  itkDebugMacro(<< "Allocate output");
  InputImageRegionType region = this->GetInput()->GetLargestPossibleRegion();

  this->GetOutput()->SetRegions(region);
  this->GetOutput()->Allocate();
  this->GetOutput()->SetOrigin(this->GetInput()->GetOrigin());
  this->GetOutput()->SetSpacing(this->GetInput()->GetSpacing());

  // Define iterators
  typedef itk::ImageRegionConstIterator<InputImageType> InputIterType;
  typedef itk::ImageRegionIterator<OutputImageType> OutputIterType;
  typedef itk::NeighborhoodIterator<OutputImageType> NeighborhoodIterType;

  InputIterType inputIt(this->GetInput(), region);
  OutputIterType outputIt(this->GetOutput(), region);

  typename NeighborhoodIterType::RadiusType radius;
  for (unsigned int i = 0; i < InputImageType::ImageDimension; i++)
    radius[i] = 1;
  NeighborhoodIterType nIt(radius, this->GetOutput(), region);

  // Set all output pixels to zero
  outputIt.GoToBegin();
  while (!outputIt.IsAtEnd())
  {
    outputIt.Set(0);
    ++outputIt;
  }

  ClusterInfo* clusters = new ClusterInfo[maxNumberOfClusters+1];

  // Initialize the list of clusters
  for (unsigned int i = 0; i <= maxNumberOfClusters; i++)
  {
    clusters[i].id = i;
    clusters[i].size = 0;
    clusters[i].order = 0;
  }

  unsigned int numberOfClusters = 0;

  // Find clusters
  itkDebugMacro(<< "Search clusters");
  inputIt.GoToBegin();
  nIt.GoToBegin();
  for ( ;!inputIt.IsAtEnd(); ++inputIt, ++nIt)
  {

    if (inputIt.Get() == 0)
      continue;

    unsigned int centerPix = (unsigned int)nIt.GetCenterPixel();
    unsigned int pixelID = clusters[centerPix].id;

    for (unsigned int i = 0; i < nIt.Size(); i++)
    {
      unsigned int neighborPix = (unsigned int)nIt.GetPixel(i);
      unsigned int neighborID = clusters[neighborPix].id;
      if (neighborID > 0)
      {
        if (pixelID > 0)
        {
          // Merge the two clusters by adjusting the ID references
          unsigned int useID = neighborID;
          unsigned int replaceID = pixelID;
          // Use the root reference
          if (clusters[pixelID].id == pixelID)
          {
            useID = pixelID;
            replaceID = neighborID;
          }
          pixelID = useID;
          if (useID != replaceID)
          {
            for (unsigned int i = 1; i <= numberOfClusters; i++)
              if (clusters[i].id == replaceID)
                clusters[i].id = useID;
          }
        }
        else
          pixelID = neighborID;
      }
    }

    if (pixelID == 0)
    {
      numberOfClusters++;
      if (numberOfClusters > maxNumberOfClusters)
      {
        numberOfClusters = maxNumberOfClusters;
      }
      pixelID = numberOfClusters;
    }

    nIt.SetCenterPixel(pixelID);
  }

  // Fix chained references
  itkDebugMacro(<< "Tracing cluster references");
  for (unsigned int i = 1; i <= numberOfClusters; i++)
  {
    if (clusters[i].id == i)
      continue;
    unsigned int j = clusters[i].id;
    while (clusters[j].id != j)
    {
      j = clusters[j].id;
    }
    clusters[i].id = j;
  }

  // Count number of root references
  unsigned int n = 0;
  for (unsigned int i = 1; i <= numberOfClusters; i++)
  {
    if (clusters[i].id == i)
      n++;
  }

  itkDebugMacro(<< "n = " << n);

  // Remap cluster labels
  itkDebugMacro(<< "Relabel 1");
  outputIt.GoToBegin();
  while (!outputIt.IsAtEnd())
  {
    if (outputIt.Get() > 0)
    {
      unsigned int theid = clusters[outputIt.Get()].id;
      outputIt.Set((OutputImagePixelType)theid);
    }
    ++outputIt;
  }

  // Compute the size of each clusters
  itkDebugMacro(<< "Compute cluster sizes");
  outputIt.GoToBegin();
  while (!outputIt.IsAtEnd())
  {
    if (outputIt.Get() > 0)
    {
      clusters[outputIt.Get()].size++;
    }
    ++outputIt;
  }

  // Insert clusters to heap and sort them by size (descending)
  Heap<ClusterInfo> cheap;
  for (unsigned int i = 0; i <= numberOfClusters; i++)
  {
    if (clusters[i].size > 0)
      cheap.Insert(clusters[i]);
  }

  // Set the cluster order (first occurence in sorted list)
  for (unsigned int i = 0; i <= numberOfClusters; i++)
  {
    if (cheap.IsEmpty())
      break;

    ClusterInfo cl = cheap.ExtractMinimum();

    unsigned int currID = cl.id;
    if (currID == 0)
      continue;

    // Only set the first occurence
    if (clusters[currID].order == 0)
      clusters[currID].order = i+1;
  }

  // Remap output labels, this time using the size ordering
  itkDebugMacro(<< "Relabel 2");
  outputIt.GoToBegin();
  while (!outputIt.IsAtEnd())
  {
    if (outputIt.Get() > 0)
    {
      unsigned int order = clusters[outputIt.Get()].order;
      OutputImagePixelType label = (OutputImagePixelType)order;
      outputIt.Set(label);
    }
    ++outputIt;
  }

  delete [] clusters;

}


#endif
