
#include "Heap.h"
#include "MSTClusteringProcess.h"

#include "muException.h"

#include "vnl/vnl_math.h"

#include <iostream>

// For sorting cluster maps based on size (descending)
class MSTCluster
{
public:
  unsigned int map;
  unsigned int size;

  MSTCluster() { this->map = 0; this->size = 0; }
  ~MSTCluster() { }

  inline MSTCluster& operator=(const MSTCluster& c)
  { this->map = c.map; this->size = c.size; return (*this); }

  inline bool operator<(const MSTCluster& c) const
  { return this->size > c.size; }
};


MSTClusteringProcess
::MSTClusteringProcess()
{
  m_NumberOfVertices = 0;

  m_MSTEdges = 0;

  m_NodeAverages = 0;

  m_SortFlag = false;
}

MSTClusteringProcess
::~MSTClusteringProcess()
{

  delete [] m_MSTEdges;
  delete [] m_NodeAverages;

}

void
MSTClusteringProcess
::SetInputVertices(const VertexList& l)
{

  m_NumberOfVertices = l.GetSize();

  delete [] m_MSTEdges;
  delete [] m_NodeAverages;

  if (m_NumberOfVertices == 0)
    return;

  m_MSTEdges = new MSTEdge[m_NumberOfVertices-1];
  m_NodeAverages = new float[m_NumberOfVertices];

  unsigned char* fringeMask = new unsigned char[m_NumberOfVertices];
  for (unsigned int k = 0; k < m_NumberOfVertices; k++)
    fringeMask[k] = 1;

  Heap<MSTEdge> vertexEdgeHeap;
  vertexEdgeHeap.Allocate(m_NumberOfVertices);

  // Insert edges connecting first vertex and the rest
  for (unsigned int k = 1; k < m_NumberOfVertices; k++)
  {
    MSTEdge e;
    e.i = k;
    e.j = 0;

    VertexType dvec = (l[0] - l[k]);
    float dist = dvec.squared_magnitude();

    e.dist = dist;
    vertexEdgeHeap.Insert(e);
  }

  unsigned int numEdges = 0;

  while (numEdges < (m_NumberOfVertices-1))
  {
    MSTEdge emin = vertexEdgeHeap.ExtractMinimum();

    // Insert vertex into tree
    fringeMask[emin.i] = 0;

    // Insert edge with minimum weight
    fringeMask[emin.j] = 0;
    m_MSTEdges[numEdges++] = emin;

    VertexType xmin = l[emin.i];

    // Go through heap and process fringe vertices
    MSTEdge* vertexEdgeArray = vertexEdgeHeap.GetElements();
    for (unsigned int k = 0; k < vertexEdgeHeap.GetNumberOfElements(); k++)
    {
      unsigned int i = vertexEdgeArray[k].i;
      if (fringeMask[i] == 0)
        continue;

      // Get distance
      VertexType dvec = (xmin - l[i]);
      float dist = dvec.squared_magnitude();

      if (dist < vertexEdgeArray[k].dist)
      {
        vertexEdgeArray[k].j = emin.i;
        vertexEdgeArray[k].dist = dist;
        vertexEdgeHeap.UpdateElementAt(k);
      }
    }
  } // while

  delete [] fringeMask;

#if 0
  // Debug
  for (int k = 0; k < numEdges; k++)
  {
    std::cout << "Edge " << k  << ": (" << m_MSTEdges[k].i << ", "
      << m_MSTEdges[k].j << ") " << m_MSTEdges[k].dist << std::endl;
  }
#endif

  // Compute node averages
  for (unsigned int k = 0; k < m_NumberOfVertices; k++)
    m_NodeAverages[k] = 0.0;

  unsigned int* countArray = new unsigned int[m_NumberOfVertices];
  for (unsigned int k = 0; k < m_NumberOfVertices; k++)
    countArray[k] = 0;

  for (unsigned int k = 0; k < (m_NumberOfVertices-1); k++)
  {
    unsigned int a = m_MSTEdges[k].i;
    unsigned int b = m_MSTEdges[k].j;

    m_NodeAverages[a] += m_MSTEdges[k].dist;
    countArray[a]++;

    m_NodeAverages[b] += m_MSTEdges[k].dist;
    countArray[b]++;
  }

  for (unsigned int k = 0; k < m_NumberOfVertices; k++)
  {
    if (countArray[k] != 0)
      m_NodeAverages[k] /= countArray[k];
  }

  delete [] countArray;

}

unsigned int
MSTClusteringProcess
::GetClusters(unsigned int* treeMap, float T)
{

  // Get number of vertices and edges
  unsigned int v = m_NumberOfVertices;
  unsigned int e = v - 1;

  // Allocate edge break flag array
  unsigned char* breakArray = new unsigned char[e];
  for (unsigned int k = 0; k < e; k++)
    breakArray[k] = 0;

  // Break edges
  unsigned int numBroken = 0;
  for (unsigned int i = 0; i < v; i++)
  {
    float thres = T * m_NodeAverages[i];

    // Never break zero length edges
    if (thres < 1e-10)
      thres = 1e-10;

    // Break the coinciding long edges
    for (unsigned int k = 0; k < e; k++)
    {
      if (breakArray[k] != 0)
        continue;

      unsigned int a = m_MSTEdges[k].i;
      unsigned int b = m_MSTEdges[k].j;

      bool incident = (i == a) || (i == b);

      if (incident && (m_MSTEdges[k].dist > thres))
      {
        breakArray[k] = 1;
        numBroken++;
      }
    }
  }

  if (numBroken == 0)
  {
    std::cerr << "No edges broken" << std::endl;
    delete [] breakArray;
    // All in one tree
    for (unsigned int k = 0; k < v; k++)
      treeMap[k] = 0;
    return 1;
  }

  // Figure out distinct trees, merge connected vertices
  for (unsigned int k = 0; k < v; k++)
    treeMap[k] = k;

  for (unsigned int i = 0; i < v; i++)
  {   

    unsigned int map1 = treeMap[i];
          
    // Check incident edges
    for (unsigned int j = 0; j < e; j++)
    {
      // Skip broken edges
      if (breakArray[j] != 0)
        continue;
    
      unsigned int a = m_MSTEdges[j].i;
      unsigned int b = m_MSTEdges[j].j; 
    
/*
      bool incident = (i == a) || (i == b);

      if (!incident)
        continue;
*/
        
      // Get the map of the other id and relabel
      unsigned int map2 = 0;
      if (i == a)
      {
        map2 = treeMap[b];
        for (unsigned int k = 0; k < v; k++)
          if (treeMap[k] == map2)
            treeMap[k] = map1;
      }
/*
      if (i == b)
      {
        map2 = treeMap[a];
        for (unsigned int k = 0; k < v; k++)
          if (treeMap[k] == map2)
            treeMap[k] = map1;
      }
*/

    } // for j

  } // for i

  delete [] breakArray;

  if (!m_SortFlag)
    return numBroken+1;

  //
  // Sort the cluster maps based on cluster size (descending)
  // Cluster 0 is the largest cluster
  //

  // Fill cluster info heap
  Heap<MSTCluster> heap;
  heap.Allocate(v);

  for (unsigned int i = 0; i < v; i++)
  {
    unsigned int s = 0;
    for (unsigned int j = 0; j < v; j++)
      if (treeMap[j] == i)
        s++;

    MSTCluster clust;
    clust.map = i;
    clust.size = s;

    if (s > 0)
      heap.Insert(clust);
  }

  // Sort based on size
  unsigned int* sortedMap = new unsigned int[v];
  unsigned int clustCount = 0;
  while (!heap.IsEmpty())
  {
    MSTCluster c = heap.ExtractMinimum();
    unsigned int m = c.map;

    for (unsigned int j = 0; j < v; j++)
      if (treeMap[j] == m)
        sortedMap[j] = clustCount;

    clustCount++;
  }

  for (unsigned int i = 0; i < v; i++)
    treeMap[i] = sortedMap[i];

  delete [] sortedMap;

  return numBroken+1;

}
