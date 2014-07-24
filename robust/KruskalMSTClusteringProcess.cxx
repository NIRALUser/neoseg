
#include "Heap.h"
#include "KruskalMSTClusteringProcess.h"

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


KruskalMSTClusteringProcess
::KruskalMSTClusteringProcess()
{
  m_NumberOfVertices = 0;

  m_MSTEdges = 0;

  m_NodeAverages = 0;

  m_SortFlag = false;
}

KruskalMSTClusteringProcess
::~KruskalMSTClusteringProcess()
{

  delete [] m_MSTEdges;
  delete [] m_NodeAverages;

}

void
KruskalMSTClusteringProcess
::SetInputVertices(const VertexList& l)
{

  m_NumberOfVertices = l.GetSize();

  delete [] m_MSTEdges;
  delete [] m_NodeAverages;

  if (m_NumberOfVertices == 0)
    return;

  m_MSTEdges = new MSTEdge[m_NumberOfVertices-1];
  m_NodeAverages = new float[m_NumberOfVertices];

  unsigned int numEdges = m_NumberOfVertices * (m_NumberOfVertices-1) / 2;

  // Create edges connecting the possible pair of vertices, insert to heap
  Heap<MSTEdge> heap;
  heap.Allocate(numEdges);
  for (unsigned int i = 0; i < m_NumberOfVertices; i++)
  {
    VertexType x = l[i];
    for (unsigned int j = i+1; j < m_NumberOfVertices; j++)
    {
      VertexType dij = (x - l[j]);

      MSTEdge e;
      e.i = i;
      e.j = j;
      e.dist = dij.magnitude();

      heap.Insert(e);
    }
  }

  // Map vertex to set (tree) it belongs to, for cycle test
  unsigned int* treeMap = new unsigned int[m_NumberOfVertices];
  for (unsigned int i = 0; i < m_NumberOfVertices; i++)
    treeMap[i] = i;

  // Number of edges in MST
  unsigned int edgeCount = 0;

  // Build MST using Kruskal's algorithm
  //
  // Edges added in ascending order
  while (!heap.IsEmpty())
  {

    MSTEdge minEdge = heap.ExtractMinimum();

    unsigned int a = minEdge.i;
    unsigned int b = minEdge.j;

    unsigned int map1 = treeMap[a];
    unsigned int map2 = treeMap[b];

    // Skip if they belong to the same tree (will form cycle)
    if (map1 == map2)
      continue;

    // Merge trees
    for (unsigned int k = 0; k < m_NumberOfVertices; k++)
    {
      if (treeMap[k] == map2)
        treeMap[k] = map1;
    }

    m_MSTEdges[edgeCount] = minEdge;
    
    edgeCount++; 

    // See if a tree is formed already
    if (edgeCount == (m_NumberOfVertices-1))
      break;
    
  } 
      
  delete [] treeMap;

  if (edgeCount != (m_NumberOfVertices-1))
  {
    std::cerr << "MST construction failed, E != (V-1)" << std::endl;
    throw "MST construction failed, E != (V-1)";
  }

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
KruskalMSTClusteringProcess
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

    // Break the coinciding long edges
    for (unsigned int k = 0; k < e; k++)
    {
      if (breakArray[k] != 0)
        continue;

      unsigned int a = m_MSTEdges[k].i;
      unsigned int b = m_MSTEdges[k].j;

      bool incident = (i == a) || (i == b);

      // Never break zero length edges
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
// TODO: FIXME:
//return whole tree with same label
    return 0;
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
      if (breakArray[j] != 0)
        continue;
    
      unsigned int a = m_MSTEdges[j].i;
      unsigned int b = m_MSTEdges[j].j; 
    
      bool incident = (i == a) || (i == b);

      if (!incident)
        continue;
        
      // Get the map of the other id
      unsigned int map2 = treeMap[b];
      if (i == b)
        map2 = treeMap[a];

      if (map1 == map2)
        continue;

      for (unsigned int k = 0; k < v; k++)
        if (treeMap[k] == map2)
          treeMap[k] = map1;
    } // for neighbors

  } // for i

  delete [] breakArray;

  if (!m_SortFlag)
    return numBroken+1;

  //
  // Sort the cluster maps based on cluster size (descending)
  // Cluster 0 is the largest cluster
  //

  // Obtain cluster info
  MSTCluster* clusters = new MSTCluster[v];
  unsigned int numNonZero = 0;
  for (unsigned int i = 0; i < v; i++)
  {
    clusters[i].map = i;
    unsigned int s = 0;
    for (unsigned int j = 0; j < v; j++)
      if (treeMap[j] == i)
        s++;
    if (s > 0)
      numNonZero++;
    clusters[i].size = s;
  }

  Heap<MSTCluster> heap;
  heap.Allocate(numNonZero);

  for (unsigned int i = 0; i < v; i++)
    if (clusters[i].size != 0)
      heap.Insert(clusters[i]);

  delete [] clusters;

  unsigned int* sortedMap = new unsigned int[v];

  for (unsigned int i = 0; i < numNonZero; i++)
  {
    MSTCluster c = heap.ExtractMinimum();
    unsigned int m = c.map;

    for (unsigned int j = 0; j < v; j++)
      if (treeMap[j] == m)
        sortedMap[j] = i;
  }

  for (unsigned int i = 0; i < v; i++)
    treeMap[i] = sortedMap[i];

  delete [] sortedMap;

  return numBroken+1;

}
