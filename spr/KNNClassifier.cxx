
#include "KNNClassifier.h"

#include "DynArray.h"
#include "Heap.h"
#include "Log.h"
#include "MersenneTwisterRNG.h"

#include "muException.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>

class TaggedDistance
{
public:
  float dist;
  unsigned char label;

  TaggedDistance() { this->label = 0; this->dist = 0; }
  ~TaggedDistance() { }

  inline TaggedDistance& operator=(const TaggedDistance& d)
  { this->dist = d.dist; this->label = d.label; return (*this); }

  inline bool operator<(const TaggedDistance& d) const
  { return this->dist < d.dist; }
};

KNNClassifier
::KNNClassifier()
{

  m_Dimension = 1;
  m_KNeighbors = 1;

  m_Prototypes = MatrixType(0, m_Dimension);
  m_CondensedPrototypes = MatrixType(0, m_Dimension);

  m_LabelRange = 0;

  m_Modified = false;

}

KNNClassifier
::~KNNClassifier()
{

}

void
KNNClassifier
::SetDimension(unsigned int d)
{
  if (d == m_Dimension)
    return;

  m_Dimension = d;
  m_Modified = true;
}

void
KNNClassifier
::SetKNeighbors(unsigned int k)
{
  if (k == m_KNeighbors)
    return;
  m_KNeighbors = k;
  m_Modified = true;
}

void
KNNClassifier
::SetTrainingData(MatrixType& protos, LabelArray& labels)
{
  if (protos.columns() != m_Dimension)
    muExceptionMacro(
      << "[KNNClassifier::SetTrainingData] Dimension mismatch");

  if (protos.rows() != labels.GetSize())
    muExceptionMacro(
      << "[KNNClassifier::SetTrainingData] Labels and prototypes size mismatch");

  m_Prototypes = protos;
  m_Labels = labels;

  unsigned int maxL = labels[0];
  for (unsigned int i = 1; i < m_Labels.GetSize(); i++)
  {
    m_Labels[i] = labels[i];
    if (labels[i] > maxL)
      maxL = labels[i];
  }

  m_LabelRange = maxL + 1;

  m_Modified = true;
}

void
KNNClassifier
::Update()
{

  unsigned int n = m_Prototypes.rows();

  unsigned char* marks = new unsigned char[n];
  for (unsigned int i = 0; i < n; i++)
    marks[i] = 1;

  // Condense prototypes by eliminating Gabriel neighbors with same labels
  for (unsigned int i = 0; i < n; i++)
  {
    VectorType a = m_Prototypes.get_row(i);

    // Find labels of Gabriel neighbors
    DynArray<unsigned char> neighborLabels;
    for (unsigned int j = 0; j < n; j++)
    {
      if (j == i)
        continue;

      VectorType a_to_center = (m_Prototypes.get_row(j) - a) * 0.5;
      VectorType center = a + a_to_center;

      float r_squared = a_to_center.squared_magnitude();

      // Brute-force search, make sure no other point in hypersphere
      bool isEdge = true;
      for (unsigned int k = 0; k < n; k++)
      {
        if (k == i || k == j)
          continue;
        float dist_squared =
          (m_Prototypes.get_row(k) - center).squared_magnitude();
        if (dist_squared < r_squared)
        {
          isEdge = false;
          break;
        }
      } // for k

      if (isEdge)
        neighborLabels.Append(m_Labels[j]);

    } // for j

    unsigned char label = m_Labels[i];

    unsigned int countMatch = 0;
    for (unsigned int q = 0; q < neighborLabels.GetSize(); q++)
      if (neighborLabels[q] == label)
        countMatch++;

    if ((countMatch >= m_KNeighbors)
        &&
        (countMatch == neighborLabels.GetSize()))
      marks[i] = 0;

  } // for i

  unsigned int m = 0;
  for (unsigned int i = 0; i < n; i++)
    if (marks[i] != 0)
      m++;

  m_CondensedPrototypes = MatrixType(m, m_Dimension);

  m_CondensedLabels.Initialize(m, 0);

  unsigned int r = 0;
  for (unsigned int i = 0; i < n; i++)
  {
    if (marks[i] == 0)
      continue;
    m_CondensedPrototypes.set_row(r, m_Prototypes.get_row(i));
    m_CondensedLabels[r] = m_Labels[i];
    r++;
  }

  delete [] marks;

  muLogMacro(
    << "KNNClassifier selected " << m << " / " << n
    << " prototypes for classification\n");

  m_Modified = false;

}

unsigned char
KNNClassifier
::Classify(const VectorType& x)
{
  if (m_Modified)
    this->Update();

  return this->_ClassifyFunction(m_CondensedPrototypes, m_CondensedLabels, x);
}

unsigned char
KNNClassifier
::ClassifyWithoutCondense(const VectorType& x)
{
  if (m_Modified)
    this->Update();

  return this->_ClassifyFunction(m_Prototypes, m_Labels, x);
}

unsigned char
KNNClassifier
::_ClassifyFunction(const MatrixType& set, const LabelArray& labels,
  const VectorType& x)
{

  unsigned int n = set.rows();
  unsigned int d = m_Dimension;

  if (x.size() != d)
    muExceptionMacro(
      << "[KNNClassifier::_ClassifyFunction] Invalid dimension");

  Heap<TaggedDistance> dheap;
  dheap.Allocate(n);

  for (unsigned int i = 0; i < n; i++)
  {
    VectorType v = set.get_row(i) - x;

    TaggedDistance d;
    d.dist = v.squared_magnitude();
    d.label = labels[i];

    dheap.Insert(d);
  }

  unsigned int* counts = new unsigned int[m_LabelRange];
  for (unsigned int i = 0; i < m_LabelRange; i++)
    counts[i] = 0;

  // Find k shortest distances
  for (unsigned int k = 0; k < m_KNeighbors; k++)
  {
    TaggedDistance d = dheap.ExtractMinimum();
    counts[d.label]++;
  }

  // Majority vote
  unsigned int imax = 0;
  unsigned int cmax = counts[0];
  for (unsigned int i = 0; i < m_LabelRange; i++)
  {
    if (counts[i] > cmax)
    {
      imax = i;
      cmax = counts[i];
    }
  }

  // Tie? Pick at random...
  DynArray<unsigned int> tieIndices;
  tieIndices.Allocate(m_LabelRange);
  tieIndices.Append(imax);

  for (unsigned int i = 0; i < m_LabelRange; i++)
  {
    if ((i != imax) && (counts[i] == cmax))
    {
      tieIndices.Append(i);
      //float r = (float)rand() / (float)RAND_MAX;
      //if (r > 0.5)
      //  imax = i;
    }
  }

  if (tieIndices.GetSize() > 1)
  {
    MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();
    unsigned int whichIndex =
      (unsigned int)
        rng->GenerateUniformIntegerUpToK(tieIndices.GetSize() - 1);
    imax = tieIndices[whichIndex];
  }

  delete [] counts;

  return imax;

}
